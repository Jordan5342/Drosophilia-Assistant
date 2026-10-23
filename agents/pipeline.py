import os
from typing import Dict, Generator, List, Optional
from agents import CLAUDE_MODEL
from agents.literature_agent import LiteratureAgent
from agents.hypothesis_agent import HypothesisAgent
from agents.gemini_critic_agent import GeminiCriticAgent

MAX_ITERATIONS = 1


def _verified_key_papers(reported_papers: List[Dict], retrieved_papers: List[Dict]) -> List[Dict]:
    """Keep only cited papers that match metadata from the literature retrieval."""
    by_pmid = {
        str(paper["pmid"]): paper
        for paper in retrieved_papers
        if paper.get("pmid")
    }
    by_title = {
        " ".join(str(paper.get("title", "")).casefold().split()): paper
        for paper in retrieved_papers
        if paper.get("title")
    }

    verified = []
    for reported in reported_papers or []:
        pmid = str(reported.get("pmid", ""))
        title = " ".join(str(reported.get("title", "")).casefold().split())
        source = by_pmid.get(pmid) if pmid else None
        source = source or (by_title.get(title) if title else None)
        if not source:
            continue

        citation = dict(source)
        if reported.get("relevance"):
            citation["relevance"] = reported["relevance"]
        verified.append(citation)
    return verified


def build_pipeline_proposal(topic: str, literature: Dict, hypothesis: Dict, critic: Dict) -> Dict:
    """Adapt one-pass agent output to the application's persisted proposal schema."""
    literature = literature or {}
    hypothesis = hypothesis or {}
    critic = critic or {}
    aims = []
    for index, aim in enumerate(hypothesis.get("specific_aims") or [], 1):
        if isinstance(aim, dict):
            aim_data = dict(aim)
            aim_data.setdefault("aim_number", index)
            aim_data.setdefault("title", f"Aim {index}")
            aim_data.setdefault("objective", aim_data.get("title", ""))
        else:
            aim_text = str(aim).strip()
            title, separator, objective = aim_text.partition("—")
            aims.append({
                "aim_number": index,
                "title": title.strip() or f"Aim {index}",
                "objective": (objective if separator else aim_text).strip(),
                "approach": hypothesis.get("experimental_approach", ""),
                "expected_outcomes": hypothesis.get("expected_outcomes", ""),
            })
            continue
        aims.append(aim_data)

    papers = literature.get("key_papers") or []
    references = []
    for number, paper in enumerate(papers, 1):
        authors = paper.get("authors") or "Unknown"
        year = paper.get("year") or "N/A"
        title = paper.get("title") or "No title"
        journal = paper.get("journal")
        citation = f"{authors} ({year}). {title}."
        if journal:
            citation += f" {journal}."
        if paper.get("pmid"):
            citation += f" PMID: {paper['pmid']}."
        references.append({
            "number": number,
            "citation": citation,
            "url": paper.get("url") or "",
            "pmid": paper.get("pmid"),
            "doi": paper.get("doi"),
            "source": paper.get("source"),
            "relevance_score": paper.get("_relevance_score"),
            "relevance": paper.get("relevance", ""),
        })

    return {
        "title": f"Research Proposal: {topic}",
        "background": literature.get("summary", ""),
        "central_hypothesis": hypothesis.get("hypothesis", ""),
        "null_hypothesis": "",
        "rationale": hypothesis.get("novelty_claim", ""),
        "specific_aims": aims,
        "experimental_approach": hypothesis.get("experimental_approach", ""),
        "controls": "",
        "timeline": [],
        "fly_lines_needed": hypothesis.get("fly_lines_needed") or [],
        "literature_findings": {
            key: literature[key]
            for key in ("genes", "pathways", "human_orthologs", "disease_context")
            if literature.get(key)
        },
        "literature_gaps": literature.get("identified_gaps") or [],
        "evidence_claims": [
            claim
            for output in (literature, hypothesis)
            for claim in (output.get("evidence_claims") or [])
            if isinstance(claim, dict)
        ],
        "references": references,
        "evidence_status": (
            "retrieved_literature" if references else "no_retrieved_literature"
        ),
        "critic_review": critic,
        "pipeline_status": "partial" if any(
            result.get("_error") for result in (literature, hypothesis, critic)
        ) else "complete",
    }


class AgentPipeline:
    def __init__(self, client, fetch_literature_fn, model: str = CLAUDE_MODEL):
        """
        client: anthropic.Anthropic instance
        fetch_literature_fn: callable matching DrosophilaAssistant.fetch_literature_for_proposal
                             signature: (topic, genes, user_clarifications) -> List[Dict]
        """
        self.lit_agent = LiteratureAgent(client, model)
        self.hyp_agent = HypothesisAgent(client, model)
        self.critic = GeminiCriticAgent(api_key=os.environ.get("GEMINI_API_KEY"))
        self.fetch_literature = fetch_literature_fn

    def run_streaming(
        self,
        query: str,
        topic: str,
        genes: List[str],
        debug: bool = False
    ) -> Generator[Dict, None, None]:
        """
        Generator yielding pipeline events as dicts: {"type": str, "data": dict}.
        Callers format these as SSE: f"event: {type}\\ndata: {json.dumps(data)}\\n\\n"
        """
        history: List[Dict] = []
        best_result: Optional[Dict] = None
        best_score = -1
        cached_papers: Optional[List[Dict]] = None
        cached_lit_output: Optional[Dict] = None

        for iteration in range(1, MAX_ITERATIONS + 1):
            critic_feedback = history[-1]["critic"] if history else None
            send_back = critic_feedback.get("send_back_to") if critic_feedback else None
            run_literature = (iteration == 1) or (send_back == "literature")

            # ── Literature phase ──────────────────────────────────────────────
            if run_literature:
                yield {"type": "agent_start", "data": {
                    "agent": "literature",
                    "iteration": iteration,
                    "label": f"Literature Agent — fetching & synthesizing papers (pass {iteration})"
                }}
                print(f"  [Pipeline pass {iteration}] Literature Agent: fetch + synthesize")
                try:
                    cached_papers = self.fetch_literature(
                        topic=topic,
                        genes=genes,
                        user_clarifications=query
                    )
                    lit_output = self.lit_agent.run(topic, cached_papers, critic_feedback)
                    lit_output["key_papers"] = _verified_key_papers(
                        lit_output.get("key_papers", []),
                        cached_papers or []
                    )
                    cached_lit_output = lit_output
                except Exception as exc:
                    print(f"  [Pipeline pass {iteration}] Literature Agent error: {exc}")
                    lit_output = {
                        "_error": str(exc),
                        "genes": genes,
                        "pathways": [],
                        "human_orthologs": [],
                        "disease_context": "",
                        "key_papers": [],
                        "identified_gaps": ["The literature step failed"],
                        "summary": f"Error: {exc}"
                    }
                    cached_lit_output = lit_output

                yield {"type": "agent_complete", "data": {
                    "agent": "literature",
                    "iteration": iteration,
                    "status": "error" if lit_output.get("_error") else "complete",
                    "error": lit_output.get("_error"),
                    "output": lit_output if debug else _lit_summary(lit_output),
                    "paper_count": len(cached_papers) if cached_papers else 0
                }}
            else:
                # Critic routed back to hypothesis only — reuse cached literature
                lit_output = cached_lit_output
                print(f"  [Pipeline pass {iteration}] Reusing cached literature")
                if debug:
                    yield {"type": "agent_skipped", "data": {
                        "agent": "literature",
                        "iteration": iteration,
                        "reason": "Critic routed back to hypothesis only — literature reused"
                    }}

            # ── Hypothesis phase ──────────────────────────────────────────────
            yield {"type": "agent_start", "data": {
                "agent": "hypothesis",
                "iteration": iteration,
                "label": f"Hypothesis Agent — designing experiment (pass {iteration})"
            }}
            print(f"  [Pipeline pass {iteration}] Hypothesis Agent")
            try:
                hyp_output = self.hyp_agent.run(lit_output, topic, critic_feedback)
            except Exception as exc:
                print(f"  [Pipeline pass {iteration}] Hypothesis Agent error: {exc}")
                hyp_output = {
                    "_error": str(exc),
                    "hypothesis": f"Error: {exc}",
                    "experimental_approach": "",
                    "fly_lines_needed": [],
                    "expected_outcomes": "",
                    "novelty_claim": "",
                    "feasibility_notes": "",
                    "specific_aims": []
                }

            yield {"type": "agent_complete", "data": {
                "agent": "hypothesis",
                "iteration": iteration,
                "status": "error" if hyp_output.get("_error") else "complete",
                "error": hyp_output.get("_error"),
                "output": hyp_output if debug else _hyp_summary(hyp_output)
            }}

            # ── Critic phase ──────────────────────────────────────────────────
            yield {"type": "agent_start", "data": {
                "agent": "critic",
                "iteration": iteration,
                "label": f"Critic Agent — evaluating proposal (pass {iteration})"
            }}
            print(f"  [Pipeline pass {iteration}] Critic Agent")
            try:
                critic_output = self.critic.run(hyp_output, lit_output, topic)
            except Exception as exc:
                print(f"  [Pipeline pass {iteration}] Critic Agent error: {exc}")
                critic_output = {
                    "_error": str(exc),
                    "scores": {},
                    "issues": [f"Critic error: {exc}"],
                    "suggestions": [],
                    "reasoning": f"Critic agent encountered an error: {exc}"
                }

            yield {"type": "agent_complete", "data": {
                "agent": "critic",
                "iteration": iteration,
                "status": "error" if critic_output.get("_error") else "complete",
                "error": critic_output.get("_error"),
                "output": critic_output,
                "verdict": critic_output.get("verdict"),
                "scores": critic_output.get("scores", {})
            }}

            # ── Track best result by total score ──────────────────────────────
            scores = critic_output.get("scores", {})
            total_score = sum(scores.values()) if scores else 0

            history.append({
                "iteration": iteration,
                "reran_literature": run_literature,
                "literature": lit_output,
                "hypothesis": hyp_output,
                "critic": critic_output,
                "total_score": total_score
            })

            if total_score > best_score:
                best_score = total_score
                best_result = {
                    "literature": cached_lit_output,
                    "hypothesis": hyp_output,
                    "critic": critic_output,
                    "score": total_score
                }

            verdict = critic_output.get("verdict")
            send_back_next = critic_output.get("send_back_to")

            # ── Routing decision ──────────────────────────────────────────────
            if verdict == "pass":
                yield {"type": "routing", "data": {
                    "iteration": iteration,
                    "verdict": "pass",
                    "send_back_to": None,
                    "message": "Proposal approved — pipeline complete"
                }}
                break

            if iteration < MAX_ITERATIONS:
                yield {"type": "routing", "data": {
                    "iteration": iteration,
                    "verdict": "fail",
                    "send_back_to": send_back_next,
                    "message": (
                        f"Routing back to {send_back_next or 'literature'} agent "
                        f"(scores: {scores})"
                    )
                }}
            else:
                yield {"type": "routing", "data": {
                    "iteration": iteration,
                    "verdict": verdict,
                    "send_back_to": None,
                    "message": _review_message(verdict, best_score)
                }}

        # ── Final result ──────────────────────────────────────────────────────
        final = best_result or {
            "literature": cached_lit_output or {},
            "hypothesis": {},
            "critic": {"scores": {}, "issues": [], "suggestions": [], "reasoning": ""},
            "score": 0
        }

        failed_agents = [
            agent
            for agent, output in (
                ("literature", final["literature"]),
                ("hypothesis", final["hypothesis"]),
                ("critic", final["critic"])
            )
            if output.get("_error")
        ]
        critic_verdict = final["critic"].get("verdict")
        passed = None
        if critic_verdict in ("pass", "fail"):
            passed = critic_verdict == "pass" and not failed_agents

        yield {"type": "pipeline_complete", "data": {
            "final_literature": final["literature"],
            "final_hypothesis": final["hypothesis"],
            "final_critic": final["critic"],
            "passed": passed,
            "status": "partial" if failed_agents else "complete",
            "failed_agents": failed_agents,
            "best_score": final.get("score", 0),
            "total_iterations": len(history),
            "iterations": history if debug else []
        }}


def _review_message(verdict: Optional[str], score: int) -> str:
    if verdict == "pass":
        return f"Proposal approved — pipeline complete (score {score}/15)"
    if verdict == "fail":
        return f"Single-pass review complete; see critic feedback (score {score}/15)"
    return (
        "Single-pass review complete. Critic scores are descriptive; "
        "no pass/fail verdict was provided."
    )


# ── Non-debug summary helpers ─────────────────────────────────────────────────

def _lit_summary(output: Dict) -> Dict:
    return {
        "genes": output.get("genes", []),
        "paper_count": len(output.get("key_papers", [])),
        "key_papers": output.get("key_papers", []),
        "identified_gaps": output.get("identified_gaps", []),
        "summary": output.get("summary", ""),
        "_error": output.get("_error")
    }


def _hyp_summary(output: Dict) -> Dict:
    return {
        "hypothesis": output.get("hypothesis", ""),
        "novelty_claim": output.get("novelty_claim", ""),
        "aims_count": len(output.get("specific_aims", [])),
        "_error": output.get("_error")
    }
