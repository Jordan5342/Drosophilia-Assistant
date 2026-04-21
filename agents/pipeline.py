from typing import Dict, Generator, List, Optional
from agents.literature_agent import LiteratureAgent
from agents.hypothesis_agent import HypothesisAgent
from agents.critic_agent import CriticAgent

MAX_ITERATIONS = 3


class AgentPipeline:
    def __init__(self, client, fetch_literature_fn):
        """
        client: anthropic.Anthropic instance
        fetch_literature_fn: callable matching DrosophilaAssistant.fetch_literature_for_proposal
                             signature: (topic, genes, user_clarifications) -> List[Dict]
        """
        self.lit_agent = LiteratureAgent(client)
        self.hyp_agent = HypothesisAgent(client)
        self.critic = CriticAgent(client)
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
                    "label": f"Literature Agent — fetching & synthesizing papers (iteration {iteration})"
                }}
                print(f"  [Pipeline i{iteration}] Literature Agent: fetch + synthesize")
                try:
                    cached_papers = self.fetch_literature(
                        topic=topic,
                        genes=genes,
                        user_clarifications=query
                    )
                    lit_output = self.lit_agent.run(topic, cached_papers, critic_feedback)
                    cached_lit_output = lit_output
                except Exception as exc:
                    print(f"  [Pipeline i{iteration}] Literature Agent error: {exc}")
                    lit_output = {
                        "_error": str(exc),
                        "genes": genes,
                        "pathways": [],
                        "human_orthologs": [],
                        "disease_context": "",
                        "key_papers": [],
                        "identified_gaps": ["Literature retrieval failed"],
                        "summary": f"Error: {exc}"
                    }
                    cached_lit_output = lit_output

                yield {"type": "agent_complete", "data": {
                    "agent": "literature",
                    "iteration": iteration,
                    "output": lit_output if debug else _lit_summary(lit_output),
                    "paper_count": len(cached_papers) if cached_papers else 0
                }}
            else:
                # Critic routed back to hypothesis only — reuse cached literature
                lit_output = cached_lit_output
                print(f"  [Pipeline i{iteration}] Reusing cached literature (critic sent back to hypothesis)")
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
                "label": f"Hypothesis Agent — designing experiment (iteration {iteration})"
            }}
            print(f"  [Pipeline i{iteration}] Hypothesis Agent")
            try:
                hyp_output = self.hyp_agent.run(lit_output, topic, critic_feedback)
            except Exception as exc:
                print(f"  [Pipeline i{iteration}] Hypothesis Agent error: {exc}")
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
                "output": hyp_output if debug else _hyp_summary(hyp_output)
            }}

            # ── Critic phase ──────────────────────────────────────────────────
            yield {"type": "agent_start", "data": {
                "agent": "critic",
                "iteration": iteration,
                "label": f"Critic Agent — evaluating proposal (iteration {iteration})"
            }}
            print(f"  [Pipeline i{iteration}] Critic Agent")
            try:
                critic_output = self.critic.run(hyp_output, lit_output, topic)
            except Exception as exc:
                print(f"  [Pipeline i{iteration}] Critic Agent error: {exc}")
                critic_output = {
                    "_error": str(exc),
                    "verdict": "fail",
                    "scores": {"novelty": 5, "relevance": 5, "feasibility": 5},
                    "issues": [f"Critic error: {exc}"],
                    "suggestions": [],
                    "send_back_to": "hypothesis",
                    "reasoning": f"Critic agent encountered an error: {exc}"
                }

            yield {"type": "agent_complete", "data": {
                "agent": "critic",
                "iteration": iteration,
                "output": critic_output,
                "verdict": critic_output.get("verdict", "fail"),
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

            verdict = critic_output.get("verdict", "fail")
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
                    "verdict": "fail",
                    "send_back_to": None,
                    "message": (
                        f"Max iterations reached. Surfacing best attempt "
                        f"(best score {best_score}/30)"
                    )
                }}

        # ── Final result ──────────────────────────────────────────────────────
        final = best_result or {
            "literature": cached_lit_output or {},
            "hypothesis": {},
            "critic": {"verdict": "fail", "scores": {}, "issues": [], "suggestions": [], "reasoning": ""},
            "score": 0
        }

        yield {"type": "pipeline_complete", "data": {
            "final_literature": final["literature"],
            "final_hypothesis": final["hypothesis"],
            "final_critic": final["critic"],
            "passed": final["critic"].get("verdict") == "pass",
            "best_score": final.get("score", 0),
            "total_iterations": len(history),
            "iterations": history if debug else []
        }}


# ── Non-debug summary helpers ─────────────────────────────────────────────────

def _lit_summary(output: Dict) -> Dict:
    return {
        "genes": output.get("genes", []),
        "paper_count": len(output.get("key_papers", [])),
        "identified_gaps": output.get("identified_gaps", []),
        "summary": output.get("summary", "")
    }


def _hyp_summary(output: Dict) -> Dict:
    return {
        "hypothesis": output.get("hypothesis", ""),
        "novelty_claim": output.get("novelty_claim", ""),
        "aims_count": len(output.get("specific_aims", []))
    }
