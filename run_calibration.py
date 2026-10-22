"""
Standalone calibration run: generates N fresh hypotheses through the real
literature -> hypothesis -> critic pipeline, scores them with an
independent-model critic (Gemini, not Claude), and saves everything to a
local JSON file for the human grading round.

Deliberately bypasses Flask, Supabase, and the iterate-until-pass routing
loop in agents/pipeline.py. For calibration data we want ONE clean
literature->hypothesis->critic pass per topic (showing natural quality
variation), not the system's self-selected best-of-3 attempt.

Usage:
    python run_calibration.py
    python run_calibration.py --n 10 --out calibration_run.json

Requires a .env file (see .env.example) with:
    ANTHROPIC_API_KEY   - for literature + hypothesis agents
    GEMINI_API_KEY       - for the independent critic
    NCBI_API_KEY          - optional, raises PubMed rate limit (not required)
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

import anthropic

from agents.literature_agent import LiteratureAgent
from agents.hypothesis_agent import HypothesisAgent
from agents.gemini_critic_agent import GeminiCriticAgent
from drosophila_assistant import ANTHROPIC_TIMEOUT_SECONDS, DrosophilaAssistant, _GENE_DB

TOPIC_GENERATION_MODEL = "claude-sonnet-4-6"

TOPIC_GENERATION_PROMPT = """\
You are helping assemble a calibration set for a Drosophila melanogaster cancer/developmental \
biology research-hypothesis pipeline.

Generate {n} distinct, varied research topics suitable for a Drosophila lab. Vary the \
biological domain (signaling pathways, tumor models, developmental processes, aging, \
neurobiology, stem cell biology, etc.) so the set is not repetitive. Each topic should be \
specific enough to drive a real literature search (i.e. not just "cancer" but something like \
"Hippo pathway dysregulation in intestinal stem cell tumorigenesis").

RESPOND WITH ONLY VALID JSON, a list of strings, no other text:
["topic 1", "topic 2", ...]
"""


def generate_topics(client: anthropic.Anthropic, n: int) -> list:
    response = client.messages.create(
        model=TOPIC_GENERATION_MODEL,
        max_tokens=1000,
        messages=[{"role": "user", "content": TOPIC_GENERATION_PROMPT.format(n=n)}],
    )
    raw = response.content[0].text.strip()
    raw = raw.strip("`")
    if raw.startswith("json"):
        raw = raw[4:].strip()
    try:
        topics = json.loads(raw)
        assert isinstance(topics, list)
        return topics[:n]
    except Exception:
        print("  ⚠️  Topic generation returned malformed JSON, raw output:")
        print(raw[:500])
        sys.exit(1)


def load_topics(path: str) -> list:
    """Topics from a previous calibration output ({"results": [{"topic": ...}]}) or a JSON list."""
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict):
        return [r["topic"] for r in data["results"] if r.get("topic")]
    return list(data)


_LEGACY_EXCLUDE = {'aging', 'gene', 'genes', 'protein', 'proteins', 'what', 'is', 'are', 'tell',
                   'me', 'about', 'a', 'n', 'the', 'of', 'in', 'sleep', 'memory', 'behavior'}


def legacy_extract_gene_names(text: str) -> list:
    """Pre-fix substring matcher ("Drosophila" -> dr), kept only to measure the fix's effect."""
    text_lower = text.lower()
    return [k for k in _GENE_DB if len(k) > 1 and k in text_lower and k not in _LEGACY_EXCLUDE][:3]


def run_single(assistant: DrosophilaAssistant, lit_agent, hyp_agent, critic_agent, topic: str,
               extract_genes=None) -> dict:
    print(f"\n{'='*70}\nTOPIC: {topic}\n{'='*70}")

    genes = (extract_genes or assistant.extract_gene_names)(topic)
    print(f"  🧬 Extracted genes: {genes}")

    print("  📚 Literature agent: fetching + synthesizing...")
    try:
        papers = assistant.fetch_literature_for_proposal(topic=topic, genes=genes, user_clarifications="")
    except Exception as exc:
        print(f"  ⚠️  Literature fetch failed: {exc}")
        papers = []

    lit_output = lit_agent.run(topic, papers)
    print(f"  ✅ Literature done. Gaps identified: {len(lit_output.get('identified_gaps', []))}")

    print("  🧪 Hypothesis agent: designing experiment...")
    hyp_output = hyp_agent.run(lit_output, topic)
    print(f"  ✅ Hypothesis: {hyp_output.get('hypothesis', '')[:100]}...")

    print("  🔍 Critic agent (Gemini): evaluating...")
    critic_output = critic_agent.run(hyp_output, lit_output, topic)
    print(f"  ✅ Scores: {critic_output.get('scores', {})}")

    return {
        "topic": topic,
        "genes": genes,
        "paper_count": len(papers),
        "papers": [{"pmid": p.get("pmid"), "title": p.get("title"), "source": p.get("source")} for p in papers],
        "literature": lit_output,
        "hypothesis": hyp_output,
        "critic": critic_output,
    }


def main():
    parser = argparse.ArgumentParser(description="Run a calibration batch of hypotheses through the pipeline.")
    parser.add_argument("--n", type=int, default=10, help="Number of hypotheses to generate (default 10)")
    parser.add_argument("--out", type=str, default="calibration_run.json", help="Output JSON file path")
    parser.add_argument("--topics-from", type=str, help="Reuse topics from a previous calibration JSON (or a JSON list) instead of generating new ones")
    parser.add_argument("--gene-extraction", choices=["fixed", "legacy"], default="fixed",
                        help="'legacy' reproduces the pre-fix substring gene matcher for before/after comparison")
    args = parser.parse_args()

    anthropic_key = os.environ.get("ANTHROPIC_API_KEY")
    gemini_key = os.environ.get("GEMINI_API_KEY")

    if not anthropic_key:
        print("❌ ANTHROPIC_API_KEY not set. Add it to your .env file.")
        sys.exit(1)
    if not gemini_key:
        print("❌ GEMINI_API_KEY not set. Add it to your .env file.")
        sys.exit(1)

    client = anthropic.Anthropic(api_key=anthropic_key, timeout=ANTHROPIC_TIMEOUT_SECONDS)
    assistant = DrosophilaAssistant(api_key=anthropic_key)
    lit_agent = LiteratureAgent(client)
    hyp_agent = HypothesisAgent(client)
    critic_agent = GeminiCriticAgent(api_key=gemini_key)

    if args.topics_from:
        topics = load_topics(args.topics_from)[:args.n]
        print(f"🧠 Reusing {len(topics)} topics from {args.topics_from}")
    else:
        print(f"🧠 Generating {args.n} varied research topics...")
        topics = generate_topics(client, args.n)
    extract_genes = legacy_extract_gene_names if args.gene_extraction == "legacy" else None
    print(f"🧬 Gene extraction: {args.gene_extraction}")
    for i, t in enumerate(topics, 1):
        print(f"  {i}. {t}")

    results = []
    for i, topic in enumerate(topics, 1):
        print(f"\n\n### Hypothesis {i}/{len(topics)} ###")
        try:
            result = run_single(assistant, lit_agent, hyp_agent, critic_agent, topic, extract_genes)
            result["index"] = i
            results.append(result)
        except Exception as exc:
            print(f"  ❌ Failed on topic '{topic}': {exc}")
            results.append({"index": i, "topic": topic, "_error": str(exc)})
        time.sleep(1)  # light pacing, stay well under any rate limits

    out_path = Path(args.out)
    with open(out_path, "w") as f:
        json.dump({"generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "gene_extraction": args.gene_extraction, "results": results}, f, indent=2)

    print(f"\n\n✅ Done. Wrote {len(results)} entries to {out_path.resolve()}")


if __name__ == "__main__":
    main()
