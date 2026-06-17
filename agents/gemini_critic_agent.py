"""
Gemini-based critic agent for the Drosophila hypothesis pipeline.

Why this exists: the original CriticAgent (agents/critic_agent.py) uses the
same model (Claude) that the HypothesisAgent uses to generate proposals.
That creates a circularity problem for ground-truth calibration -- the
critic may share blind spots with the generator, and there's no independent
signal to check the rubric itself against. Using a different model family
(Gemini) as critic removes that shared-model confound, on top of the human
grading layer added separately.

Rubric for this calibration round (replaces feasibility/rigor/gap_specificity,
1-10, pass/fail):
    - hallucination_defense : does the hypothesis avoid asserting fabricated
      genes/papers/facts, AND does it correctly separate established fact
      from speculation (hedges appropriately, doesn't overclaim)?
    - genetic_feasibility   : are the fly lines, crosses, and genetic tools
      real, available, and workable in a standard timeline?
    - control_reasoning     : are appropriate controls (genetic, statistical,
      environmental) specified and justified?

Scale: integer 1-5 for each dimension. No pass/fail verdict -- this round is
purely descriptive calibration data for human graders, not a routing gate.
"""

import concurrent.futures
import json
import os
from typing import Dict

from google import genai
from google.genai import types

from agents import parse_json_response

# Gemini 3 Flash is the current recommended free-tier model as of June 2026.
# If your account/region resolves a different string, override via env var.
MODEL = os.environ.get("GEMINI_CRITIC_MODEL", "gemini-3-flash")

REVIEWER_PROMPTS = {
    "hallucination_defense": """\
You are a Drosophila geneticist fact-checking a research proposal for fabricated or \
misrepresented information.
You care ONLY about whether every factual claim in this proposal is real and \
appropriately hedged.

Red flags you look for:
- Gene names, FBgn IDs, or BDSC stock numbers that do not actually exist
- Citations to papers, findings, or pathway interactions that are invented or misattributed
- Claims stated as established fact when they are actually speculative or unverified
- A "novelty" or "gap" claim that is factually wrong because the topic is already well-published
- Human ortholog or disease-relevance claims that overstate what is actually known

Green flags:
- Specific claims are traceable to real, named prior work or established pathway biology
- Uncertain or speculative claims are explicitly flagged as such, not stated as fact
- The proposal stays within what is plausible given standard Drosophila biology

You do NOT care about feasibility of execution or control design. Only: is this true, \
and is uncertainty handled honestly?

RESPOND WITH ONLY VALID JSON:
{
  "score": <integer 1-5>,
  "issues": ["specific fabrication or overclaiming problem 1"],
  "suggestions": ["specific fix 1"],
  "reasoning": "1-2 sentences on hallucination defense only"
}""",

    "genetic_feasibility": """\
You are a Drosophila lab manager reviewing a proposal for genetic and technical feasibility.
You care ONLY about whether the genetic tools and experimental plan can actually be executed.

Red flags you look for:
- Fly lines that don't exist at BDSC or are not realistically obtainable
- Cross schemes that are genetically impossible, miscount generations, or take far longer \
than the timeline implies
- GAL4/UAS, balancer, or marker logic that doesn't make genetic sense
- Assays requiring equipment a standard fly lab wouldn't have
- Sample sizes unrealistic for the stated timeline

You do NOT care about whether claims are factually true (that's a separate reviewer's job) \
or about control design specifically. Only: can a grad student actually execute this cross \
and assay as described?

RESPOND WITH ONLY VALID JSON:
{
  "score": <integer 1-5>,
  "issues": ["specific feasibility problem 1"],
  "suggestions": ["specific fix 1"],
  "reasoning": "1-2 sentences on genetic feasibility only"
}""",

    "control_reasoning": """\
You are a biostatistician and experimental design reviewer focused entirely on controls.
You care ONLY about whether the proposed controls are sufficient to make the result \
interpretable.

Red flags you look for:
- Missing negative controls or an inappropriate positive control
- No genetic background control (e.g. GAL4-only or UAS-only controls absent)
- Confounds not addressed (mating status, vial density, temperature, time of day)
- Single-sex experiments without justification
- No mention of how results will be distinguished from off-target or genetic background effects

Green flags:
- Explicit genetic controls (driver-alone, responder-alone) are named
- Confounds are anticipated and addressed
- The logic connecting "what we measure" to "what we conclude" is sound

You do NOT care about whether the science itself is novel or true. Only: if this experiment \
runs exactly as described, will the controls let you trust the result?

RESPOND WITH ONLY VALID JSON:
{
  "score": <integer 1-5>,
  "issues": ["specific control problem 1"],
  "suggestions": ["specific fix 1"],
  "reasoning": "1-2 sentences on control reasoning only"
}""",
}


class GeminiCriticAgent:
    def __init__(self, api_key: str = None):
        self.client = genai.Client(api_key=api_key or os.environ.get("GEMINI_API_KEY"))

    def run(self, hypothesis_output: Dict, literature_output: Dict, topic: str) -> Dict:
        user_content = f"RESEARCH TOPIC: {topic}\n\n"
        user_content += "LITERATURE SYNTHESIS:\n"
        user_content += json.dumps(literature_output, indent=2)
        user_content += "\n\nHYPOTHESIS & DESIGN:\n"
        user_content += json.dumps(hypothesis_output, indent=2)
        user_content += "\n\nEvaluate this proposal from your specific perspective. Return JSON only."

        def call_reviewer(name_and_prompt):
            name, prompt = name_and_prompt
            response = self.client.models.generate_content(
                model=MODEL,
                contents=user_content,
                config=types.GenerateContentConfig(
                    system_instruction=prompt,
                    max_output_tokens=800,
                    temperature=0.2,
                ),
            )
            raw = (response.text or "").strip()
            result = parse_json_response(raw, f"gemini_critic_{name}")
            result["_reviewer"] = name
            return name, result

        reviews = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            for name, review in executor.map(call_reviewer, REVIEWER_PROMPTS.items()):
                reviews[name] = review

        scores = {
            "hallucination_defense": reviews.get("hallucination_defense", {}).get("score", 3),
            "genetic_feasibility": reviews.get("genetic_feasibility", {}).get("score", 3),
            "control_reasoning": reviews.get("control_reasoning", {}).get("score", 3),
        }

        all_issues = []
        all_suggestions = []
        for review in reviews.values():
            all_issues.extend(review.get("issues", []))
            all_suggestions.extend(review.get("suggestions", []))

        reasoning = " | ".join(
            f"{name.upper()}: {r.get('reasoning', '')}"
            for name, r in reviews.items()
            if r.get("reasoning")
        )

        # No pass/fail verdict for this calibration round -- scores are
        # descriptive output for human graders, not a routing gate.
        return {
            "scores": scores,
            "issues": all_issues,
            "suggestions": all_suggestions,
            "reasoning": reasoning,
            "reviewer_details": reviews,
            "critic_model": MODEL,
        }