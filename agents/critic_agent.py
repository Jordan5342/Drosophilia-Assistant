import concurrent.futures
import json
from typing import Dict

import anthropic

from agents import parse_json_response

MODEL = "claude-sonnet-4-20250514"

REVIEWER_PROMPTS = {
    "feasibility": """\
You are a Drosophila lab manager reviewing a proposal for technical feasibility.
You care ONLY about whether this can actually be executed.

Red flags you look for:
- Fly lines that don't exist at BDSC or are not publicly available
- Cross schemes that are genetically impossible or take >6 months to establish
- Assays that require equipment a standard fly lab wouldn't have
- Sample sizes that are unrealistic for the timeline stated
- No mention of controls or balancer chromosomes where needed

You do NOT care about novelty or disease relevance. Only: can a grad student actually do this?

RESPOND WITH ONLY VALID JSON:
{
  "score": <integer 1-10>,
  "issues": ["specific feasibility problem 1", "specific feasibility problem 2"],
  "suggestions": ["specific fix 1"],
  "reasoning": "1-2 sentences on feasibility only"
}""",

    "novelty": """\
You are a senior Drosophila researcher who has reviewed hundreds of grant proposals.
You care ONLY about whether this hypothesis advances the field.

Red flags you look for:
- Hypothesis that essentially repeats published work with a minor variation
- "We will characterize X" without a clear mechanistic prediction
- No engagement with why existing results are insufficient
- Novelty claim that doesn't actually reference a specific gap in the provided literature
- Aims that are descriptive rather than mechanistic

You do NOT care about feasibility or disease framing. Only: does this move the science forward?

RESPOND WITH ONLY VALID JSON:
{
  "score": <integer 1-10>,
  "issues": ["specific novelty problem 1"],
  "suggestions": ["specific fix 1"],
  "reasoning": "1-2 sentences on novelty only"
}""",

    "rigor": """\
You are a biostatistician and experimental design reviewer.
You care ONLY about whether this experiment will produce interpretable results.

Red flags you look for:
- No stated sample size or power justification
- Missing negative controls or inappropriate positive controls
- Confounds that aren't addressed (e.g. mating status, vial density, temperature variation)
- Statistical test not specified or wrong test for the data type
- Expected outcomes that aren't actually measurable or are too vague to falsify
- Single-sex experiments without justification

You do NOT care about novelty or broader impact. Only: will the data be interpretable?

RESPOND WITH ONLY VALID JSON:
{
  "score": <integer 1-10>,
  "issues": ["specific rigor problem 1"],
  "suggestions": ["specific fix 1"],
  "reasoning": "1-2 sentences on rigor only"
}""",
}


class CriticAgent:
    def __init__(self, client: anthropic.Anthropic):
        self.client = client

    def run(self, hypothesis_output: Dict, literature_output: Dict, topic: str) -> Dict:
        user_content = f"RESEARCH TOPIC: {topic}\n\n"
        user_content += "LITERATURE SYNTHESIS:\n"
        user_content += json.dumps(literature_output, indent=2)
        user_content += "\n\nHYPOTHESIS & DESIGN:\n"
        user_content += json.dumps(hypothesis_output, indent=2)
        user_content += "\n\nEvaluate this proposal from your specific perspective. Return JSON only."

        def call_reviewer(name_and_prompt):
            name, prompt = name_and_prompt
            response = self.client.messages.create(
                model=MODEL,
                max_tokens=800,
                system=prompt,
                messages=[{"role": "user", "content": user_content}],
            )
            result = parse_json_response(response.content[0].text.strip(), f"critic_{name}")
            result["_reviewer"] = name
            return name, result

        reviews = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            for name, review in executor.map(call_reviewer, REVIEWER_PROMPTS.items()):
                reviews[name] = review

        scores = {
            "novelty":     reviews.get("novelty",     {}).get("score", 5),
            "feasibility": reviews.get("feasibility", {}).get("score", 5),
            "rigor":       reviews.get("rigor",       {}).get("score", 5),
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

        verdict = "pass" if all(v >= 6 for v in scores.values()) else "fail"

        send_back_to = None
        if verdict == "fail":
            send_back_to = "hypothesis"
            novelty_issues = " ".join(reviews.get("novelty", {}).get("issues", []))
            if "literature" in novelty_issues.lower() or "papers" in novelty_issues.lower():
                send_back_to = "literature"

        return {
            "verdict": verdict,
            "scores": scores,
            "issues": all_issues,
            "suggestions": all_suggestions,
            "send_back_to": send_back_to,
            "reasoning": reasoning,
            "reviewer_details": reviews,
        }
