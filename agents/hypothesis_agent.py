import json
from typing import Dict, List, Optional
import anthropic
from agents import parse_json_response

MODEL = "claude-sonnet-4-20250514"

SYSTEM_PROMPT = """\
You are a Drosophila melanogaster experimental design agent.
Design a specific, testable hypothesis and experimental plan based on a literature synthesis.

You will receive:
- A literature synthesis JSON from the Literature Agent
- The original research topic
- (Optional) Critic feedback requesting specific improvements — you MUST address every issue

RESPOND WITH ONLY VALID JSON — no text before or after the JSON object.

Required schema:
{
  "hypothesis": "A specific, falsifiable, quantified hypothesis (e.g., 'Overexpression of foxo in the fat body will extend median lifespan by ≥20% in mated female w1118 flies under 25°C standard conditions')",
  "experimental_approach": "Detailed experimental strategy: GAL4/UAS system details, tissue specificity, assay types, controls (3-5 sentences)",
  "fly_lines_needed": [
    {
      "line": "Specific genotype in Drosophila notation (e.g., w*; UAS-foxo/+)",
      "purpose": "Role in the experiment",
      "source": "BDSC stock number or lab source"
    }
  ],
  "expected_outcomes": "Specific, quantified expected results with measurable endpoints and statistical thresholds",
  "novelty_claim": "What makes this hypothesis novel — cite specific gaps from the literature synthesis",
  "feasibility_notes": "Technical feasibility, estimated timeline (months), potential pitfalls",
  "specific_aims": [
    "Aim 1: Short title — specific objective",
    "Aim 2: Short title — specific objective",
    "Aim 3: Short title — specific objective"
  ]
}

RULES:
- Hypothesis MUST be falsifiable and quantified (include numbers and conditions)
- Fly lines must use proper Drosophila genotype notation
- Novelty must explicitly reference gaps identified in the literature synthesis
- If critic feedback is provided, address EVERY issue before anything else
- Specific aims must be discrete, independently testable units
"""


class HypothesisAgent:
    def __init__(self, client: anthropic.Anthropic):
        self.client = client

    def run(
        self,
        literature_output: Dict,
        topic: str,
        critic_feedback: Optional[Dict] = None
    ) -> Dict:
        user_content = f"RESEARCH TOPIC: {topic}\n\n"
        user_content += "LITERATURE SYNTHESIS:\n"
        user_content += json.dumps(literature_output, indent=2)
        user_content += "\n\n"

        if critic_feedback and not critic_feedback.get("_parse_error"):
            user_content += "CRITIC FEEDBACK — ADDRESS EVERY ISSUE LISTED BELOW:\n"
            issues = critic_feedback.get("issues", [])
            suggestions = critic_feedback.get("suggestions", [])
            user_content += f"Issues:\n" + "\n".join(f"  - {i}" for i in issues) + "\n"
            if suggestions:
                user_content += f"Suggestions:\n" + "\n".join(f"  - {s}" for s in suggestions) + "\n"
            user_content += f"Routing reason: {critic_feedback.get('reasoning', '')}\n\n"

        user_content += "Design a specific, testable hypothesis and experimental plan. Return JSON only."

        response = self.client.messages.create(
            model=MODEL,
            max_tokens=2000,
            system=SYSTEM_PROMPT,
            messages=[{"role": "user", "content": user_content}]
        )

        return parse_json_response(response.content[0].text.strip(), "hypothesis_agent")
