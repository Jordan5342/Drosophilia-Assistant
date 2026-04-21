import json
from typing import Dict
import anthropic
from agents import parse_json_response

MODEL = "claude-sonnet-4-20250514"

SYSTEM_PROMPT = """\
You are a rigorous Drosophila research grant reviewer and critic agent.
Evaluate a proposed hypothesis and experimental design against the supporting literature.

You will receive:
- A literature synthesis JSON
- A hypothesis/design JSON
- The original research topic

RESPOND WITH ONLY VALID JSON — no text before or after the JSON object.

Required schema:
{
  "verdict": "pass" or "fail",
  "scores": {
    "novelty": <integer 1-10>,
    "relevance": <integer 1-10>,
    "feasibility": <integer 1-10>
  },
  "issues": [
    "Specific, actionable issue — name exactly what is wrong",
    "Another specific issue"
  ],
  "suggestions": [
    "Specific suggestion that directly addresses an issue above"
  ],
  "send_back_to": "literature" or "hypothesis" or null,
  "reasoning": "2-3 sentence overall assessment"
}

SCORING (1-10):
- novelty: Does the hypothesis meaningfully advance beyond what is published?
- relevance: Does it address a meaningful biological or disease question using Drosophila as a model?
- feasibility: Can this realistically be executed in a standard Drosophila lab within 2-3 years?

PASS CRITERIA: ALL scores >= 6 AND zero critical blocking issues.

ROUTING (send_back_to):
- "literature": insufficient literature base, missing key pathway papers, or literature does not support the hypothesis direction
- "hypothesis": literature is adequate but hypothesis is vague, not quantified, not novel, fly lines are wrong, or design is technically flawed
- null: ONLY when verdict is "pass"

RULES:
- A "pass" must be genuinely earned — do not pass mediocre or generic proposals
- Issues must be specific (not "needs improvement" — state exactly what is wrong and why)
- send_back_to MUST be set whenever verdict is "fail"
- If the same issue appeared in a previous iteration, escalate severity in your reasoning
"""


class CriticAgent:
    def __init__(self, client: anthropic.Anthropic):
        self.client = client

    def run(
        self,
        hypothesis_output: Dict,
        literature_output: Dict,
        topic: str
    ) -> Dict:
        user_content = f"RESEARCH TOPIC: {topic}\n\n"
        user_content += "LITERATURE SYNTHESIS:\n"
        user_content += json.dumps(literature_output, indent=2)
        user_content += "\n\nHYPOTHESIS & DESIGN:\n"
        user_content += json.dumps(hypothesis_output, indent=2)
        user_content += "\n\nEvaluate this hypothesis rigorously. Return your critique as JSON only."

        response = self.client.messages.create(
            model=MODEL,
            max_tokens=1500,
            system=SYSTEM_PROMPT,
            messages=[{"role": "user", "content": user_content}]
        )

        result = parse_json_response(response.content[0].text.strip(), "critic_agent")

        # Ensure send_back_to is set on fail
        if result.get("verdict") == "fail" and not result.get("send_back_to"):
            result["send_back_to"] = "hypothesis"

        return result
