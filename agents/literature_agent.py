import json
from typing import Dict, List, Optional
import anthropic
from agents import parse_json_response

MODEL = "claude-sonnet-4-20250514"

SYSTEM_PROMPT = """\
You are a Drosophila melanogaster literature synthesis agent.
Analyze provided research papers and return a structured JSON synthesis.

You will receive:
- A research topic/query
- Paper abstracts from PubMed and FlyBase
- (Optional) Critic feedback from a previous iteration requesting specific improvements

RESPOND WITH ONLY VALID JSON — no text before or after the JSON object.

Required schema:
{
  "genes": ["list of Drosophila gene symbols relevant to the topic"],
  "pathways": ["list of relevant signaling pathways"],
  "human_orthologs": [
    {"fly_gene": "foxo", "human_gene": "FOXO3", "disease": "Longevity, Cancer"}
  ],
  "disease_context": "How this topic connects to human disease (2-3 sentences)",
  "key_papers": [
    {
      "pmid": "12345678",
      "title": "Full paper title",
      "authors": "Author et al.",
      "year": "2022",
      "url": "https://pubmed.ncbi.nlm.nih.gov/12345678/",
      "relevance": "One sentence explaining why this paper is critical"
    }
  ],
  "identified_gaps": [
    "Specific gap 1 — what is unknown and why it matters",
    "Specific gap 2"
  ],
  "summary": "2-3 sentence synthesis of the current state of the field"
}

RULES:
- Only include genes and papers from the PROVIDED literature
- Gaps must be specific and actionable (not generic "more research needed")
- If critic feedback is provided, directly address every issue raised
- Human orthologs must include disease relevance
"""


class LiteratureAgent:
    def __init__(self, client: anthropic.Anthropic):
        self.client = client

    def run(
        self,
        topic: str,
        papers: List[Dict],
        critic_feedback: Optional[Dict] = None
    ) -> Dict:
        user_content = f"RESEARCH TOPIC: {topic}\n\n"

        if critic_feedback and not critic_feedback.get("_parse_error"):
            user_content += "CRITIC FEEDBACK FROM PREVIOUS ITERATION (address every issue):\n"
            user_content += json.dumps(critic_feedback, indent=2)
            user_content += "\n\n"

        if papers:
            user_content += f"PROVIDED LITERATURE ({len(papers)} papers):\n\n"
            for i, p in enumerate(papers, 1):
                user_content += f"[{i}] {p.get('authors', 'Unknown')} ({p.get('year', 'N/A')})\n"
                user_content += f"    Title: {p.get('title', 'No title')}\n"
                user_content += f"    PMID: {p.get('pmid', 'N/A')}\n"
                user_content += f"    URL: {p.get('url', '')}\n"
                user_content += f"    Abstract: {p.get('abstract', 'No abstract available')}\n"
                user_content += f"    Source: {p.get('source', 'Unknown')}\n\n"
        else:
            user_content += (
                "No papers were retrieved for this topic. "
                "Synthesize from your knowledge of Drosophila biology but clearly "
                "flag every claim as unverified and maximize the identified_gaps list.\n\n"
            )

        user_content += "Return your synthesis as JSON only."

        response = self.client.messages.create(
            model=MODEL,
            max_tokens=2000,
            system=SYSTEM_PROMPT,
            messages=[{"role": "user", "content": user_content}]
        )

        return parse_json_response(response.content[0].text.strip(), "literature_agent")
