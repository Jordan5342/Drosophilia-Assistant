import anthropic
import re
from typing import Optional, Dict, List, Tuple

# Intent keywords that trigger research planning mode (conservative)
PLANNING_INTENT_STRONG = [
    'write a proposal', 'create a proposal', 'build a proposal', 'draft a proposal',
    'design an experiment', 'create an experiment', 'plan an experiment',
    'generate a hypothesis', 'create a hypothesis', 'formulate a hypothesis',
    'help me research', 'research plan', 'experimental design',
    'specific aims', 'research project', 'study design',
    'help me test', 'how would i test', 'how would i study',
    'plan research', 'create a study'
]

PLANNING_INTENT_CONTEXTUAL = [
    'expand on this', 'build on this', 'go further', 'next steps',
    'how to test this', 'experiment around', 'proposal around',
    'hypothesis around', 'study this further'
]

PROPOSAL_SECTIONS = [
    'title',
    'background',
    'central_hypothesis',
    'specific_aims',
    'experimental_approach',
    'potential_pitfalls',
    'timeline',
    'references'
]


def detect_planning_intent(message: str, has_prior_context: bool = False) -> Tuple[bool, str]:
    """
    Detect if message has research planning intent.
    Returns (is_planning, confidence) where confidence is 'strong' or 'contextual'.
    Conservative: only fires on clear signals.
    """
    msg_lower = message.lower()

    for phrase in PLANNING_INTENT_STRONG:
        if phrase in msg_lower:
            return True, 'strong'

    if has_prior_context:
        for phrase in PLANNING_INTENT_CONTEXTUAL:
            if phrase in msg_lower:
                return True, 'contextual'

    return False, ''


class ResearchPlanner:
    def __init__(self, anthropic_client: anthropic.Anthropic):
        self.client = anthropic_client
        self.current_proposal: Optional[Dict] = None
        self.proposal_context: Dict = {}  # topic, genes, papers from prior chat

    def set_context_from_chat(self, topic: str, genes: List[str], papers: List[Dict]):
        """Called by the assistant to seed the planner with chat context."""
        self.proposal_context = {
            'topic': topic,
            'genes': genes,
            'papers': papers
        }

    def extract_topic_from_message(self, message: str) -> str:
        """Pull the core research topic out of the user's message."""
        # Strip planning keywords to get the topic
        clean = message.lower()
        for phrase in PLANNING_INTENT_STRONG + PLANNING_INTENT_CONTEXTUAL:
            clean = clean.replace(phrase, '')

        # Remove filler words
        fillers = ['for', 'about', 'on', 'regarding', 'around', 'this', 'the', 'a', 'an', 'me', 'please', 'can you', 'could you', 'i want to', 'i would like to']
        for filler in fillers:
            clean = re.sub(r'\b' + filler + r'\b', '', clean)

        topic = clean.strip().strip(',').strip()
        return topic if topic else message.strip()

    def generate_clarifying_questions(self, topic: str, genes: List[str]) -> str:
        """
        Ask 1-2 targeted questions before generating the proposal.
        Only fires when topic is ambiguous or critical info is missing.
        """
        gene_str = ', '.join(genes) if genes else 'the gene(s) of interest'

        prompt = f"""A researcher wants to build a Drosophila research proposal about: "{topic}"
Identified genes: {gene_str}

Generate 1-2 SHORT clarifying questions (not more) that would meaningfully improve the proposal.
Focus on: (1) specific phenotype or biological question they care about, (2) available techniques or lab resources if unclear.
Do NOT ask if the info is already obvious from the topic.
Format as a brief friendly message, not a list. Keep it under 60 words total."""

        response = self.client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=200,
            messages=[{"role": "user", "content": prompt}]
        )
        return response.content[0].text

    def generate_proposal(
        self,
        topic: str,
        genes: List[str],
        papers: List[Dict],
        user_clarifications: str = "",
        conversation_history: List[Dict] = None
    ) -> Dict:
        """
        Generate a full structured research proposal.
        Returns a dict with all sections for rendering and docx export.
        """

        # Format papers for the prompt
        paper_context = ""
        if papers:
            paper_context = "\n\nAVAILABLE LITERATURE (use these as real citations):\n"
            for i, p in enumerate(papers[:10], 1):
                paper_context += f"[{i}] {p.get('authors', 'Unknown')} ({p.get('year', 'N/A')}). {p.get('title', 'No title')}. PMID: {p.get('pmid', 'N/A')}. URL: {p.get('url', '')}\n"
                if p.get('abstract'):
                    paper_context += f"    Abstract: {p['abstract'][:300]}...\n"

        gene_str = ', '.join(genes) if genes else 'the relevant gene(s)'
        clarification_str = f"\nAdditional researcher context: {user_clarifications}" if user_clarifications else ""

        # Build conversation context summary if available
        chat_context = ""
        if conversation_history:
            recent = conversation_history[-6:]  # last 3 exchanges
            chat_context = "\nPRIOR CONVERSATION CONTEXT:\n"
            for msg in recent:
                role = "Researcher" if msg['role'] == 'user' else "Assistant"
                content = msg['content']
                if '======' in content:
                    content = content.split('\n\n')[0]
                chat_context += f"{role}: {content[:300]}\n"

        system_prompt = """You are an expert Drosophila researcher and scientific writing specialist helping grad students and postdocs build rigorous research proposals.

Generate a complete, scientifically rigorous generic academic research proposal in JSON format.

CRITICAL RULES:
1. Cite the provided papers using [1], [2], etc. notation throughout
2. If a section has thin literature, explicitly flag it as: "⚠️ LITERATURE GAP: [description]" — these are research opportunities
3. Be specific about Drosophila techniques (GAL4/UAS, CRISPR, behavioral assays, genetic screens, etc.)
4. The hypothesis must be testable and falsifiable
5. Pitfalls must be specific and include mitigation strategies
6. Timeline must be realistic for academic research

Return ONLY valid JSON with this exact structure:
{
  "title": "Full descriptive proposal title",
  "background": "2-3 paragraphs covering what is known, key findings, and why this matters. Cite papers heavily. Flag any gaps.",
  "central_hypothesis": "One clear, testable hypothesis statement",
  "null_hypothesis": "The corresponding null hypothesis",
  "rationale": "Why this hypothesis is worth testing (2-3 sentences)",
  "specific_aims": [
    {
      "aim_number": 1,
      "title": "Short aim title",
      "objective": "What this aim will accomplish",
      "approach": "Specific experimental methods",
      "expected_outcomes": "What you expect to find and why",
      "potential_pitfalls": "What could go wrong",
      "mitigation": "How to address those pitfalls"
    }
  ],
  "experimental_approach": "Overview of the overall experimental strategy, including key Drosophila strains, tools, and readouts",
  "controls": "Key positive and negative controls",
  "timeline": [
    {"phase": "Phase 1 (Months 1-4)", "milestones": "Key milestones for this phase"},
    {"phase": "Phase 2 (Months 5-10)", "milestones": "Key milestones"},
    {"phase": "Phase 3 (Months 11-16)", "milestones": "Key milestones"}
  ],
  "literature_gaps": ["Gap 1 description", "Gap 2 description"],
  "references": [
    {"number": 1, "citation": "Author et al. (Year). Title. Journal. PMID: XXXXX", "url": "https://..."}
  ]
}"""

        user_prompt = f"""Create a research proposal for the following:

Topic: {topic}
Key genes/proteins: {gene_str}{clarification_str}
{paper_context}
{chat_context}

Generate 2-3 specific aims. Make the proposal scientifically rigorous and grounded in the provided literature. Flag literature gaps explicitly."""

        response = self.client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=4000,
            system=system_prompt,
            messages=[{"role": "user", "content": user_prompt}]
        )

        raw = response.content[0].text

        # Parse JSON - strip any markdown fencing
        clean = re.sub(r'^```(?:json)?\s*', '', raw.strip(), flags=re.MULTILINE)
        clean = re.sub(r'\s*```$', '', clean.strip(), flags=re.MULTILINE)

        proposal = None
        try:
            import json
            proposal = json.loads(clean)
        except Exception as e:
            print(f"  ⚠️  JSON parse error: {e}")
            # Return a minimal error structure
            proposal = {
                "title": f"Research Proposal: {topic}",
                "background": raw,
                "central_hypothesis": "See full text above",
                "null_hypothesis": "",
                "rationale": "",
                "specific_aims": [],
                "experimental_approach": "",
                "controls": "",
                "timeline": [],
                "literature_gaps": [],
                "references": []
            }

        self.current_proposal = proposal
        return proposal

    def refine_proposal(
        self,
        refinement_request: str,
        conversation_history: List[Dict]
    ) -> Tuple[Dict, str]:
        """
        Refine the current proposal based on user feedback.
        Returns (updated_proposal, summary_of_changes).
        """
        if not self.current_proposal:
            return None, "No proposal exists yet to refine."

        import json

        system_prompt = """You are refining a Drosophila research proposal based on researcher feedback.
        
Return a JSON object with TWO fields:
1. "updated_proposal": the complete updated proposal JSON (same structure as the original)
2. "changes_summary": a brief 1-2 sentence plain English description of what changed

Keep all sections that weren't mentioned. Only modify what the researcher asks for.
Return ONLY valid JSON, no markdown fencing."""

        user_prompt = f"""Current proposal:
{json.dumps(self.current_proposal, indent=2)}

Researcher's refinement request: "{refinement_request}"

Update the proposal accordingly and return the updated version with a changes summary."""

        response = self.client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=4000,
            system=system_prompt,
            messages=[{"role": "user", "content": user_prompt}]
        )

        raw = response.content[0].text
        clean = re.sub(r'^```(?:json)?\s*', '', raw.strip(), flags=re.MULTILINE)
        clean = re.sub(r'\s*```$', '', clean.strip(), flags=re.MULTILINE)

        try:
            result = json.loads(clean)
            updated = result.get('updated_proposal', self.current_proposal)
            summary = result.get('changes_summary', 'Proposal updated.')
            self.current_proposal = updated
            return updated, summary
        except Exception as e:
            print(f"  ⚠️  Refinement parse error: {e}")
            return self.current_proposal, "Refinement applied (see updated proposal)."

    def format_proposal_for_chat(self, proposal: Dict) -> str:
        """Format proposal as readable markdown for the chat interface."""
        md = []

        md.append(f"## 📄 {proposal.get('title', 'Research Proposal')}\n")

        md.append("### Background & Significance")
        md.append(proposal.get('background', '') + "\n")

        md.append("### Central Hypothesis")
        md.append(f"**H₁:** {proposal.get('central_hypothesis', '')}")
        if proposal.get('null_hypothesis'):
            md.append(f"**H₀:** {proposal.get('null_hypothesis', '')}")
        if proposal.get('rationale'):
            md.append(f"\n*Rationale:* {proposal.get('rationale', '')}")
        md.append("")

        aims = proposal.get('specific_aims', [])
        if aims:
            md.append("### Specific Aims")
            for aim in aims:
                md.append(f"\n**Aim {aim.get('aim_number', '?')}: {aim.get('title', '')}**")
                md.append(f"*Objective:* {aim.get('objective', '')}")
                md.append(f"*Approach:* {aim.get('approach', '')}")
                md.append(f"*Expected Outcomes:* {aim.get('expected_outcomes', '')}")
                md.append(f"*Potential Pitfalls:* {aim.get('potential_pitfalls', '')} — *Mitigation:* {aim.get('mitigation', '')}")

        if proposal.get('experimental_approach'):
            md.append("\n### Experimental Approach")
            md.append(proposal['experimental_approach'])

        if proposal.get('controls'):
            md.append("\n### Controls")
            md.append(proposal['controls'])

        gaps = proposal.get('literature_gaps', [])
        if gaps:
            md.append("\n### ⚠️ Literature Gaps (Research Opportunities)")
            for gap in gaps:
                md.append(f"- {gap}")

        timeline = proposal.get('timeline', [])
        if timeline:
            md.append("\n### Timeline")
            for phase in timeline:
                md.append(f"**{phase.get('phase', '')}:** {phase.get('milestones', '')}")

        refs = proposal.get('references', [])
        if refs:
            md.append("\n### References")
            for ref in refs:
                url = ref.get('url', '')
                citation = ref.get('citation', '')
                num = ref.get('number', '?')
                if url:
                    md.append(f"[{num}] {citation} — [Link]({url})")
                else:
                    md.append(f"[{num}] {citation}")

        return "\n".join(md)

    def get_proposal_for_export(self) -> Optional[Dict]:
        """Return current proposal for docx export."""
        return self.current_proposal