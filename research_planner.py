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
        Ask thorough clarifying questions before generating the proposal.
        More questions = better proposal. Err on the side of asking more.
        """
        gene_str = ', '.join(genes) if genes else 'not yet specified'

        prompt = f"""A grad student or postdoc wants to build a Drosophila research proposal about: "{topic}"
Identified genes so far: {gene_str}

You need to ask clarifying questions to write a specific, non-generic proposal. Ask about ALL of the following that are not already clear from the topic:
1. The specific phenotype or biological readout (e.g. lifespan? sleep duration? locomotor decline? metabolic rate? memory?)
2. Which tissue or cell type is the focus (e.g. fat body, neurons, gut, muscle, whole animal?)
3. What genetic tools they have access to (specific GAL4 lines, UAS-RNAi stocks, CRISPR capability, specific mutant stocks?)
4. Whether they have any preliminary data or a key paper already in mind that should anchor the proposal
5. The scope of the project (lab rotation, thesis chapter, grant application?)

Format as a friendly conversational paragraph asking 3-4 of these questions. More information leads to a better proposal so do not be shy about asking. Be direct and specific. Keep it under 120 words."""

        response = self.client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=300,
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

        # Assess literature quality to calibrate specificity
        has_strong_literature = len(papers) >= 3
        literature_warning = ""
        if not has_strong_literature:
            literature_warning = """
IMPORTANT: The literature provided is thin (fewer than 3 papers). This means:
- Do NOT pad the background with generic textbook descriptions of pathways
- Be explicit about what is actually known vs. what is assumed
- Flag major gaps honestly rather than writing around them
- The background should be shorter and more honest rather than longer and generic
- Ask the reader to consult additional literature in the gaps section
"""

        system_prompt = f"""You are an expert Drosophila researcher and scientific writing specialist helping grad students and postdocs build rigorous, specific research proposals.

Generate a complete research proposal in JSON format. The proposal must be specific and grounded — not generic.

CRITICAL RULES:
1. Cite the provided papers using [1], [2], etc. throughout. If fewer than 3 papers are provided, acknowledge this explicitly.
2. NEVER pad sections with generic pathway descriptions (e.g. "The insulin signaling pathway is conserved across species..."). Every sentence must be specific to the research question.
3. If literature is thin, flag it honestly: "⚠️ LITERATURE GAP: [specific gap]" and keep that section shorter rather than filling it with textbook knowledge.
4. Experimental approaches must name SPECIFIC tools: exact GAL4 drivers, specific UAS lines, specific assays with parameters (e.g. "CAFE assay measuring food intake at days 5, 15, 30").
5. The hypothesis must be falsifiable and specific — not "X affects aging" but "Loss of X in adult fat body neurons will extend median lifespan by >15% through upregulation of autophagy markers Atg1 and Atg8."
6. Pitfalls must be specific to THIS experiment, not generic lab advice.
7. Timeline must be realistic for the stated scope.
{literature_warning}

Return ONLY valid JSON with this exact structure:
{{
  "title": "Full descriptive proposal title",
  "background": "2-3 paragraphs covering what is specifically known about this gene/phenotype combination. Cite papers. Flag gaps honestly. Do NOT use generic pathway descriptions.",
  "central_hypothesis": "One specific, testable, falsifiable hypothesis with predicted magnitude of effect",
  "null_hypothesis": "The corresponding null hypothesis",
  "rationale": "Why THIS specific question is worth asking now — what gap does it fill? (2-3 sentences)",
  "specific_aims": [
    {{
      "aim_number": 1,
      "title": "Short aim title",
      "objective": "Precise objective for this aim",
      "approach": "Specific methods with named tools, strains, assay parameters, and sample sizes",
      "expected_outcomes": "Specific predicted results with measurable endpoints",
      "potential_pitfalls": "Pitfalls specific to this aim and approach",
      "mitigation": "Specific mitigation strategies"
    }}
  ],
  "experimental_approach": "Overview naming specific Drosophila strains, GAL4 drivers, UAS lines, and primary readouts",
  "controls": "Specific positive and negative controls for each major experiment",
  "timeline": [
    {{"phase": "Phase 1 (Months 1-4)", "milestones": "Specific milestones"}},
    {{"phase": "Phase 2 (Months 5-10)", "milestones": "Specific milestones"}},
    {{"phase": "Phase 3 (Months 11-16)", "milestones": "Specific milestones"}}
  ],
  "literature_gaps": ["Specific gap 1 with why it matters", "Specific gap 2"],
  "references": [
    {{"number": 1, "citation": "Author et al. (Year). Title. Journal. PMID: XXXXX", "url": "https://..."}}
  ]
}}"""

        context_quality = "STRONG" if has_strong_literature else "LIMITED"
        user_prompt = f"""Create a research proposal for the following:

Topic: {topic}
Key genes/proteins: {gene_str}{clarification_str}
Literature context quality: {context_quality} ({len(papers)} papers retrieved)
{paper_context}
{chat_context}

REQUIREMENTS:
- Generate 2-3 specific aims with concrete experimental approaches
- If literature context is LIMITED, be honest and specific about what is unknown rather than padding with generic pathway descriptions
- Every experimental approach must name specific strains, drivers, or assay parameters
- The hypothesis must include a predicted measurable outcome
- Flag real literature gaps — do not write around them with textbook knowledge"""


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