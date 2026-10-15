import anthropic
import json
import re
from typing import Any, Optional, Dict, List, Set, Tuple
from evidence_verifier import format_evidence_review, verify_claims

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


EXPERIMENT_DESIGN_PATTERNS = [
    r'design\s+(the\s+)?experiment\s+(for\s+)?aim\s*(\d+)',
    r'(do|run|build|create|show|generate)\s+(the\s+)?experiment\s+(for\s+)?aim\s*(\d+)',
    r'aim\s*(\d+)\s+experiment',
    r'(protocol|methods?|how\s+to\s+(run|do|conduct))\s+(for\s+)?aim\s*(\d+)',
    r'design\s+aim\s*(\d+)',
    r'(experiment|protocol)\s+design\s+(for\s+)?aim\s*(\d+)',
]

def detect_experiment_design_intent(message: str) -> Tuple[bool, int]:
    """
    Detect if user wants to design an experiment for a specific aim.
    Returns (is_experiment_design, aim_number).
    """
    import re
    msg_lower = message.lower()
    for pattern in EXPERIMENT_DESIGN_PATTERNS:
        m = re.search(pattern, msg_lower)
        if m:
            # Extract aim number from any group
            for g in m.groups():
                if g and g.isdigit():
                    return True, int(g)
            return True, 1  # default to aim 1 if no number found
    return False, 0


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


def _canonical_references(
    papers: List[Dict],
    citation_numbers: Optional[Set[int]] = None
) -> List[Dict]:
    """Build a numbered bibliography from retrieved records, not model output."""
    references = []
    for number, paper in enumerate(papers, 1):
        if citation_numbers is not None and number not in citation_numbers:
            continue
        pmid = str(paper.get("pmid", "")).strip()
        title = paper.get("title") or "No title"
        authors = paper.get("authors") or "Unknown"
        year = paper.get("year") or "N/A"
        citation = f"{authors} ({year}). {title}."
        if paper.get("journal"):
            citation += f" {paper['journal']}."
        if pmid:
            citation += f" PMID: {pmid}."
        if paper.get("doi"):
            citation += f" DOI: {paper['doi']}."

        url = paper.get("url") or (
            f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid.isdigit() else ""
        )
        reference = {
            "number": number,
            "citation": citation,
            "url": url,
            "pmid": pmid or None
        }
        for key in ("doi", "source", "_relevance_score"):
            if paper.get(key) is not None:
                reference["relevance_score" if key == "_relevance_score" else key] = paper[key]
        references.append(reference)
    return references


def paper_citation_token(paper: Dict, index: int) -> str:
    pmid = str(paper.get("pmid", "")).strip()
    return f"PMID:{pmid}" if pmid else f"SOURCE:{index}"


def validate_answer_citations(answer: str, papers: List[Dict]) -> str:
    """Replace model-authored reference lists with citations from retrieved records."""
    answer = re.sub(
        r"(?ims)^\s*(?:#{1,3}\s*)?(?:\*\*)?references(?:\*\*)?\s*:?\s*$[\s\S]*$",
        "",
        answer,
    ).rstrip()
    if not papers:
        answer = (
            "**Evidence limitation:** No relevant papers were retrieved. "
            "Any general background below is unverified and is not a "
            "paper-specific finding.\n\n"
            + answer
        )
    token_numbers = {
        paper_citation_token(paper, index): index
        for index, paper in enumerate(papers, 1)
    }
    cited_numbers = set()
    invalid_markers = []

    def replace_token(match):
        if match.group(3):
            invalid_markers.append(match.group(3))
            return ""
        token = f"{match.group(1)}:{match.group(2)}"
        number = token_numbers.get(token)
        if number is None:
            invalid_markers.append(token)
            return ""
        cited_numbers.add(number)
        return f"[{number}]"

    answer = re.sub(
        r"\[(?:(PMID|SOURCE):([A-Za-z0-9._-]+)|(\d+))\]",
        replace_token,
        answer,
    )

    if papers and not cited_numbers:
        answer += (
            "\n\n**Evidence limitation:** Retrieved papers were available, but the "
            "answer did not cite any of them. Treat its factual claims as unverified."
        )
    if invalid_markers:
        answer += (
            "\n\n**Evidence limitation:** Invalid citation marker(s) were removed "
            "because they did not match retrieved papers."
        )

    references = _canonical_references(papers, cited_numbers)
    if references:
        answer += "\n\n### References\n"
        for reference in references:
            answer += (
                f"\n[{reference['number']}] {reference['citation']} "
                f"[Source]({reference['url']})"
            )
    return answer


def _proposal_citation_numbers(proposal: Dict) -> Set[int]:
    numbers = set()
    for reference in proposal.get("references", []) or []:
        try:
            numbers.add(int(reference["number"]))
        except (KeyError, TypeError, ValueError):
            continue

    for text in re.findall(r'"(?:[^"\\]|\\.)*"', json.dumps(proposal)):
        numbers.update(int(number) for number in re.findall(r"\[(\d+)\]", text))
    return numbers


def _remove_unverified_citation_markers(value: Any, valid_numbers: Set[int]) -> Any:
    if isinstance(value, str):
        return re.sub(
            r"\[(\d+)\]",
            lambda match: match.group(0) if int(match.group(1)) in valid_numbers else "",
            value
        )
    if isinstance(value, list):
        return [_remove_unverified_citation_markers(item, valid_numbers) for item in value]
    if isinstance(value, dict):
        return {
            key: _remove_unverified_citation_markers(item, valid_numbers)
            for key, item in value.items()
        }
    return value


class ResearchPlanner:
    def __init__(self, anthropic_client: anthropic.Anthropic):
        self.client = anthropic_client
        self.current_proposal: Optional[Dict] = None
        self.proposal_context: Dict = {}  # topic, genes, papers from prior chat

    def set_context_from_chat(self, topic: str, genes: List[str], papers: List[Dict], clarifications: str = ""):
        """Called by the assistant to seed the planner with chat context."""
        self.proposal_context = {
            'topic': topic,
            'genes': genes,
            'papers': papers,
            'clarifications': clarifications
        }

    # Scope/admin words that should never be treated as biological topics
    SCOPE_WORDS = [
        'rotation', 'rotational', 'rotate', 'rotating',
        'thesis', 'dissertation', 'chapter',
        'grant', 'application', 'proposal',
        'project', 'study', 'research',
        'rotation project', 'thesis chapter', 'grant application',
        'short project', 'lab rotation',
        'all capabilities', 'all tools', 'all techniques',
        'access to all', 'full access', 'full capabilities',
    ]

    def strip_scope_words(self, text: str) -> str:
        """Remove scope/admin words that describe project type, not biology."""
        clean = text.lower()
        for word in self.SCOPE_WORDS:
            clean = re.sub(r'\b' + re.escape(word) + r'\b', '', clean)
        return clean.strip()

    def extract_topic_from_message(self, message: str) -> str:
        """Pull the core research topic out of the user's message."""
        # Strip planning keywords to get the topic
        clean = message.lower()
        for phrase in PLANNING_INTENT_STRONG + PLANNING_INTENT_CONTEXTUAL:
            clean = clean.replace(phrase, '')

        # Remove scope/admin words first
        clean = self.strip_scope_words(clean)

        # Remove filler words
        fillers = ['for', 'about', 'on', 'regarding', 'around', 'this', 'the', 'a', 'an', 'me', 'please', 'can you', 'could you', 'i want to', 'i would like to']
        for filler in fillers:
            clean = re.sub(r'\b' + filler + r'\b', '', clean)

        topic = clean.strip().strip(',').strip()
        return topic if topic else message.strip()

    def generate_clarifying_questions(self, topic: str, genes: List[str],
                                       conversation_history: List[Dict] = None) -> Optional[str]:
        """
        Ask clarifying questions before generating the proposal.
        Returns None if the existing context is already sufficient.
        """
        gene_str = ', '.join(genes) if genes else 'not yet specified'

        context_summary = ""
        if conversation_history:
            recent = conversation_history[-6:]
            lines = []
            for msg in recent:
                role = "Researcher" if msg['role'] == 'user' else "Assistant"
                content = msg['content']
                if '======' in content:
                    content = content.split('\n\n')[0]
                lines.append(f"{role}: {content[:200]}")
            context_summary = "\n".join(lines)

        prompt = f"""A grad student wants to build a Drosophila research proposal about: "{topic}"
Identified genes: {gene_str}

Recent conversation (use this to avoid asking things already answered):
{context_summary if context_summary else "No prior conversation."}

Ask ONLY about things not already clear from the topic and conversation above:
1. Specific phenotype/readout (lifespan, locomotor decline, memory, sleep, etc.)
2. Tissue or cell type focus (fat body, neurons, gut, muscle, whole animal)
3. Genetic tools available (specific GAL4 lines, UAS-RNAi, CRISPR)
4. Any preliminary data or anchor paper

If the topic already specifies 3 or more of these, respond with exactly: SUFFICIENT
Otherwise ask 2-3 focused questions in a friendly paragraph under 100 words."""

        response = self.client.messages.create(
            model="claude-sonnet-4-6",
            max_tokens=300,
            messages=[{"role": "user", "content": prompt}]
        )
        result = response.content[0].text.strip()
        return None if result.startswith("SUFFICIENT") else result

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
IMPORTANT: The retrieved literature is limited (fewer than 3 papers). This means:
- Do not invent additional sources or rely on unprovided sources from memory
- Keep literature claims narrowly tied to the supplied abstracts and cite each factual claim
- Explicitly label any inference or proposed experimental assumption as unverified
- State when the available evidence cannot support a conclusion and identify it as a gap
"""

        system_prompt = f"""You are an expert Drosophila researcher and scientific writing specialist helping grad students and postdocs build rigorous, specific research proposals.

Generate a complete research proposal in JSON format. The proposal must be specific and grounded — not generic.

CRITICAL RULES:
1. Cite papers from the PROVIDED LITERATURE LIST using numeric citations [1], [2], etc. The number must match the paper's position in that list. Do NOT invent citations or add references that are not in the provided list.
2. Cite each factual statement about prior research immediately; a real citation is not proof that the cited abstract supports a claim.
3. Do NOT include unsourced factual statements as established findings. If the provided literature does not support a claim, mark the gap or state that evidence is insufficient.
4. Clearly distinguish sourced findings, biological inferences, hypotheses, and proposed experimental choices. Label inferences and proposed choices explicitly.
5. Never claim that a strain, stock number, interaction, result, or novelty claim has been verified unless the provided literature supports that exact claim.
6. Experimental approaches may propose specific tools and assays, but label them as proposals requiring independent validation, not established facts.
7. If literature is thin, flag it honestly and keep literature claims narrow rather than filling gaps from memory.
8. Pitfalls must be specific to THIS experiment, not generic lab advice.
9. Timeline must be realistic for the stated scope.
10. Include every material factual claim in "evidence_claims" as an exact quote from the proposal text, with an evidence_type of finding, inference, general_background, proposal, or prediction, and source_ids containing the exact PMID/SOURCE identifiers for cited papers. Include unsupported factual claims too; do not omit claims because evidence is weak.
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
  "evidence_claims": [
    {{"claim": "Exact claim text from another proposal field", "evidence_type": "finding|inference|general_background|proposal|prediction", "source_ids": ["PMID or SOURCE token"]}}
  ],
  "references": [
    {{"number": 1, "citation": "Author et al. (Year). Title. PMID: XXXXX", "url": "https://..."}}
  ]
}}"""

        context_quality = "STRONG" if has_strong_literature else "LIMITED"
        # Build a strict topic anchor to prevent drift
        topic_anchor = f"""STRICT TOPIC ANCHOR — DO NOT DEVIATE:
The researcher explicitly asked about: "{topic}"
Researcher clarifications: "{user_clarifications if user_clarifications else 'none provided'}"
Genes mentioned by researcher: {gene_str}

CRITICAL: The proposal MUST be about the topic and genes stated above.
Do NOT invent new genes, pathways, or research directions not mentioned by the researcher or present in the provided literature.
If the literature pulls in tangential topics, ignore them — stay focused on what the researcher asked for.
If you do not have enough specific literature on the exact topic, say so explicitly in the background rather than pivoting to a related topic you know more about."""

        user_prompt = f"""Create a research proposal for the following:

{topic_anchor}

Literature context quality: {context_quality} ({len(papers)} papers retrieved)
{paper_context}
{chat_context}

REQUIREMENTS:
- Generate 2-3 specific aims directly addressing the stated topic
- Stay anchored to the researcher's stated topic — do not drift to related but different questions
- If literature context is LIMITED for the exact topic, be honest rather than padding or pivoting
- Every experimental approach must name specific strains, drivers, or assay parameters
- The hypothesis must include a predicted measurable outcome
- Flag real literature gaps — do not write around them with invented research directions"""


        response = self.client.messages.create(
            model="claude-sonnet-4-6",
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

        citation_numbers = _proposal_citation_numbers(proposal)
        available_numbers = set(range(1, min(len(papers), 10) + 1))
        valid_numbers = citation_numbers & available_numbers
        proposal = _remove_unverified_citation_markers(proposal, valid_numbers)
        proposal["references"] = _canonical_references(papers[:10], valid_numbers)
        proposal["evidence_status"] = (
            "retrieved_literature" if papers else "no_retrieved_literature"
        )
        proposal["evidence_review"] = verify_claims(
            json.dumps({
                key: value
                for key, value in proposal.items()
                if key not in {"references", "evidence_review", "evidence_claims"}
            }, ensure_ascii=False),
            papers[:10],
            self.client,
            structured_claims=proposal.get("evidence_claims"),
        )
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
            model="claude-sonnet-4-6",
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
            updated["references"] = self.current_proposal.get("references", [])
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

        if proposal.get('evidence_status') == 'no_retrieved_literature':
            md.append(
                "\n> **Evidence limitation:** No papers were retrieved. "
                "Treat literature-based claims in this proposal as unverified."
            )
        elif proposal.get('evidence_status') == 'retrieved_literature':
            md.append(
                "\n> Literature grounding is limited to the retrieved records "
                "and abstracts. Citation identity is checked, but each source "
                "should be consulted to confirm claim-level support."
            )

        if proposal.get('evidence_review'):
            md.append("\n" + format_evidence_review(proposal["evidence_review"]))

        findings = proposal.get('literature_findings') or {}
        for label, key in (
            ("Genes", "genes"),
            ("Pathways", "pathways"),
            ("Human Orthologs", "human_orthologs"),
            ("Disease Context", "disease_context"),
        ):
            value = findings.get(key)
            if value:
                if isinstance(value, dict):
                    value = "; ".join(f"{name}: {detail}" for name, detail in value.items())
                elif isinstance(value, list):
                    value = "; ".join(
                        f"{item.get('name', item.get('gene', ''))}: "
                        f"{item.get('description', item.get('function', ''))}".strip(": ")
                        if isinstance(item, dict) else str(item)
                        for item in value
                    )
                md.append(f"\n### {label}")
                md.append(str(value))

        fly_lines = proposal.get('fly_lines_needed') or []
        if fly_lines:
            md.append("\n### Fly Lines Needed")
            for line in fly_lines:
                if isinstance(line, dict):
                    detail = f"`{line.get('line', 'Unspecified')}`"
                    if line.get('purpose'):
                        detail += f" — {line['purpose']}"
                    if line.get('source'):
                        detail += f" (Source: {line['source']})"
                    md.append(f"- {detail}")
                else:
                    md.append(f"- {line}")

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

        review = proposal.get('critic_review') or {}
        if review:
            scores = review.get('scores') or {}
            if scores:
                md.append("\n### Independent Review")
                md.append("**Critic scores:** " + " · ".join(
                    f"{name}: {score}/5" for name, score in scores.items()
                ))
            for label, key in (("Issues", "issues"), ("Suggestions", "suggestions")):
                items = review.get(key) or []
                if items:
                    md.append(f"\n**{label}:**")
                    md.extend(f"- {item}" for item in items)
            if review.get('reasoning'):
                md.append(f"\n*{review['reasoning']}*")

        refs = proposal.get('references', [])
        if refs:
            md.append("\n### References")
            for ref in refs:
                url = ref.get('url', '')
                citation = ref.get('citation', '')
                num = ref.get('number', '?')
                if url:
                    details = []
                    if ref.get('source'):
                        details.append(f"source: {ref['source']}")
                    if isinstance(ref.get('relevance_score'), (int, float)):
                        details.append(
                            f"TF-IDF text-match heuristic: {ref['relevance_score']:.2f}"
                        )
                    note = f" ({'; '.join(details)})" if details else ""
                    doi_link = (
                        f" — [DOI](https://doi.org/{ref['doi']})"
                        if ref.get('doi') else ""
                    )
                    md.append(f"[{num}] {citation} — [Link]({url}){doi_link}{note}")
                else:
                    md.append(f"[{num}] {citation}")

        return "\n".join(md)


    def generate_experiment_design(self, aim: Dict, proposal: Dict) -> Dict:
        """
        Generate a detailed experiment design for a specific aim.
        Returns a structured dict with all protocol sections.
        """
        import json

        aim_num = aim.get('aim_number', 1)
        aim_title = aim.get('title', '')

        # Pull key context from the proposal
        genes = []
        for ref in proposal.get('references', []):
            pass  # not needed here

        system_prompt = """You are an expert Drosophila researcher generating a detailed, bench-ready experimental protocol for a specific aim from a research proposal.

Generate a complete experiment design in JSON format. This should be specific enough that a graduate student can execute it directly.

CRITICAL RULES:
1. Cross schemes must specify exact genotypes, which sex carries which chromosome, and how many generations
2. All timepoints must be specific (exact days, not "early" or "late")
3. Data collection sheet must be a structured table format (use pipe | notation for columns)
4. Go/no-go criteria must have specific measurable thresholds
5. Statistical analysis must name the exact test and software
6. Be specific about Drosophila husbandry details (vial size, fly density, food type, temperature)

Return ONLY valid JSON with this exact structure:
{
  "aim_number": 1,
  "aim_title": "Short title",
  "cross_scheme": {
    "overview": "Brief description of crossing strategy",
    "parental_genotypes": [
      {"line": "Line name", "genotype": "Full genotype notation", "source": "Stock center or lab stock", "notes": "Any special handling"}
    ],
    "generations": [
      {"generation": "P (Parental)", "cross": "♀ Genotype A × ♂ Genotype B", "instructions": "Step by step what to do", "timing": "How long, when to collect"}
    ],
    "expected_progeny": "Description of progeny classes and how to identify experimental vs control flies",
    "balancer_notes": "Any notes about balancer chromosomes or selection"
  },
  "vial_setup": {
    "vial_type": "Standard 25mm vial or 50ml bottle",
    "fly_density": "X flies per vial",
    "food": "Standard cornmeal-yeast-agar or specific recipe",
    "temperature": "25°C",
    "light_cycle": "12:12 LD",
    "humidity": "~60%",
    "flipping_schedule": "Every X days",
    "special_conditions": "Any drug treatments, dietary modifications etc"
  },
  "timepoints": [
    {"day": 0, "action": "What to do", "what_to_collect": "Samples/measurements", "notes": "Any special considerations"}
  ],
  "data_collection_sheet": {
    "description": "What this sheet tracks",
    "columns": ["Col1", "Col2", "Col3"],
    "example_row": ["Example value 1", "Example value 2", "Example value 3"],
    "scoring_criteria": "How to score ambiguous cases"
  },
  "statistical_analysis": {
    "primary_test": "Test name and why",
    "software": "R / Prism / Python",
    "sample_size_justification": "Why n=X is sufficient",
    "censoring_criteria": "When to censor a fly (escaped, contamination etc)",
    "multiple_comparisons": "How to handle if applicable"
  },
  "go_no_go_criteria": [
    {"checkpoint": "What to check", "timing": "When (e.g. Day 5)", "threshold": "Specific measurable threshold", "if_fail": "What to do if it fails"}
  ],
  "troubleshooting": [
    {"problem": "Common problem", "likely_cause": "Why it happens", "solution": "What to do"}
  ],
  "estimated_duration": "Total time from cross setup to data collection complete",
  "materials_needed": ["Item 1", "Item 2"]
}"""

        user_prompt = f"""Generate a detailed experiment design for this specific aim:

AIM {aim_num}: {aim_title}
Objective: {aim.get('objective', '')}
Approach (from proposal): {aim.get('approach', '')}
Expected Outcomes: {aim.get('expected_outcomes', '')}
Potential Pitfalls: {aim.get('potential_pitfalls', '')}
Mitigation: {aim.get('mitigation', '')}

PROPOSAL CONTEXT:
Title: {proposal.get('title', '')}
Experimental Approach: {proposal.get('experimental_approach', '')}
Controls: {proposal.get('controls', '')}

Generate a complete bench-ready protocol. Be specific about genotypes, timing, and thresholds."""

        response = self.client.messages.create(
            model="claude-sonnet-4-6",
            max_tokens=4000,
            system=system_prompt,
            messages=[{"role": "user", "content": user_prompt}]
        )

        raw = response.content[0].text
        clean = re.sub(r'^```(?:json)?\s*', '', raw.strip(), flags=re.MULTILINE)
        clean = re.sub(r'\s*```$', '', clean.strip(), flags=re.MULTILINE)

        try:
            design = json.loads(clean)
        except Exception as e:
            print(f"  ⚠️  Experiment design JSON parse error: {e}")
            design = {
                "aim_number": aim_num,
                "aim_title": aim_title,
                "error": "Failed to parse structured design",
                "raw": raw
            }

        return design

    def format_experiment_design(self, design: Dict) -> str:
        """Format experiment design as markdown for display."""
        md = []

        md.append(f"## 🧪 Experiment Design — Aim {design.get('aim_number', '?')}: {design.get('aim_title', '')}")
        md.append("")

        if design.get('error'):
            md.append(design.get('raw', 'Error generating design.'))
            return "\n".join(md)

        # Estimated duration
        if design.get('estimated_duration'):
            md.append(f"⏱️ **Estimated Duration:** {design['estimated_duration']}")
            md.append("")

        # Cross scheme
        cs = design.get('cross_scheme', {})
        if cs:
            md.append("### 🧬 Cross Scheme")
            md.append(cs.get('overview', ''))
            md.append("")

            parents = cs.get('parental_genotypes', [])
            if parents:
                md.append("**Parental Lines:**")
                for p in parents:
                    md.append(f"- **{p.get('line', '')}:** `{p.get('genotype', '')}` — {p.get('source', '')}  ")
                    if p.get('notes'):
                        md.append(f"  *{p['notes']}*")
                md.append("")

            gens = cs.get('generations', [])
            if gens:
                md.append("**Crossing Generations:**")
                for g in gens:
                    md.append(f"**{g.get('generation', '')}:** {g.get('cross', '')}")
                    md.append(f"- {g.get('instructions', '')}")
                    md.append(f"- *Timing: {g.get('timing', '')}*")
                md.append("")

            if cs.get('expected_progeny'):
                md.append(f"**Expected Progeny:** {cs['expected_progeny']}")
            if cs.get('balancer_notes'):
                md.append(f"**Balancer Notes:** {cs['balancer_notes']}")
            md.append("")

        # Vial setup
        vs = design.get('vial_setup', {})
        if vs:
            md.append("### 🫙 Vial Setup & Husbandry")
            for key, label in [
                ('vial_type', 'Vial Type'), ('fly_density', 'Fly Density'),
                ('food', 'Food'), ('temperature', 'Temperature'),
                ('light_cycle', 'Light Cycle'), ('humidity', 'Humidity'),
                ('flipping_schedule', 'Flipping Schedule'), ('special_conditions', 'Special Conditions')
            ]:
                if vs.get(key):
                    md.append(f"- **{label}:** {vs[key]}")
            md.append("")

        # Timepoints
        timepoints = design.get('timepoints', [])
        if timepoints:
            md.append("### 📅 Timepoints & Sampling Schedule")
            for tp in timepoints:
                md.append(f"**Day {tp.get('day', '?')}:** {tp.get('action', '')}")
                if tp.get('what_to_collect'):
                    md.append(f"- *Collect/Measure:* {tp['what_to_collect']}")
                if tp.get('notes'):
                    md.append(f"- *Note:* {tp['notes']}")
            md.append("")

        # Data collection sheet
        dcs = design.get('data_collection_sheet', {})
        if dcs:
            md.append("### 📊 Data Collection Sheet")
            md.append(dcs.get('description', ''))
            cols = dcs.get('columns', [])
            example = dcs.get('example_row', [])
            if cols:
                md.append("| " + " | ".join(cols) + " |")
                md.append("| " + " | ".join(["---"] * len(cols)) + " |")
                if example:
                    md.append("| " + " | ".join(str(e) for e in example) + " |")
            if dcs.get('scoring_criteria'):
                md.append(f"\n*Scoring criteria:* {dcs['scoring_criteria']}")
            md.append("")

        # Statistical analysis
        sa = design.get('statistical_analysis', {})
        if sa:
            md.append("### 📈 Statistical Analysis")
            for key, label in [
                ('primary_test', 'Primary Test'), ('software', 'Software'),
                ('sample_size_justification', 'Sample Size'), ('censoring_criteria', 'Censoring Criteria'),
                ('multiple_comparisons', 'Multiple Comparisons')
            ]:
                if sa.get(key):
                    md.append(f"- **{label}:** {sa[key]}")
            md.append("")

        # Go/no-go criteria
        gng = design.get('go_no_go_criteria', [])
        if gng:
            md.append("### 🚦 Go/No-Go Criteria")
            for g in gng:
                md.append(f"**{g.get('checkpoint', '')}** *(Check: {g.get('timing', '')})*")
                md.append(f"- Threshold: {g.get('threshold', '')}")
                md.append(f"- If fail: {g.get('if_fail', '')}")
            md.append("")

        # Troubleshooting
        ts = design.get('troubleshooting', [])
        if ts:
            md.append("### 🔧 Troubleshooting")
            for t in ts:
                md.append(f"**Problem:** {t.get('problem', '')}")
                md.append(f"- *Cause:* {t.get('likely_cause', '')}")
                md.append(f"- *Solution:* {t.get('solution', '')}")
            md.append("")

        # Materials
        materials = design.get('materials_needed', [])
        if materials:
            md.append("### 🧴 Materials Needed")
            for m in materials:
                md.append(f"- {m}")

        return "\n".join(md)

    def get_proposal_for_export(self) -> Optional[Dict]:
        """Return current proposal for docx export."""
        return self.current_proposal