import anthropic
import functools
import html
import os
import urllib.request
from Bio import Entrez
import json
import requests
from pathlib import Path
from typing import Optional, Dict, List
import re
from evidence_verifier import (
    format_evidence_review,
    parse_structured_chat_response,
    verify_and_revise_claims,
)
from research_planner import (
    ResearchPlanner,
    detect_experiment_design_intent,
    detect_planning_intent,
    paper_citation_token,
    validate_answer_citations,
)

# Configure Entrez
Entrez.email = "jordanszt@icloud.com"
Entrez.api_key = os.environ.get("NCBI_API_KEY", "")

# Bio.Entrez calls urlopen() with no timeout, so one stalled NCBI request hangs the
# whole request (or calibration run) forever. Give it a deadline instead.
NCBI_TIMEOUT_SECONDS = 30
Entrez.urlopen = functools.partial(urllib.request.urlopen, timeout=NCBI_TIMEOUT_SECONDS)

# Abstracts were cut at 500 characters, which hid later sentences from the answer model and
# the evidence check (supported claims came back "uncertain"). Keep them whole, with a cap
# well above typical abstract length to bound prompt size.
MAX_ABSTRACT_CHARS = 4000


def _clean_abstract(abstract: str) -> str:
    # PubMed abstracts carry markup ("The <i>Drosophila</i> leg"); quotes copied from them come
    # back as plain text and failed the evidence check's verbatim match.
    abstract = html.unescape(re.sub(r"<[^>]+>", "", abstract))
    return abstract[:MAX_ABSTRACT_CHARS] + '...' if len(abstract) > MAX_ABSTRACT_CHARS else abstract


# Anthropic's default is 10 minutes per call; fail sooner so the UI can report it.
ANTHROPIC_TIMEOUT_SECONDS = 300.0

# Load full FlyBase gene database at module startup
_GENE_DB = {}
_gene_db_path = Path(__file__).parent / 'gene_db.json'

if _gene_db_path.exists():
    with open(_gene_db_path) as f:
        _GENE_DB = json.load(f)
    print(f"✅ Loaded gene database: {len(_GENE_DB)} entries")
else:
    print("⚠️  gene_db.json not found — gene lookup will be limited")

# Many FlyBase symbols are also everyday words ("in", "how", "for", "path").
# These are skipped unless the text explicitly marks them as a gene.
_COMMON_WORDS = {
    'a', 'about', 'after', 'all', 'also', 'an', 'and', 'any', 'are', 'as', 'at', 'be', 'been',
    'before', 'between', 'both', 'but', 'by', 'can', 'could', 'do', 'does', 'during', 'each',
    'for', 'from', 'had', 'has', 'have', 'he', 'her', 'him', 'his', 'how', 'if', 'in', 'into',
    'is', 'it', 'its', 'may', 'me', 'might', 'more', 'most', 'my', 'no', 'not', 'of', 'on',
    'or', 'other', 'our', 'out', 'over', 'so', 'some', 'such', 'than', 'that', 'the', 'their',
    'them', 'then', 'there', 'these', 'they', 'this', 'those', 'through', 'to', 'under', 'up',
    'us', 'was', 'we', 'were', 'what', 'when', 'where', 'which', 'while', 'who', 'why', 'will',
    'with', 'within', 'would', 'you', 'your', 'tell', 'n', 'e', 'g', 'eg', 'ie', 'al', 'et',
    'vs', 'via', 'dr', 'mr', 'st', 'fig',
    # research vocabulary
    'aging', 'gene', 'genes', 'protein', 'proteins', 'sleep', 'memory', 'behavior', 'cell',
    'cells', 'path', 'pathway', 'role', 'set', 'red', 'step', 'side', 'link', 'run', 'key',
    'model', 'test', 'term', 'list', 'lab', 'net', 'level', 'levels', 'cycle', 'clock',
    'broad', 'dorsal', 'cortex', 'expansion', 'expanded', 'drop', 'flower', 'diamond',
    'arrow', 'cut', 'loss', 'gain', 'fat', 'body', 'wing', 'eye', 'gut', 'brain', 'stem',
    'new', 'novel', 'using', 'effect', 'effects', 'fate', 'growth', 'stress', 'signal',
    'signaling', 'adult', 'larval', 'larva', 'fly', 'flies', 'study', 'mutant', 'mutants',
}

# Single-word FlyBase gene names ("Notch", "hippo", "yorkie") mapped to their database keys.
_GENE_NAME_INDEX = {}
for _key, _entry in _GENE_DB.items():
    _name = (_entry.get('name') or '').lower()
    if len(_name) >= 4 and ' ' not in _name and _name not in _COMMON_WORDS:
        _GENE_NAME_INDEX.setdefault(_name, _key)


def search_flybase_fixed(gene_name: str) -> Optional[Dict]:
    gene_lower = gene_name.strip().lower()

    entry = _GENE_DB.get(gene_lower)
    if entry:
        print(f"  ✅ Found: {entry['symbol']} ({entry['fbgn']}) [Local DB]")
        return {
            'symbol': entry['symbol'],
            'name': entry.get('name', entry['symbol']),
            'fbgn': entry['fbgn'],
            'summary': f"Gene symbol: {entry['symbol']}.",
            'synonyms': [],
            'url': entry['url']
        }

    print(f"  ℹ️  '{gene_name}' not found in gene database")
    return None


def _pubmed_doi(article: Dict) -> str:
    for identifier in article.get('PubmedData', {}).get('ArticleIdList', []):
        if getattr(identifier, 'attributes', {}).get('IdType') == 'doi':
            return str(identifier)
    return ""


def _publication_year(article_data: Dict) -> str:
    dates = article_data.get('ArticleDate') or []
    if dates and dates[0].get('Year'):
        return str(dates[0]['Year'])
    pub_date = article_data.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
    if pub_date.get('Year'):
        return str(pub_date['Year'])
    medline_date = str(pub_date.get('MedlineDate', ''))
    year_match = re.search(r'\b(19|20)\d{2}\b', medline_date)
    return year_match.group(0) if year_match else 'N/A'


class DrosophilaAssistant:
    def __init__(self, api_key):
        self.client = anthropic.Anthropic(api_key=api_key, timeout=ANTHROPIC_TIMEOUT_SECONDS)
        self.conversation_history = []
        self.mentioned_genes = set()
        self.last_topic = ""
        self.last_papers = []
        self.last_genes = []

        # Research planner — integrated but separate
        self.planner = ResearchPlanner(self.client)
        self.planning_mode = False          # Are we currently in a planning session?
        self.awaiting_clarification = False # Did we ask a clarifying question?

    def determine_paper_count(self, user_query: str) -> int:
        query_lower = user_query.lower()
        broad = ['top', 'list', 'overview', 'review', 'genes in', 'pathways', 'mechanisms', 'role of', 'involved in', 'best', 'major']
        specific = ['what is', 'tell me about', 'how does', 'function of']
        if any(ind in query_lower for ind in broad):
            return 8
        if any(ind in query_lower for ind in specific):
            return 5
        return 6

    def get_flybase_publications(self, fbgn: str, max_results: int = 10) -> List[Dict]:
        """Get publications for a gene via NCBI eUtils gene→pubmed link."""
        print(f"  📚 Fetching publications for {fbgn} via NCBI...")
        try:
            # Step 1: resolve FBgn → NCBI Gene ID
            handle = Entrez.esearch(db="gene", term=f"{fbgn}[Gene ID] AND 7227[taxid]", retmax=1)
            record = Entrez.read(handle)
            handle.close()
            gene_ids = record.get("IdList", [])

            # Fall back to symbol search if FBgn cross-ref not found
            if not gene_ids:
                symbol = next(
                    (v['symbol'] for k, v in _GENE_DB.items() if v.get('fbgn') == fbgn),
                    None
                )
                if symbol:
                    handle = Entrez.esearch(db="gene", term=f"{symbol}[gene] AND 7227[taxid]", retmax=1)
                    record = Entrez.read(handle)
                    handle.close()
                    gene_ids = record.get("IdList", [])

            if not gene_ids:
                print(f"    ⚠️  No NCBI Gene ID found for {fbgn}")
                return []

            ncbi_gene_id = gene_ids[0]

            # Step 2: gene → pubmed elink (curated NCBI gene-publication links)
            handle = Entrez.elink(dbfrom="gene", db="pubmed", id=ncbi_gene_id)
            link_records = Entrez.read(handle)
            handle.close()

            pmids = []
            for linkset in link_records:
                for link_db in linkset.get("LinkSetDb", []):
                    if link_db.get("LinkName") == "gene_pubmed":
                        pmids = [str(l["Id"]) for l in link_db.get("Link", [])]
                        break
                if pmids:
                    break

            if not pmids:
                print(f"    ℹ️  No linked publications for {fbgn}, falling back to PubMed search")
                return self._fallback_pubmed_search(fbgn, max_results)

            pmids = pmids[:max_results]

            # Step 3: fetch abstracts
            handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            publications = []
            for article in records.get('PubmedArticle', []):
                try:
                    medline = article['MedlineCitation']
                    article_data = medline['Article']
                    title = article_data.get('ArticleTitle', 'No title')
                    abstract = article_data.get('Abstract', {}).get('AbstractText', [''])
                    abstract = ' '.join(str(a) for a in abstract) if isinstance(abstract, list) else str(abstract)
                    authors = []
                    for author in article_data.get('AuthorList', [])[:3]:
                        if 'LastName' in author:
                            authors.append(f"{author['LastName']} {author.get('Initials', '')}")
                    pmid = str(medline['PMID'])
                    year = _publication_year(article_data)
                    publications.append({
                        'title': title,
                        'authors': ', '.join(authors) + ' et al.' if authors else 'Unknown',
                        'year': year,
                        'journal': article_data.get('Journal', {}).get('Title', ''),
                        'doi': _pubmed_doi(article),
                        'pmid': pmid,
                        'abstract': _clean_abstract(abstract),
                        'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        'source': 'FlyBase-linked PubMed',
                    })
                except Exception:
                    continue

            print(f"  ✅ Retrieved {len(publications)} publications for {fbgn}")
            return publications

        except Exception as e:
            print(f"  ⚠️  NCBI lookup failed for {fbgn}: {e}")
            return self._fallback_pubmed_search(fbgn, max_results)

    def _fallback_pubmed_search(self, fbgn: str, max_results: int) -> List[Dict]:
        try:
            entry = next((v for v in _GENE_DB.values() if v.get('fbgn') == fbgn), None)
            gene_name = entry['symbol'] if entry else fbgn
            handle = Entrez.esearch(db="pubmed", term=f"Drosophila AND {gene_name}", retmax=max_results + 2, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            pmids = record.get("IdList", [])[:max_results]
            if not pmids:
                return []
            handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            publications = []
            for article in records.get('PubmedArticle', []):
                try:
                    medline = article['MedlineCitation']
                    article_data = medline['Article']
                    title = article_data.get('ArticleTitle', 'No title')
                    abstract = article_data.get('Abstract', {}).get('AbstractText', [''])
                    abstract = ' '.join(str(a) for a in abstract) if isinstance(abstract, list) else str(abstract)
                    authors = []
                    for author in article_data.get('AuthorList', [])[:3]:
                        if 'LastName' in author:
                            authors.append(f"{author['LastName']} {author.get('Initials', '')}")
                    pmid = str(medline['PMID'])
                    year = _publication_year(article_data)
                    publications.append({
                        'title': title,
                        'authors': ', '.join(authors) + ' et al.' if authors else 'Unknown',
                        'year': year,
                        'journal': article_data.get('Journal', {}).get('Title', ''),
                        'doi': _pubmed_doi(article),
                        'pmid': pmid,
                        'abstract': _clean_abstract(abstract),
                        'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        'source': 'PubMed',
                    })
                except Exception:
                    continue
            return publications
        except Exception:
            return []

    def filter_relevant_papers(self, papers: List[Dict], topic: str, genes: List[str], user_clarifications: str = "") -> List[Dict]:
        """
        Filter and rank papers for relevance using:
        1. Hard filter: must mention Drosophila (drops off-organism papers)
        2. TF-IDF cosine similarity: ranks remaining papers by lexical
           similarity to the researcher's query — discriminates by rare
           meaningful terms rather than raw keyword counts.
        """
        if not papers:
            return []

        from sklearn.feature_extraction.text import TfidfVectorizer
        from sklearn.metrics.pairwise import cosine_similarity
        import numpy as np

        drosophila_terms = {'drosophila', 'melanogaster', 'fruit fly', 'flies', 'drosophilid'}

        # Hard filter: must mention Drosophila OR come from FlyBase
        # FlyBase papers are already organism-specific by definition
        drosophila_papers = []
        for paper in papers:
            title = (paper.get('title', '') or '').lower()
            abstract = (paper.get('abstract', '') or '').lower()
            combined = title + ' ' + abstract
            from_flybase = paper.get('source', '').startswith('FlyBase')
            if from_flybase or any(term in combined for term in drosophila_terms):
                drosophila_papers.append(paper)
            else:
                print(f"    ❌ Dropped (no Drosophila): {paper.get('title', '')[:60]}")

        if not drosophila_papers:
            return []

        # Build query from topic + clarifications + gene names
        # Weight genes more heavily by repeating them
        gene_str = ' '.join(genes * 3)  # repeat separated gene terms for emphasis
        query = ' '.join(part for part in (topic, user_clarifications, gene_str) if part).strip()

        # Build corpus: query + all paper texts
        paper_texts = []
        for p in drosophila_papers:
            title = p.get('title', '') or ''
            abstract = p.get('abstract', '') or ''
            paper_texts.append(f"{title} {abstract}")

        corpus = [query] + paper_texts

        # TF-IDF vectorize — use sublinear_tf to reduce dominance of
        # very frequent terms, min_df=1 since corpus is small
        try:
            vectorizer = TfidfVectorizer(
                sublinear_tf=True,
                min_df=1,
                stop_words='english',
                ngram_range=(1, 2)  # unigrams + bigrams catch "stem cell", "fat body" etc
            )
            tfidf_matrix = vectorizer.fit_transform(corpus)

            # Query is index 0, papers are 1..n
            query_vec = tfidf_matrix[0]
            paper_vecs = tfidf_matrix[1:]

            similarities = cosine_similarity(query_vec, paper_vecs)[0]

            # Attach scores and sort
            for i, paper in enumerate(drosophila_papers):
                paper['_relevance_score'] = float(similarities[i])

            drosophila_papers.sort(key=lambda p: p.get('_relevance_score', 0), reverse=True)

            # No threshold — keep all Drosophila papers, just reorder by TF-IDF score.
            # The hard Drosophila filter above already ensures organism relevance.
            # TF-IDF's job here is ranking, not filtering.
            print(f"    ℹ️  TF-IDF ranked {len(drosophila_papers)} papers")

        except Exception as e:
            print(f"    ⚠️  TF-IDF scoring failed ({e}), using unranked Drosophila papers")

        dropped_total = len(papers) - len(drosophila_papers)
        print(f"  🔍 Relevance filter: kept {len(drosophila_papers)}/{len(papers)} papers ({dropped_total} dropped)")

        return drosophila_papers

    # Gene-family nicknames that PubMed matches poorly ("FOX" also matches authors named Fox),
    # mapped to the term papers use.
    _GENE_FAMILY_TERMS = {
        'fox': 'forkhead',
        'hox': 'homeotic',
    }
    _CHAT_FILLER_WORDS = {
        'explain', 'describe', 'know', 'please', 'give', 'overview', 'show', 'find',
        'drosophila', 'melanogaster', 'family', 'families',
        'abt', 'pls', 'plz', 'u', 'ur', 'thx', 'wat', 'wht', 'hw', 'info',
    }

    def build_chat_search_query(self, message: str) -> str:
        """Turn a chat question into PubMed search terms ("tell me about FOX genes" -> "forkhead").

        PubMed ANDs every word, so conversational filler otherwise returns nothing.
        """
        family_terms, terms = [], []
        for word in re.findall(r"[A-Za-z0-9][\w\-]*", message):
            lower = word.lower()
            # Searches are already restricted to Drosophila; "dros..." also drops misspellings
            # like "drosophilia", which PubMed would otherwise require to match.
            if lower in _COMMON_WORDS or lower in self._CHAT_FILLER_WORDS or lower.startswith("dros"):
                continue
            family = self._GENE_FAMILY_TERMS.get(lower)
            if family:
                family_terms.append(family)
            else:
                terms.append(word)
        # Recognized family terms lead, so callers that drop trailing terms keep them longest.
        return " ".join(family_terms + terms) or message

    def search_pubmed(self, query, max_results=5):
        try:
            full_query = f"Drosophila AND ({query})"
            handle = Entrez.esearch(db="pubmed", term=full_query, retmax=max_results + 3, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            id_list = record.get("IdList", [])
            if not id_list:
                return []

            handle = Entrez.efetch(db="pubmed", id=id_list, rettype="abstract", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            papers = []
            for article in records.get('PubmedArticle', [])[:max_results]:
                medline = article['MedlineCitation']
                article_data = medline['Article']
                title = article_data.get('ArticleTitle', 'No title')
                abstract = article_data.get('Abstract', {}).get('AbstractText', [''])
                abstract = ' '.join(str(a) for a in abstract) if isinstance(abstract, list) else abstract
                authors = []
                for author in article_data.get('AuthorList', [])[:3]:
                    if 'LastName' in author:
                        authors.append(f"{author['LastName']} {author.get('Initials', '')}")
                pmid = str(medline['PMID'])
                year = _publication_year(article_data)
                papers.append({
                    'title': title,
                    'authors': ', '.join(authors) + ' et al.' if authors else 'Unknown',
                    'abstract': _clean_abstract(abstract),
                    'pmid': pmid,
                    'year': year,
                    'journal': article_data.get('Journal', {}).get('Title', ''),
                    'doi': _pubmed_doi(article),
                    'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    'source': 'PubMed'
                })

            return papers[:max_results]
        except Exception as e:
            print(f"  ❌ PubMed error: {e}")
            return []

    def extract_gene_names(self, text: str) -> List[str]:
        """Return up to 3 gene database keys mentioned in text, in order of appearance.

        Matches whole tokens only (never substrings, so "Drosophila" no longer yields "dr").
        Symbols that are common words count only when marked as a gene ("the how gene").
        Full gene names ("Notch", "Yorkie") count when capitalized, since many are plain words.
        """
        marked = set()
        for match in re.findall(r'(?:\bgenes?\s+([\w\-:()]+)|([\w\-:()]+)\s+genes?\b)', text, re.IGNORECASE):
            marked.add((match[0] or match[1]).lower())

        potential_genes = []
        for match in re.finditer(r"[A-Za-z0-9][\w\-:().'\[\]]*", text):
            token = match.group(0)
            for candidate in (token, token.rstrip(".,;:!?)]'\"")):
                lower = candidate.lower()
                if lower in _GENE_DB and len(lower) > 1:
                    if lower not in _COMMON_WORDS or lower in marked:
                        potential_genes.append(lower)
                    break
                if candidate[:1].isupper() and lower in _GENE_NAME_INDEX:
                    potential_genes.append(_GENE_NAME_INDEX[lower])
                    break

        unique = list(dict.fromkeys(potential_genes))
        return unique[:3]

    def format_papers(self, papers):
        if not papers:
            return ""
        formatted = "\n" + "="*70 + "\n"
        formatted += "RELEVANT PUBLICATIONS - CITE THESE IN YOUR RESPONSE\n"
        formatted += "="*70 + "\n\n"
        for i, paper in enumerate(papers, 1):
            token = paper_citation_token(paper, i)
            formatted += f"[{token}] {paper['authors']} ({paper['year']})\n"
            formatted += f"    Title: {paper['title']}\n"
            formatted += f"    URL: {paper['url']}\n"
            formatted += f"    Abstract: {paper['abstract']}\n"
            formatted += f"    Source: {paper.get('source', 'Unknown')}\n"
            if paper.get('journal'):
                formatted += f"    Journal: {paper['journal']}\n"
            if paper.get('doi'):
                formatted += f"    DOI: https://doi.org/{paper['doi']}\n"
            if isinstance(paper.get('_relevance_score'), (int, float)):
                formatted += (
                    "    TF-IDF text-match heuristic (not evidence quality): "
                    f"{paper['_relevance_score']:.2f}\n"
                )
            formatted += "\n"
        formatted += "="*70 + "\n"
        formatted += "INSTRUCTIONS: Cite a supported claim using the exact bracketed identifier shown above, e.g. [PMID:12345678].\n"
        formatted += "Do not create a References section; the application will add verified source metadata.\n"
        formatted += "="*70 + "\n\n"
        return formatted

    # Maps keyword sets → focused PubMed query (no "aging" unless the topic is aging)
    _TOPIC_QUERY_MAP = [
        (['lifespan', 'longevity', 'aging', 'ageing'],
            "Drosophila lifespan extension insulin signaling FOXO"),
        (['locomotor', 'climbing', 'movement', 'motor'],
            "Drosophila locomotor activity climbing assay neuromuscular"),
        (['fat body', 'adipose', 'metabolic', 'metabolism', 'lipid'],
            "Drosophila fat body metabolism lipid storage"),
        (['neuron', 'brain', 'neurodegeneration', 'cognitive', 'memory', 'neural'],
            "Drosophila neurodegeneration memory learning brain"),
        (['sleep', 'circadian', 'rhythm'],
            "Drosophila sleep circadian clock regulation"),
        (['gut', 'intestin', 'microbiome'],
            "Drosophila intestine gut stem cell homeostasis"),
        (['muscle', 'sarcopenia', 'flight'],
            "Drosophila muscle development flight sarcopenia"),
        (['stress', 'oxidative', 'ros', 'reactive oxygen'],
            "Drosophila oxidative stress ROS resistance"),
        (['autophagy', 'atg', 'lysosome'],
            "Drosophila autophagy ATG protein degradation"),
        (['tor', 'rapamycin', 'mtor', 'torc'],
            "Drosophila TOR signaling rapamycin nutrient sensing"),
        (['foxo', 'forkhead'],
            "Drosophila FOXO transcription factor stress"),
        (['inr', 'insulin', 'dilp', 'igf'],
            "Drosophila insulin signaling InR DILP"),
        (['notch', 'delta', 'serrate', 'jagged'],
            "Drosophila Notch signaling lateral inhibition"),
        (['wnt', 'wingless', 'armadillo', 'frizzled'],
            "Drosophila Wingless Wnt signaling development"),
        (['hedgehog', 'patched', 'smoothened'],
            "Drosophila Hedgehog signaling patterning"),
        (['jak', 'stat', 'hopscotch'],
            "Drosophila JAK STAT signaling immunity"),
        (['immunity', 'immune', 'toll', 'imd', 'infection', 'antimicrobial'],
            "Drosophila innate immunity Toll IMD antimicrobial"),
        (['apoptosis', 'cell death', 'reaper', 'caspase'],
            "Drosophila apoptosis cell death caspase"),
        (['stem cell', 'niche', 'progenitor'],
            "Drosophila stem cell niche self-renewal"),
        (['cancer', 'tumor', 'proliferation', 'neoplasm'],
            "Drosophila tumor cancer cell proliferation"),
        (['mitochondria', 'mitochondrial', 'electron transport', 'oxidative phosphorylation'],
            "Drosophila mitochondria electron transport chain"),
        (['epigenetic', 'chromatin', 'histone', 'methylation'],
            "Drosophila epigenetics chromatin remodeling"),
        (['development', 'patterning', 'morphogen'],
            "Drosophila development patterning body plan"),
    ]

    def _build_topic_queries(self, combined: str, topic: str) -> List[str]:
        def matches(kw: str) -> bool:
            return bool(re.search(r'\b' + re.escape(kw) + r'\b', combined))

        queries = [
            query for keywords, query in self._TOPIC_QUERY_MAP
            if any(matches(kw) for kw in keywords)
        ]
        if not queries:
            queries = [f"Drosophila {topic}"]
        return queries

    def fetch_literature_for_proposal(self, topic: str, genes: List[str], user_clarifications: str = "") -> List[Dict]:
        """
        Dedicated literature fetch for proposal generation.
        Runs multiple targeted PubMed searches to get 10-15 relevant papers.
        Much more thorough than the normal chat paper fetch.
        """
        print(f"  📚 Fetching dedicated proposal literature...")
        all_papers = []
        seen_pmids = set()

        def add_papers(new_papers):
            for p in new_papers:
                pmid = p.get('pmid', '')
                if pmid and pmid not in seen_pmids:
                    seen_pmids.add(pmid)
                    all_papers.append(p)

        # Search 1: Gene-specific FlyBase lookups
        for gene in genes[:3]:
            gene_info = search_flybase_fixed(gene)
            if gene_info:
                pubs = self.get_flybase_publications(gene_info['fbgn'], max_results=6)
                add_papers(pubs)
                print(f"    FlyBase {gene}: {len(pubs)} papers")

        # Search 2: Gene + topic keywords
        topic_keywords = " ".join(topic.split()[:4])
        for gene in genes[:2]:
            papers = self.search_pubmed(f"{gene} {topic_keywords}", max_results=5)
            add_papers(papers)
        if not genes:
            add_papers(self.search_pubmed(f"Drosophila {topic_keywords}", max_results=5))

        # Search 3: Topic-specific searches
        clean_clarifs = self.planner.strip_scope_words(user_clarifications) if user_clarifications else ""
        combined = f"{topic} {clean_clarifs}".lower()
        topic_queries = self._build_topic_queries(combined, topic)

        for query in topic_queries[:4]:  # limit to 4 topic searches
            papers = self.search_pubmed(query, max_results=4)
            add_papers(papers)
            print(f"    PubMed '{query[:50]}': {len(papers)} papers")

        # Search 4: Review papers scoped to the actual topic
        review_query = f"Drosophila {' '.join(topic.split()[:3])} review"
        review_papers = self.search_pubmed(review_query, max_results=3)
        add_papers(review_papers)

        # Filter for relevance before returning — drop non-Drosophila and off-topic papers
        filtered = self.filter_relevant_papers(all_papers, topic, genes, user_clarifications)

        print(f"  ✅ Papers after relevance filter: {len(filtered)}/{len(all_papers)}")
        return filtered[:15]  # cap at 15 to avoid token overflow

    def handle_planning_request(
        self,
        user_message: str,
        force: bool = False,
        auto_revise_claims: bool = True,
    ) -> Dict:
        """
        Handle a research planning request. Returns a dict with:
        - type: 'clarification' | 'proposal' | 'refinement'
        - content: the text to show the user
        - proposal: the structured proposal dict (if type == 'proposal')
        - ready_for_export: bool
        """
        print(f"\n{'='*70}\n🔬 PLANNING MODE: {user_message[:60]}\n{'='*70}")

        # If we already have a proposal, this is a refinement request
        if self.planner.current_proposal and not force:
            print("  ↻ Refining existing proposal...")
            updated, summary = self.planner.refine_proposal(
                user_message,
                self.conversation_history,
                auto_revise_claims=auto_revise_claims,
            )
            formatted = self.planner.format_proposal_for_chat(updated)
            return {
                'type': 'refinement',
                'content': f"✅ {summary}\n\n---\n\n{formatted}",
                'proposal': updated,
                'ready_for_export': True
            }

        # Extract topic — prefer last_topic if the message is a planning trigger phrase
        # rather than a real topic description (e.g. "plan research", "build a proposal")
        extracted_topic = self.planner.extract_topic_from_message(user_message)
        if len(extracted_topic.split()) <= 2 and self.last_topic:
            # Extracted topic is too short — likely a planning trigger phrase, use last_topic
            topic = self.last_topic
            print(f"  📌 Using cached topic: {topic[:60]}")
        else:
            topic = extracted_topic or self.last_topic or user_message[:100]

        genes = self.extract_gene_names(user_message)
        if not genes and self.last_genes:
            genes = self.last_genes

        papers = self.last_papers if self.last_papers else []

        # Fix 2: Preserve prior clarification answers across the thin-literature block
        # If we already have stored clarifications from a previous attempt, carry them forward
        stored_clarifications = self.planner.proposal_context.get('clarifications', '')

        # If we're awaiting clarification, use this message as the clarification and generate
        if self.awaiting_clarification:
            self.awaiting_clarification = False
            print("  📝 Received clarification — fetching dedicated literature...")

            proposal_topic = self.planner.proposal_context.get('topic', topic)
            proposal_genes = self.planner.proposal_context.get('genes', genes)

            # Run a thorough dedicated literature search using clarification context
            proposal_papers = self.fetch_literature_for_proposal(
                topic=proposal_topic,
                genes=proposal_genes,
                user_clarifications=user_message
            )

            # Merge with any cached papers, deduplicate by PMID
            cached = self.planner.proposal_context.get('papers', [])
            seen = {p.get('pmid') for p in proposal_papers if p.get('pmid')}
            for p in cached:
                if p.get('pmid') not in seen:
                    proposal_papers.append(p)
                    seen.add(p.get('pmid'))

            MIN_PAPERS = 4
            if len(proposal_papers) < MIN_PAPERS:
                # Not enough literature — warn user but SAVE their clarifications
                # so they don't have to answer questions again after chatting more
                found = len(proposal_papers)
                self.planner.set_context_from_chat(
                    proposal_topic, proposal_genes, proposal_papers,
                    clarifications=user_message  # store their answers
                )
                return {
                    'type': 'clarification',
                    'content': (
                        f"⚠️ **Literature too thin to generate a strong proposal.**\n\n"
                        f"I searched PubMed and FlyBase but only found **{found} paper(s)** "
                        f"on *{proposal_topic}*. A well-cited proposal needs at least {MIN_PAPERS}.\n\n"
                        f"**Your answers have been saved** — you won't need to answer the questions again.\n\n"
                        f"A few options:\n"
                        f"- **Chat with me first** (e.g. *'What is known about FOXO in fat body aging?'*) "
                        f"to pull more papers, then click **Plan Research** again\n"
                        f"- **Share specific papers** you know — paste a PMID or title\n"
                        f"- **Broaden your topic** slightly if it's very niche\n\n"
                        f"What would you like to do?"
                    ),
                    'proposal': None,
                    'ready_for_export': False
                }

            print(f"  ✅ Generating proposal with {len(proposal_papers)} papers...")
            # Strip scope/admin words from clarifications before using for topic/literature
            clean_clarifications = self.planner.strip_scope_words(user_message)
            print(f"  🧹 Cleaned clarifications: '{clean_clarifications[:80]}'")

            # Update the planner context with the enriched paper set + clarifications
            self.planner.set_context_from_chat(proposal_topic, proposal_genes, proposal_papers, clarifications=clean_clarifications)

            proposal = self.planner.generate_proposal(
                topic=proposal_topic,
                genes=proposal_genes,
                papers=proposal_papers,
                user_clarifications=clean_clarifications,
                conversation_history=self.conversation_history
            )
            formatted = self.planner.format_proposal_for_chat(proposal)
            return {
                'type': 'proposal',
                'content': formatted,
                'proposal': proposal,
                'ready_for_export': True
            }

        # If we already have stored clarifications from a previous attempt
        # (e.g. user answered questions, hit thin-literature block, then chatted more)
        # skip asking again and go straight to generation with the stored answers
        # Strip scope words from stored clarifications before reuse
        if stored_clarifications:
            stored_clarifications = self.planner.strip_scope_words(stored_clarifications)
        if stored_clarifications and len(stored_clarifications.split()) > 3:
            print(f"  ♻️  Reusing stored clarifications: {stored_clarifications[:60]}...")
            self.planner.set_context_from_chat(topic, genes, papers, clarifications=stored_clarifications)
            proposal_papers = self.fetch_literature_for_proposal(
                topic=topic,
                genes=genes,
                user_clarifications=stored_clarifications
            )
            cached = papers
            seen = {p.get('pmid') for p in proposal_papers if p.get('pmid')}
            for p in cached:
                if p.get('pmid') not in seen:
                    proposal_papers.append(p)
                    seen.add(p.get('pmid'))

            MIN_PAPERS = 4
            if len(proposal_papers) < MIN_PAPERS:
                return {
                    'type': 'clarification',
                    'content': (
                        f"⚠️ **Still finding limited literature** ({len(proposal_papers)} papers found).\n\n"
                        f"Try asking me a more specific question about your topic first, "
                        f"then click Plan Research again."
                    ),
                    'proposal': None,
                    'ready_for_export': False
                }

            self.planner.set_context_from_chat(topic, genes, proposal_papers, clarifications=stored_clarifications)
            proposal = self.planner.generate_proposal(
                topic=topic,
                genes=genes,
                papers=proposal_papers,
                user_clarifications=stored_clarifications,
                conversation_history=self.conversation_history,
                auto_revise_claims=auto_revise_claims,
            )
            formatted = self.planner.format_proposal_for_chat(proposal)
            return {
                'type': 'proposal',
                'content': formatted,
                'proposal': proposal,
                'ready_for_export': True
            }

        # No prior clarifications — check if existing context is sufficient
        print("  ❓ Checking if clarifying questions are needed...")
        question = self.planner.generate_clarifying_questions(topic, genes, self.conversation_history)

        if question is None:
            print("  ✅ Sufficient context — skipping clarification")
            self.planner.set_context_from_chat(topic, genes, papers)
            proposal_papers = self.fetch_literature_for_proposal(
                topic=topic, genes=genes, user_clarifications=topic
            )
            self.planner.set_context_from_chat(topic, genes, proposal_papers)
            proposal = self.planner.generate_proposal(
                topic=topic, genes=genes, papers=proposal_papers,
                conversation_history=self.conversation_history,
                auto_revise_claims=auto_revise_claims,
            )
            formatted = self.planner.format_proposal_for_chat(proposal)
            return {
                'type': 'proposal',
                'content': formatted,
                'proposal': proposal,
                'ready_for_export': True
            }

        self.planner.set_context_from_chat(topic, genes, papers)
        self.awaiting_clarification = True
        return {
            'type': 'clarification',
            'content': f"🔬 **Research Planner activated!**\n\nBefore I build your proposal, {question}",
            'proposal': None,
            'ready_for_export': False
        }

    def chat(
        self,
        user_message: str,
        force_planning: bool = False,
        auto_revise_claims: bool = True,
    ) -> Dict:
        """
        Main chat function. Returns dict with:
        - response: text to display
        - is_planning: bool
        - proposal: dict or None
        - ready_for_export: bool
        """
        print(f"\n{'='*70}\nQuery: {user_message[:80]}\n{'='*70}")

        # Check for experiment design intent — redirect to panel buttons
        msg_lower = user_message.lower()
        experiment_design_phrases = [
            'do the experiment', 'design the experiment', 'run aim', 'do aim',
            'experiment for aim', 'design aim', 'protocol for aim', 'how to run aim'
        ]
        if any(phrase in msg_lower for phrase in experiment_design_phrases) and self.planner.current_proposal:
            aims = self.planner.current_proposal.get('specific_aims', [])
            aim_list = ', '.join([f"Aim {a.get('aim_number')}: {a.get('title', '')}" for a in aims])
            return {
                'response': f"💡 To design a bench-ready protocol for a specific aim, click the aim buttons in the **Research Proposal panel** on the right.\n\nYour current proposal has:\n{aim_list}\n\nClick **Aim 1**, **Aim 2**, or **Aim 3** in the panel to generate a full experiment design including cross scheme, timepoints, data collection sheet, and go/no-go criteria.",
                'is_planning': False,
                'proposal': None,
                'ready_for_export': False
            }

        # Check for experiment design intent first (most specific)
        is_experiment_design, aim_number = detect_experiment_design_intent(user_message)
        if is_experiment_design and self.planner.current_proposal:
            print(f"  🧪 Experiment design intent detected for Aim {aim_number}")
            aim = next((a for a in self.planner.current_proposal.get('specific_aims', [])
                       if a.get('aim_number') == aim_number), None)
            if aim:
                design = self.planner.generate_experiment_design(aim, self.planner.current_proposal)
                formatted = self.planner.format_experiment_design(design)
                self.conversation_history.append({"role": "user", "content": user_message})
                self.conversation_history.append({"role": "assistant", "content": formatted})
                return {
                    'response': formatted,
                    'is_planning': True,
                    'proposal': None,
                    'ready_for_export': False,
                    'experiment_design': True,
                    'aim_number': aim_number
                }
            else:
                # Aim not found — tell user which aims exist
                aims = [a.get('aim_number') for a in self.planner.current_proposal.get('specific_aims', [])]
                return {
                    'response': f"I don't see Aim {aim_number} in the current proposal. Available aims: {aims}. Click the aim button in the panel or ask for one of those.",
                    'is_planning': False,
                    'proposal': None,
                    'ready_for_export': False
                }

        # Check for planning intent
        has_prior_context = bool(self.last_topic or self.last_papers)
        is_planning, confidence = detect_planning_intent(user_message, has_prior_context)

        # If awaiting clarification, route to planner if message looks like
        # an answer (short, descriptive) rather than a new research question.
        # NOTE: do NOT require self.planning_mode — it isn't set yet on the turn
        # when the user answers clarifying questions (it gets set after the call).
        is_clarification_answer = (
            self.awaiting_clarification
            and not is_planning  # not a new explicit planning request
            and len(user_message.split()) < 80  # reasonably short answer
            and not user_message.strip().startswith('?')  # doesn't start as a question
        )
        print(f"  🔍 Routing: awaiting={self.awaiting_clarification}, is_planning={is_planning}, force={force_planning}, is_clarification_answer={is_clarification_answer}, words={len(user_message.split())}")

        if force_planning or is_planning or is_clarification_answer:
            self.planning_mode = True
            result = self.handle_planning_request(
                user_message,
                force=force_planning,
                auto_revise_claims=auto_revise_claims,
            )

            # Add to history
            self.conversation_history.append({"role": "user", "content": user_message})
            self.conversation_history.append({"role": "assistant", "content": result['content']})
            if len(self.conversation_history) > 12:
                self.conversation_history = self.conversation_history[-12:]

            return {
                'response': result['content'],
                'is_planning': True,
                'proposal': result.get('proposal'),
                'ready_for_export': result.get('ready_for_export', False)
            }

        # ── Normal chat mode ──────────────────────────────────────────────────
        self.planning_mode = False
        # Only reset awaiting_clarification if this is clearly a new unrelated question
        # (i.e. not a short answer that just failed the is_clarification_answer check)
        if not self.awaiting_clarification:
            pass  # already false, nothing to reset
        else:
            # Keep awaiting_clarification=True only if message is very short
            # (user typed something ambiguous). Reset if it's a full research question.
            if len(user_message.split()) > 15:
                self.awaiting_clarification = False
        all_papers = []

        genes = self.extract_gene_names(user_message)
        if genes:
            print(f"  Found genes: {genes}")
            self.last_genes = genes
            for gene in genes:
                gene_info = search_flybase_fixed(gene)
                if gene_info:
                    num_papers = self.determine_paper_count(user_message)
                    pubs = self.get_flybase_publications(gene_info['fbgn'], max_results=num_papers)
                    if pubs:
                        all_papers.extend(pubs)

        num_needed = self.determine_paper_count(user_message)
        if len(all_papers) < num_needed:
            print(f"  Supplementing with PubMed...")
            # PubMed ANDs every term, so one typo or rare word narrows the results to a few
            # stray matches; drop trailing terms and keep searching until enough papers match.
            terms = self.build_chat_search_query(user_message).split()
            pubmed_papers = []
            while terms and len(pubmed_papers) < num_needed - len(all_papers):
                seen = {p.get('pmid') for p in all_papers + pubmed_papers}
                pubmed_papers += [
                    p for p in self.search_pubmed(
                        " ".join(terms), max_results=num_needed - len(all_papers) - len(pubmed_papers)
                    )
                    if p.get('pmid') not in seen
                ]
                terms = terms[:-1]
            if pubmed_papers:
                all_papers.extend(pubmed_papers)

        # Filter chat papers for relevance too
        if all_papers:
            all_papers = self.filter_relevant_papers(all_papers, user_message, genes)

        # Cache for planner
        if all_papers:
            self.last_papers = all_papers
        # Only update last_topic if this isn't a planning trigger phrase
        _is_plan, _ = detect_planning_intent(user_message, False)
        if not _is_plan and len(user_message.split()) > 3:
            self.last_topic = user_message[:100]

        publication_context = self.format_papers(all_papers) if all_papers else ""
        enhanced_message = user_message
        if publication_context:
            enhanced_message += "\n\n" + publication_context

        self.conversation_history.append({"role": "user", "content": enhanced_message})
        if len(self.conversation_history) > 10:
            self.conversation_history = self.conversation_history[-10:]

        system_prompt = """You are a specialized AI assistant for Drosophila melanogaster (fruit fly) research.

Your expertise includes:
- Drosophila genetics, development, and molecular biology
- Gene nomenclature and function
- Experimental techniques in fly research
- Developmental biology and signaling pathways

EVIDENCE AND UNCERTAINTY RULES:
- Treat the current message's RELEVANT PUBLICATIONS block as the only retrieved evidence for this answer. Prior chat text is not a source.
- Put the exact bracketed source identifier shown in the current publications block immediately after each literature-supported claim (for example, [PMID:12345678]). These identifiers are unique to the current retrieved records and avoid confusing citations from earlier turns.
- Never invent papers, PMID/DOI values, author names, results, or citation numbers. The application will add source details from retrieved records; do not write a References section yourself.
- If a point is an inference rather than directly reported evidence, label it as an inference and cite the evidence it builds on. Clearly label proposed experiments and predictions as proposals, not established findings.
- If a claim is not supported by the supplied title or abstract, say the evidence is insufficient or describe it as unverified instead of stating it as fact.
- If no publications were retrieved, state that no relevant papers were found. You may offer brief general background only if it is explicitly labeled unverified and do not make paper-specific claims.
- Do not force citations onto every sentence or claim multiple papers support unless the evidence actually warrants it.
- Return a JSON object with "answer" (the user-facing prose) and "evidence_claims" (an array of material claims). Each claim object must contain an exact contiguous quote from answer, evidence_type (finding, inference, general_background, proposal, or prediction), and source_ids (the exact PMID/SOURCE identifiers cited by that claim, or []).
- Include every material factual claim in evidence_claims, including unsupported claims; do not omit one because it lacks a citation. Do not include conversational transitions or advice.

When discussing genes or research:
- Emphasize genetic pathways and mechanisms
- Mention human orthologues and conservation
- Discuss functional relationships and interactions
- Always cite relevant papers provided

At the end of responses, if the topic could support a research proposal, add a brief note:
"💡 *Want to turn this into a research proposal? Click **Plan Research** or ask me to build a proposal.*"
"""

        print("🤖 Calling Claude...")
        response = self.client.messages.create(
            model="claude-sonnet-4-6",
            max_tokens=6000,
            system=system_prompt,
            messages=self.conversation_history
        )

        assistant_message = response.content[0].text
        assistant_message, structured_claims = parse_structured_chat_response(
            assistant_message
        )
        response_content, evidence_review = verify_and_revise_claims(
            assistant_message,
            all_papers,
            self.client,
            structured_claims=structured_claims,
            content={
                "answer": assistant_message,
                "evidence_claims": structured_claims or [],
            },
            auto_revise=auto_revise_claims,
        )
        assistant_message = response_content["answer"]
        assistant_message = validate_answer_citations(assistant_message, all_papers)
        assistant_message += "\n\n" + format_evidence_review(evidence_review)
        self.conversation_history.append({"role": "assistant", "content": assistant_message})

        print("✅ Response generated\n" + "="*70 + "\n")
        return {
            'response': assistant_message,
            'is_planning': False,
            'proposal': None,
            'ready_for_export': False
        }

    def reset_conversation(self):
        self.conversation_history = []
        self.mentioned_genes = set()
        self.last_topic = ""
        self.last_papers = []
        self.last_genes = []
        self.planning_mode = False
        self.awaiting_clarification = False
        self.planner = ResearchPlanner(self.client)