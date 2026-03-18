import anthropic
import os
from Bio import Entrez
import json
import requests
from typing import Optional, Dict, List
import re
from research_planner import ResearchPlanner, detect_planning_intent

# Configure Entrez
Entrez.email = "jordanszt@icloud.com"

# KNOWN GENES DATABASE - 73 common Drosophila genes
KNOWN_GENES = {
    # Insulin/TOR pathway
    'inr': ('InR', 'Insulin-like receptor', 'FBgn0283499'),
    'insulin receptor': ('InR', 'Insulin-like receptor', 'FBgn0283499'),
    'foxo': ('foxo', 'forkhead box, sub-group O', 'FBgn0038197'),
    'dfoxo': ('foxo', 'forkhead box, sub-group O', 'FBgn0038197'),
    'tor': ('Tor', 'Target of rapamycin', 'FBgn0021796'),
    'pten': ('Pten', 'Pten', 'FBgn0026379'),
    'akt': ('Akt1', 'Akt1', 'FBgn0010379'),
    'akt1': ('Akt1', 'Akt1', 'FBgn0010379'),
    'pi3k': ('Pi3K92E', 'Pi3 kinase 92E', 'FBgn0015279'),
    'dilp2': ('Dilp2', 'Drosophila insulin-like peptide 2', 'FBgn0036046'),
    # Notch pathway
    'n': ('N', 'Notch', 'FBgn0004647'),
    'notch': ('N', 'Notch', 'FBgn0004647'),
    'dl': ('Dl', 'Delta', 'FBgn0000463'),
    'delta': ('Dl', 'Delta', 'FBgn0000463'),
    'ser': ('Ser', 'Serrate', 'FBgn0004197'),
    'serrate': ('Ser', 'Serrate', 'FBgn0004197'),
    # Hedgehog pathway
    'hh': ('hh', 'hedgehog', 'FBgn0004644'),
    'hedgehog': ('hh', 'hedgehog', 'FBgn0004644'),
    'ptc': ('ptc', 'patched', 'FBgn0003892'),
    'patched': ('ptc', 'patched', 'FBgn0003892'),
    'smo': ('smo', 'smoothened', 'FBgn0003444'),
    'smoothened': ('smo', 'smoothened', 'FBgn0003444'),
    'ci': ('ci', 'cubitus interruptus', 'FBgn0004859'),
    # Wnt/Wingless pathway
    'wg': ('wg', 'wingless', 'FBgn0004009'),
    'wingless': ('wg', 'wingless', 'FBgn0004009'),
    'arm': ('arm', 'armadillo', 'FBgn0000117'),
    'armadillo': ('arm', 'armadillo', 'FBgn0000117'),
    'fz': ('fz', 'frizzled', 'FBgn0001085'),
    'frizzled': ('fz', 'frizzled', 'FBgn0001085'),
    # EGFR pathway
    'egfr': ('Egfr', 'Epidermal growth factor receptor', 'FBgn0003731'),
    'ras': ('Ras85D', 'Ras oncogene at 85D', 'FBgn0003204'),
    'ras85d': ('Ras85D', 'Ras oncogene at 85D', 'FBgn0003204'),
    # JAK/STAT pathway
    'hop': ('hop', 'hopscotch', 'FBgn0004864'),
    'hopscotch': ('hop', 'hopscotch', 'FBgn0004864'),
    'stat': ('Stat92E', 'Signal-transducer and activator of transcription protein at 92E', 'FBgn0016917'),
    'stat92e': ('Stat92E', 'Signal-transducer and activator of transcription protein at 92E', 'FBgn0016917'),
    # Segmentation genes
    'eve': ('eve', 'even-skipped', 'FBgn0000606'),
    'even-skipped': ('eve', 'even-skipped', 'FBgn0000606'),
    'ftz': ('ftz', 'fushi tarazu', 'FBgn0001078'),
    'en': ('en', 'engrailed', 'FBgn0000577'),
    'engrailed': ('en', 'engrailed', 'FBgn0000577'),
    'hb': ('hb', 'hunchback', 'FBgn0001180'),
    'hunchback': ('hb', 'hunchback', 'FBgn0001180'),
    'kr': ('Kr', 'Kruppel', 'FBgn0001325'),
    'kruppel': ('Kr', 'Kruppel', 'FBgn0001325'),
    # Maternal genes
    'bcd': ('bcd', 'bicoid', 'FBgn0000166'),
    'bicoid': ('bcd', 'bicoid', 'FBgn0000166'),
    'nos': ('nos', 'nanos', 'FBgn0002962'),
    'nanos': ('nos', 'nanos', 'FBgn0002962'),
    'osk': ('osk', 'oskar', 'FBgn0003028'),
    'oskar': ('osk', 'oskar', 'FBgn0003028'),
    # Aging/longevity
    'mth': ('mth', 'methuselah', 'FBgn0002774'),
    'methuselah': ('mth', 'methuselah', 'FBgn0002774'),
    'sir2': ('Sir2', 'Sir2', 'FBgn0010309'),
    'sod': ('Sod1', 'Superoxide dismutase 1', 'FBgn0003462'),
    'sod1': ('Sod1', 'Superoxide dismutase 1', 'FBgn0003462'),
    'cat': ('Cat', 'Catalase', 'FBgn0000251'),
    'catalase': ('Cat', 'Catalase', 'FBgn0000251'),
    # Apoptosis
    'p53': ('p53', 'p53', 'FBgn0039044'),
    'dmp53': ('p53', 'p53', 'FBgn0039044'),
    'rpr': ('rpr', 'reaper', 'FBgn0011706'),
    'reaper': ('rpr', 'reaper', 'FBgn0011706'),
    'hid': ('hid', 'head involution defective', 'FBgn0003997'),
    # Circadian / Sleep
    'per': ('per', 'period', 'FBgn0003068'),
    'period': ('per', 'period', 'FBgn0003068'),
    'tim': ('tim', 'timeless', 'FBgn0014396'),
    'timeless': ('tim', 'timeless', 'FBgn0014396'),
    'clk': ('Clk', 'Clock', 'FBgn0023076'),
    'clock': ('Clk', 'Clock', 'FBgn0023076'),
    'cyc': ('cyc', 'cycle', 'FBgn0023094'),
    'cycle': ('cyc', 'cycle', 'FBgn0023094'),
    'sss': ('sss', 'sleepless', 'FBgn0038528'),
    'sleepless': ('sss', 'sleepless', 'FBgn0038528'),
    'shaker': ('Sh', 'Shaker', 'FBgn0003380'),
    'sh': ('Sh', 'Shaker', 'FBgn0003380'),
    # Classic markers
    'w': ('w', 'white', 'FBgn0003996'),
    'white': ('w', 'white', 'FBgn0003996'),
    'y': ('y', 'yellow', 'FBgn0004034'),
    'yellow': ('y', 'yellow', 'FBgn0004034'),
}


def search_flybase_fixed(gene_name: str) -> Optional[Dict]:
    """Optimized FlyBase search: known genes first, then quick lookup fallback."""
    try:
        gene_name = gene_name.strip()
        gene_lower = gene_name.lower()

        if gene_lower in KNOWN_GENES:
            symbol, name, fbgn = KNOWN_GENES[gene_lower]
            print(f"  ✅ Found: {symbol} ({fbgn}) [Database]")
            return {
                'symbol': symbol,
                'name': name,
                'fbgn': fbgn,
                'summary': f"Gene symbol: {symbol}.",
                'synonyms': [],
                'url': f"https://flybase.org/reports/{fbgn}"
            }

        try:
            session = requests.Session()
            session.headers.update({
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            })
            response = session.get(
                f"https://flybase.org/search?query={gene_name}",
                timeout=8,
                allow_redirects=True
            )
            fbgn_match = re.search(r'FBgn\d{7}', response.text)
            if fbgn_match:
                fbgn = fbgn_match.group(0)
                print(f"  ✅ Found: {gene_name} ({fbgn}) [Lookup]")
                return {
                    'symbol': gene_name,
                    'name': f"{gene_name}",
                    'fbgn': fbgn,
                    'summary': f"Gene symbol: {gene_name}.",
                    'synonyms': [],
                    'url': f"https://flybase.org/reports/{fbgn}"
                }
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError, requests.exceptions.ProxyError):
            pass

        print(f"  ℹ️  '{gene_name}' not found")
        return None

    except Exception as e:
        print(f"  ⚠️  Error: {e}")
        return None


class DrosophilaAssistant:
    def __init__(self, api_key):
        self.client = anthropic.Anthropic(api_key=api_key)
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
        try:
            print(f"  📚 Fetching publications for {fbgn}...")
            url = f"https://flybase.org/reports/{fbgn}"
            try:
                session = requests.Session()
                session.headers.update({
                    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
                })
                max_retries = 3
                for attempt in range(max_retries):
                    response = session.get(url, timeout=8, allow_redirects=True)
                    if response.status_code == 200:
                        break
                    elif response.status_code == 202:
                        if attempt < max_retries - 1:
                            import time
                            wait_time = 2 ** attempt
                            print(f"    ⏳ FlyBase processing (202), waiting {wait_time}s...")
                            time.sleep(wait_time)
                        else:
                            return self._fallback_pubmed_search(fbgn, max_results)
                    else:
                        return self._fallback_pubmed_search(fbgn, max_results)
            except (requests.exceptions.Timeout, requests.exceptions.ConnectionError, requests.exceptions.ProxyError):
                return self._fallback_pubmed_search(fbgn, max_results)

            pmid_pattern = r'/pubmed/(\d+)'
            pmids = re.findall(pmid_pattern, response.text, re.IGNORECASE)
            pmids = list(dict.fromkeys(pmids))[:max_results]

            if not pmids:
                return self._fallback_pubmed_search(fbgn, max_results)

            try:
                handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="xml")
                records = Entrez.read(handle)
                handle.close()

                publications = []
                for i, article in enumerate(records.get('PubmedArticle', [])):
                    try:
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
                        year = 'N/A'
                        if 'ArticleDate' in article_data and article_data['ArticleDate']:
                            year = article_data['ArticleDate'][0].get('Year', 'N/A')
                        publications.append({
                            'title': title,
                            'authors': ', '.join(authors) + ' et al.' if authors else 'Unknown',
                            'year': year,
                            'pmid': pmid,
                            'abstract': abstract[:500] + '...' if len(abstract) > 500 else abstract,
                            'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                            'source': 'FlyBase',
                            'flybase_rank': i + 1
                        })
                    except Exception as e:
                        continue

                if publications:
                    print(f"  ✅ Retrieved {len(publications)} publications from FlyBase")
                return publications

            except Exception as e:
                return [{'pmid': pmid, 'title': 'Publication', 'authors': 'See PubMed', 'year': 'N/A',
                         'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", 'source': 'FlyBase'} for pmid in pmids[:max_results]]

        except Exception as e:
            return self._fallback_pubmed_search(fbgn, max_results)

    def _fallback_pubmed_search(self, fbgn: str, max_results: int) -> List[Dict]:
        try:
            gene_name = fbgn
            for key, (symbol, name, fbn) in KNOWN_GENES.items():
                if fbn == fbgn:
                    gene_name = symbol
                    break

            full_query = f"Drosophila AND {gene_name}"
            handle = Entrez.esearch(db="pubmed", term=full_query, retmax=max_results + 2, sort="relevance")
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
                    abstract = ' '.join(str(a) for a in abstract) if isinstance(abstract, list) else abstract
                    authors = []
                    for author in article_data.get('AuthorList', [])[:3]:
                        if 'LastName' in author:
                            authors.append(f"{author['LastName']} {author.get('Initials', '')}")
                    pmid = str(medline['PMID'])
                    year = 'N/A'
                    if 'ArticleDate' in article_data and article_data['ArticleDate']:
                        year = article_data['ArticleDate'][0].get('Year', 'N/A')
                    publications.append({
                        'title': title,
                        'authors': ', '.join(authors) + ' et al.' if authors else 'Unknown',
                        'year': year,
                        'pmid': pmid,
                        'abstract': abstract[:500] + '...' if len(abstract) > 500 else abstract,
                        'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        'source': 'FlyBase (via PubMed)'
                    })
                except Exception:
                    continue

            return publications
        except Exception as e:
            return []

    def filter_relevant_papers(self, papers: List[Dict], topic: str, genes: List[str], user_clarifications: str = "") -> List[Dict]:
        """
        Filter papers for relevance before passing to the planner.
        Drops papers that don't mention Drosophila or the core topic.
        Two-pass: hard filter (must mention Drosophila), then soft score.
        """
        if not papers:
            return []

        combined_context = f"{topic} {user_clarifications}".lower()

        # Build keyword sets
        drosophila_terms = {'drosophila', 'melanogaster', 'fruit fly', 'flies', 'drosophilid'}

        # Topic keywords extracted from context
        topic_keywords = set()
        keyword_map = {
            'aging': ['aging', 'ageing', 'lifespan', 'longevity', 'senescence', 'gerontology'],
            'autophagy': ['autophagy', 'atg', 'lysosome', 'autophagosome'],
            'intestin': ['intestin', 'gut', 'midgut', 'epithelium', 'stem cell'],
            'foxo': ['foxo', 'forkhead', 'transcription factor'],
            'insulin': ['insulin', 'igf', 'inr', 'dilp', 'irs'],
            'tor': ['tor', 'rapamycin', 'mtor', 's6k'],
            'locomotor': ['locomotor', 'climbing', 'motor', 'movement'],
            'sleep': ['sleep', 'circadian', 'rhythm'],
            'neuron': ['neuron', 'neural', 'brain', 'cognitive'],
            'fat body': ['fat body', 'adipose', 'lipid'],
            'muscle': ['muscle', 'sarcopenia', 'myopathy'],
            'stress': ['oxidative stress', 'ros', 'paraquat', 'hydrogen peroxide'],
        }

        for key, terms in keyword_map.items():
            if any(k in combined_context for k in [key] + terms[:2]):
                topic_keywords.update(terms)

        gene_keywords = {g.lower() for g in genes}

        filtered = []
        for paper in papers:
            title = (paper.get('title', '') or '').lower()
            abstract = (paper.get('abstract', '') or '').lower()
            combined = title + ' ' + abstract

            # Hard filter: must mention Drosophila or fly model
            has_drosophila = any(term in combined for term in drosophila_terms)
            if not has_drosophila:
                print(f"    ❌ Dropped (no Drosophila): {paper.get('title', '')[:60]}")
                continue

            # Soft score: how many topic keywords appear
            topic_score = sum(1 for kw in topic_keywords if kw in combined)
            gene_score = sum(2 for g in gene_keywords if g in combined)  # genes weighted higher
            total_score = topic_score + gene_score

            paper['_relevance_score'] = total_score
            filtered.append(paper)

        # Sort by relevance score descending
        filtered.sort(key=lambda p: p.get('_relevance_score', 0), reverse=True)

        dropped = len(papers) - len(filtered)
        if dropped > 0:
            print(f"  🔍 Relevance filter: kept {len(filtered)}/{len(papers)} papers ({dropped} dropped)")

        return filtered

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
                year = 'N/A'
                if 'ArticleDate' in article_data and article_data['ArticleDate']:
                    year = article_data['ArticleDate'][0].get('Year', 'N/A')
                papers.append({
                    'title': title,
                    'authors': ', '.join(authors) + ' et al.' if authors else 'Unknown',
                    'abstract': abstract[:500] + '...' if len(abstract) > 500 else abstract,
                    'pmid': pmid,
                    'year': year,
                    'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    'source': 'PubMed'
                })

            return papers[:max_results]
        except Exception as e:
            print(f"  ❌ PubMed error: {e}")
            return []

    def extract_gene_names(self, text: str) -> List[str]:
        potential_genes = []
        exclude = {'aging', 'gene', 'genes', 'protein', 'proteins', 'what', 'is', 'are', 'tell',
                   'me', 'about', 'a', 'n', 'the', 'of', 'in', 'sleep', 'memory', 'behavior'}
        text_lower = text.lower()

        for gene_key in KNOWN_GENES.keys():
            if len(gene_key) > 1 and gene_key in text_lower and gene_key not in exclude:
                potential_genes.append(gene_key)

        pattern = r'(?:gene\s+([a-zA-Z]{2,})|([a-zA-Z]{2,})\s+gene)'
        matches = re.findall(pattern, text, re.IGNORECASE)
        for match in matches:
            gene = match[0] or match[1]
            if gene.lower() not in exclude and len(gene) > 1:
                potential_genes.append(gene)

        seen = set()
        unique = []
        for g in potential_genes:
            g_lower = g.lower()
            if g_lower not in seen and g_lower not in exclude:
                unique.append(g_lower)
                seen.add(g_lower)

        return unique[:3]

    def format_papers(self, papers):
        if not papers:
            return ""
        formatted = "\n" + "="*70 + "\n"
        formatted += "RELEVANT PUBLICATIONS - CITE THESE IN YOUR RESPONSE\n"
        formatted += "="*70 + "\n\n"
        for i, paper in enumerate(papers, 1):
            formatted += f"[{i}] {paper['authors']} ({paper['year']})\n"
            formatted += f"    Title: {paper['title']}\n"
            formatted += f"    URL: {paper['url']}\n"
            formatted += f"    Abstract: {paper['abstract']}\n"
            formatted += f"    Source: {paper.get('source', 'Unknown')}\n\n"
        formatted += "="*70 + "\n"
        formatted += "INSTRUCTIONS: Reference these papers in your answer using [1], [2], etc.\n"
        formatted += "Include the URL links and create a References section at the end\n"
        formatted += "="*70 + "\n\n"
        return formatted

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

        # Search 2: Gene + aging
        for gene in genes[:2]:
            papers = self.search_pubmed(f"{gene} aging lifespan longevity", max_results=5)
            add_papers(papers)

        # Search 3: Topic-specific searches
        topic_queries = []

        # Extract key biological terms from topic + clarifications
        combined = f"{topic} {user_clarifications}".lower()

        if any(w in combined for w in ['lifespan', 'longevity', 'aging', 'ageing']):
            topic_queries.append("Drosophila lifespan extension insulin signaling FOXO")
        if any(w in combined for w in ['locomotor', 'climbing', 'movement', 'motor']):
            topic_queries.append("Drosophila locomotor aging decline")
        if any(w in combined for w in ['fat body', 'adipose', 'metabolic']):
            topic_queries.append("Drosophila fat body aging metabolism")
        if any(w in combined for w in ['neuron', 'brain', 'cognitive', 'memory', 'neural']):
            topic_queries.append("Drosophila neuronal aging brain")
        if any(w in combined for w in ['sleep']):
            topic_queries.append("Drosophila sleep aging circadian")
        if any(w in combined for w in ['gut', 'intestin']):
            topic_queries.append("Drosophila gut aging intestine stem cell")
        if any(w in combined for w in ['muscle', 'sarcopenia']):
            topic_queries.append("Drosophila muscle aging sarcopenia")
        if any(w in combined for w in ['stress', 'oxidative', 'ros']):
            topic_queries.append("Drosophila oxidative stress resistance aging")
        if any(w in combined for w in ['autophagy', 'atg']):
            topic_queries.append("Drosophila autophagy aging lifespan")
        if any(w in combined for w in ['tor', 'rapamycin', 'mtor']):
            topic_queries.append("Drosophila TOR signaling aging rapamycin")
        if any(w in combined for w in ['foxo', 'forkhead']):
            topic_queries.append("Drosophila FOXO transcription factor lifespan")
        if any(w in combined for w in ['inr', 'insulin', 'dilp', 'igf']):
            topic_queries.append("Drosophila insulin signaling InR aging lifespan")

        # Default broad search if no specific matches
        if not topic_queries:
            topic_queries = [f"Drosophila {topic} aging"]

        for query in topic_queries[:4]:  # limit to 4 topic searches
            papers = self.search_pubmed(query, max_results=4)
            add_papers(papers)
            print(f"    PubMed '{query[:50]}': {len(papers)} papers")

        # Search 4: Key review papers — always useful for proposals
        review_papers = self.search_pubmed("Drosophila aging review mechanisms longevity", max_results=3)
        add_papers(review_papers)

        # Filter for relevance before returning — drop non-Drosophila and off-topic papers
        filtered = self.filter_relevant_papers(all_papers, topic, genes, user_clarifications)

        print(f"  ✅ Papers after relevance filter: {len(filtered)}/{len(all_papers)}")
        return filtered[:15]  # cap at 15 to avoid token overflow

    def handle_planning_request(self, user_message: str, force: bool = False) -> Dict:
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
                self.conversation_history
            )
            formatted = self.planner.format_proposal_for_chat(updated)
            return {
                'type': 'refinement',
                'content': f"✅ {summary}\n\n---\n\n{formatted}",
                'proposal': updated,
                'ready_for_export': True
            }

        # Extract topic
        topic = self.planner.extract_topic_from_message(user_message)
        if not topic and self.last_topic:
            topic = self.last_topic

        genes = self.extract_gene_names(user_message)
        if not genes and self.last_genes:
            genes = self.last_genes

        papers = self.last_papers if self.last_papers else []

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
                # Not enough literature — warn user instead of generating a thin proposal
                found = len(proposal_papers)
                return {
                    'type': 'clarification',
                    'content': (
                        f"⚠️ **Literature too thin to generate a strong proposal.**\n\n"
                        f"I searched PubMed and FlyBase but only found **{found} paper(s)** "
                        f"relevant to *{proposal_topic}*. A well-cited proposal needs at least {MIN_PAPERS}.\n\n"
                        f"A few options:\n"
                        f"- **Chat with me first** about your topic (e.g. *'What is known about FOXO in fat body aging?'*) "
                        f"to pull more papers, then click Plan Research\n"
                        f"- **Share specific papers** you already know — paste a PMID or title and I'll incorporate them\n"
                        f"- **Broaden your topic** slightly if it's very niche\n\n"
                        f"What would you like to do?"
                    ),
                    'proposal': None,
                    'ready_for_export': False
                }

            print(f"  ✅ Generating proposal with {len(proposal_papers)} papers...")
            # Update the planner context with the enriched paper set + clarifications
            self.planner.set_context_from_chat(proposal_topic, proposal_genes, proposal_papers, clarifications=user_message)

            proposal = self.planner.generate_proposal(
                topic=proposal_topic,
                genes=proposal_genes,
                papers=proposal_papers,
                user_clarifications=user_message,
                conversation_history=self.conversation_history
            )
            formatted = self.planner.format_proposal_for_chat(proposal)
            return {
                'type': 'proposal',
                'content': formatted,
                'proposal': proposal,
                'ready_for_export': True
            }

        # Always ask clarifying questions on first proposal generation.
        # Only skip if we already asked and are now receiving the answer (awaiting_clarification).
        # Having genes + papers is not sufficient — we need researcher-specific context
        # to avoid generic proposals.
        print("  ❓ Asking clarifying questions for better specificity...")
        question = self.planner.generate_clarifying_questions(topic, genes)
        self.planner.set_context_from_chat(topic, genes, papers)
        self.awaiting_clarification = True
        return {
            'type': 'clarification',
            'content': f"🔬 **Research Planner activated!**\n\nBefore I build your proposal, {question}",
            'proposal': None,
            'ready_for_export': False
        }

    def chat(self, user_message: str, force_planning: bool = False) -> Dict:
        """
        Main chat function. Returns dict with:
        - response: text to display
        - is_planning: bool
        - proposal: dict or None
        - ready_for_export: bool
        """
        print(f"\n{'='*70}\nQuery: {user_message[:80]}\n{'='*70}")

        # Check for planning intent
        has_prior_context = bool(self.last_topic or self.last_papers)
        is_planning, confidence = detect_planning_intent(user_message, has_prior_context)

        # If awaiting clarification, only route to planner if message looks like
        # an answer (short, descriptive) rather than a new research question.
        is_clarification_answer = (
            self.awaiting_clarification
            and self.planning_mode
            and not is_planning  # not a new planning request
            and len(user_message.split()) < 60  # reasonably short
            and '?' not in user_message[:20]  # doesn't start as a question
        )

        if force_planning or is_planning or is_clarification_answer:
            self.planning_mode = True
            result = self.handle_planning_request(user_message, force=force_planning)

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
        self.awaiting_clarification = False  # reset if user went back to normal chat
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
            pubmed_papers = self.search_pubmed(user_message, max_results=num_needed - len(all_papers))
            if pubmed_papers:
                all_papers.extend(pubmed_papers)

        # Filter chat papers for relevance too
        if all_papers:
            all_papers = self.filter_relevant_papers(all_papers, user_message, genes)

        # Cache for planner
        if all_papers:
            self.last_papers = all_papers
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

CRITICAL: When you see "RELEVANT PUBLICATIONS" section with papers listed:
1. YOU MUST cite these papers throughout your answer
2. Reference them like: [1], [2], [3], etc.
3. Create a "References" section at the end listing all papers with links
4. Format: "[1] Author et al. (Year). Title. [Link](URL)"

IMPORTANT INSTRUCTIONS:
- Always cite the papers provided - DO NOT ignore them
- If papers are provided, you MUST reference them multiple times
- Include clickable links to papers in the References section
- Be confident citing these papers - they were found for this specific query
- Include a References section at the end of EVERY response with papers and links
- If no papers are provided, just answer from your knowledge

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
            model="claude-sonnet-4-20250514",
            max_tokens=3000,
            system=system_prompt,
            messages=self.conversation_history
        )

        assistant_message = response.content[0].text
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