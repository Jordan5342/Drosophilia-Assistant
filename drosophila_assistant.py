import anthropic
import os
from Bio import Entrez
import json
import requests
from typing import Optional, Dict, List
import re

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
    
    # Circadian
    'per': ('per', 'period', 'FBgn0003068'),
    'period': ('per', 'period', 'FBgn0003068'),
    'tim': ('tim', 'timeless', 'FBgn0014396'),
    'timeless': ('tim', 'timeless', 'FBgn0014396'),
    'clk': ('Clk', 'Clock', 'FBgn0023076'),
    'clock': ('Clk', 'Clock', 'FBgn0023076'),
    
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
        
        # Method 1: Check known genes (instant, no network call)
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
        
        # Method 2: Quick web lookup (timeout 8s)
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
    
    def determine_paper_count(self, user_query: str) -> int:
        """Determine number of papers needed"""
        query_lower = user_query.lower()
        broad = ['top', 'list', 'overview', 'review', 'genes in', 'pathways', 'mechanisms', 'role of', 'involved in', 'best', 'major']
        specific = ['what is', 'tell me about', 'how does', 'function of']
        
        if any(ind in query_lower for ind in broad):
            return 8
        if any(ind in query_lower for ind in specific):
            return 5
        return 6
    
    def get_flybase_publications(self, fbgn: str, max_results: int = 10) -> List[Dict]:
        """Fetch publications from FlyBase with proxy handling and fallback"""
        try:
            print(f"  📚 Fetching publications for {fbgn}...")
            url = f"https://flybase.org/reports/{fbgn}"
            
            # Try with explicit session and headers to avoid proxy issues
            try:
                session = requests.Session()
                session.headers.update({
                    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
                })
                response = session.get(url, timeout=5, allow_redirects=True)
                
                # Handle 202 Accepted status (page still processing)
                if response.status_code == 202:
                    print(f"    ⚠️  FlyBase still processing (202), retrying...")
                    import time
                    time.sleep(1)
                    response = session.get(url, timeout=5, allow_redirects=True)
                
                if response.status_code != 200:
                    print(f"    ⚠️  FlyBase returned status {response.status_code}, using fallback")
                    return self._fallback_pubmed_search(fbgn, max_results)
                
            except (requests.exceptions.Timeout, requests.exceptions.ConnectionError, requests.exceptions.ProxyError):
                print(f"    ⚠️  FlyBase connection failed, falling back to direct PubMed search")
                return self._fallback_pubmed_search(fbgn, max_results)
            
            # Extract PMIDs from page - use the pattern that actually works
            pmid_pattern = r'/pubmed/(\d+)'
            pmids = re.findall(pmid_pattern, response.text, re.IGNORECASE)
            
            # Remove duplicates while preserving order
            pmids = list(dict.fromkeys(pmids))[:max_results]
            
            print(f"    Found {len(pmids)} PMIDs on FlyBase page")
            
            if not pmids:
                print(f"    ℹ️  No PMIDs found on page, falling back to PubMed search")
                return self._fallback_pubmed_search(fbgn, max_results)
            
            # Try to fetch from PubMed
            try:
                print(f"    Fetching details from PubMed for {len(pmids)} articles...")
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
                        print(f"    ⚠️  Error parsing article {i+1}: {e}")
                        continue
                
                if publications:
                    print(f"  ✅ Retrieved {len(publications)} publications from FlyBase")
                return publications
                
            except Exception as e:
                print(f"    ⚠️  PubMed fetch error: {e}")
                # Fallback: return basic PMID links
                print(f"    Returning {len(pmids)} basic links as fallback")
                return [{'pmid': pmid, 'title': f'Publication', 'authors': 'See PubMed', 'year': 'N/A', 'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", 'source': 'FlyBase'} for pmid in pmids[:max_results]]
            
        except Exception as e:
            print(f"    ⚠️  Error fetching publications: {e}")
            return self._fallback_pubmed_search(fbgn, max_results)
    
    def _fallback_pubmed_search(self, fbgn: str, max_results: int) -> List[Dict]:
        """Fallback: search PubMed directly using gene name from FBgn"""
        try:
            # Extract gene name from known genes database or use FBgn
            gene_name = fbgn  # Default
            for key, (symbol, name, fbn) in KNOWN_GENES.items():
                if fbn == fbgn:
                    gene_name = symbol
                    break
            
            print(f"    Searching PubMed for '{gene_name}' (Drosophila)...")
            full_query = f"Drosophila AND {gene_name}"
            
            handle = Entrez.esearch(db="pubmed", term=full_query, retmax=max_results + 2, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            pmids = record.get("IdList", [])[:max_results]
            
            if not pmids:
                print(f"    No results found")
                return []
            
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
                        'source': 'FlyBase (via PubMed)'
                    })
                except Exception:
                    continue
            
            if publications:
                print(f"  ✅ Retrieved {len(publications)} publications via PubMed")
            return publications
            
        except Exception as e:
            print(f"    ⚠️  Fallback search failed: {e}")
            return []
    
    def search_pubmed(self, query, max_results=5):
        """Simplified PubMed search"""
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
        """Extract gene names from query - improved to avoid false positives"""
        potential_genes = []
        exclude = {'aging', 'gene', 'genes', 'protein', 'proteins', 'what', 'is', 'are', 'tell', 'me', 'about', 'a', 'n', 'the', 'of', 'in'}
        
        text_lower = text.lower()
        
        # First: check for explicit gene mentions in known database (only multi-char genes)
        for gene_key in KNOWN_GENES.keys():
            if len(gene_key) > 1 and gene_key in text_lower and gene_key not in exclude:
                potential_genes.append(gene_key)
        
        # Second: look for patterns like "gene X" or "X gene"
        pattern = r'(?:gene\s+([a-zA-Z]{2,})|([a-zA-Z]{2,})\s+gene)'
        matches = re.findall(pattern, text, re.IGNORECASE)
        for match in matches:
            gene = match[0] or match[1]
            if gene.lower() not in exclude and len(gene) > 1:
                potential_genes.append(gene)
        
        # Remove duplicates while preserving order
        seen = set()
        unique = []
        for g in potential_genes:
            g_lower = g.lower()
            if g_lower not in seen and g_lower not in exclude:
                unique.append(g_lower)
                seen.add(g_lower)
        
        # Limit to 3 genes per query
        return unique[:3]
    
    def format_papers(self, papers):
        """Format papers for Claude with explicit links and citation format"""
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
        formatted += f"INSTRUCTIONS: Reference these papers in your answer using [1], [2], etc.\n"
        formatted += f"Include the URL links and create a References section at the end\n"
        formatted += "="*70 + "\n\n"
        return formatted
    
    def chat(self, user_message):
        """Main chat function"""
        print(f"\n{'='*70}\nQuery: {user_message[:80]}...\n{'='*70}")
        
        all_papers = []
        
        # Extract and lookup genes
        genes = self.extract_gene_names(user_message)
        if genes:
            print(f"  Found genes: {genes}")
            for gene in genes:
                gene_info = search_flybase_fixed(gene)
                if gene_info:
                    num_papers = self.determine_paper_count(user_message)
                    pubs = self.get_flybase_publications(gene_info['fbgn'], max_results=num_papers)
                    if pubs:
                        all_papers.extend(pubs)
        
        # Supplement with PubMed if needed
        num_needed = self.determine_paper_count(user_message)
        if len(all_papers) < num_needed:
            print(f"  Supplementing with PubMed...")
            pubmed_papers = self.search_pubmed(user_message, max_results=num_needed - len(all_papers))
            if pubmed_papers:
                all_papers.extend(pubmed_papers)
        
        # Format papers (if any)
        publication_context = self.format_papers(all_papers) if all_papers else ""
        
        # Build message
        enhanced_message = user_message
        if publication_context:
            enhanced_message += "\n\n" + publication_context
        
        self.conversation_history.append({"role": "user", "content": enhanced_message})
        if len(self.conversation_history) > 10:
            self.conversation_history = self.conversation_history[-10:]
        
        # System prompt that emphasizes citing papers
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

Remember: The papers in "RELEVANT PUBLICATIONS" are provided for YOU to cite and use!"""
        
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
        return assistant_message
    
    def reset_conversation(self):
        """Clear history"""
        self.conversation_history = []
        self.mentioned_genes = set()