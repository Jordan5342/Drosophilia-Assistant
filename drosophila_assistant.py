import anthropic
import os
from Bio import Entrez
import json
import requests
from typing import Optional, Dict, List, Tuple
import re
from collections import Counter

# Configure Entrez (PubMed API)
Entrez.email = "jordanszt@icloud.com"  # Required by NCBI

class DrosophilaAssistant:
    def __init__(self, api_key):
        self.client = anthropic.Anthropic(api_key=api_key)
        self.conversation_history = []
        self.mentioned_genes = set()  # Track genes discussed in session
        self.mentioned_pathways = set()  # Track pathways discussed
    
    def determine_paper_count(self, user_query: str) -> int:
        """Determine optimal number of papers based on query type"""
        query_lower = user_query.lower()
        
        # Broad/review queries need more papers
        broad_indicators = ['top', 'list', 'overview', 'review', 'genes in', 'pathways', 
                           'mechanisms', 'role of', 'involved in', 'best', 'major']
        
        # Specific queries need fewer papers
        specific_indicators = ['what is', 'tell me about', 'how does', 'function of']
        
        # Check for broad query
        if any(indicator in query_lower for indicator in broad_indicators):
            return 15  # More papers for filtering
        
        # Check for specific query
        if any(indicator in query_lower for indicator in specific_indicators):
            return 8  # Fewer papers for specific questions
        
        # Default: moderate
        return 10
    
    def calculate_relevance_score(self, paper: Dict, query_terms: List[str]) -> float:
        """
        Calculate relevance score for a paper based on query terms.
        Returns score between 0 and 1.
        """
        title = paper.get('title', '').lower()
        abstract = paper.get('abstract', '').lower()
        
        # Combine title and abstract, weight title more heavily
        title_weight = 3.0
        abstract_weight = 1.0
        
        score = 0.0
        max_score = 0.0
        
        for term in query_terms:
            term_lower = term.lower()
            term_score = 0.0
            
            # Count occurrences in title (weighted more)
            title_count = title.count(term_lower)
            term_score += title_count * title_weight
            
            # Count occurrences in abstract
            abstract_count = abstract.count(term_lower)
            term_score += abstract_count * abstract_weight
            
            score += term_score
            # Maximum possible score for this term (appearing 3 times in title, 5 in abstract)
            max_score += (3 * title_weight + 5 * abstract_weight)
        
        # Normalize to 0-1 range
        if max_score > 0:
            return min(score / max_score, 1.0)
        return 0.0
    
    def filter_and_rank_papers(self, papers: List[Dict], query: str, top_n: int = 8) -> List[Dict]:
        """
        Filter and rank papers by relevance, returning top N results.
        """
        if not papers:
            return []
        
        # Extract key terms from query
        query_terms = self.extract_key_terms(query)
        
        print(f"  🎯 Filtering with key terms: {', '.join(query_terms)}")
        
        # Calculate relevance score for each paper
        scored_papers = []
        for paper in papers:
            score = self.calculate_relevance_score(paper, query_terms)
            paper['relevance_score'] = score
            scored_papers.append(paper)
        
        # Sort by relevance score (highest first)
        scored_papers.sort(key=lambda x: x['relevance_score'], reverse=True)
        
        # Filter out very low relevance papers (score < 0.05)
        filtered_papers = [p for p in scored_papers if p['relevance_score'] >= 0.05]
        
        # Return top N papers
        top_papers = filtered_papers[:top_n]
        
        if top_papers:
            print(f"  ✅ Filtered to {len(top_papers)} most relevant papers")
            print(f"  📊 Relevance scores: {[f'{p[\"relevance_score\"]:.2f}' for p in top_papers[:3]]}")
        
        return top_papers
    
    def extract_key_terms(self, query: str) -> List[str]:
        """Extract key scientific terms from query for relevance scoring"""
        # Remove common stop words
        stop_words = {'the', 'a', 'an', 'and', 'or', 'but', 'in', 'on', 'at', 'to', 'for',
                     'of', 'with', 'by', 'from', 'up', 'about', 'into', 'through', 'during',
                     'what', 'which', 'who', 'when', 'where', 'why', 'how', 'is', 'are',
                     'was', 'were', 'be', 'been', 'being', 'do', 'does', 'did', 'tell',
                     'me', 'give', 'show', 'find', 'top', 'best', 'most'}
        
        # Tokenize and filter
        words = re.findall(r'\b\w+\b', query.lower())
        key_terms = [w for w in words if w not in stop_words and len(w) > 2]
        
        return key_terms
    
    def reformulate_query_for_pubmed(self, user_query: str) -> str:
        """Use Claude to reformulate user queries into optimal PubMed search terms"""
        try:
            print(f"  🤖 Using AI to reformulate query...")
            
            # Use Claude to intelligently reformulate the query
            response = self.client.messages.create(
                model="claude-sonnet-4-20250514",
                max_tokens=150,
                temperature=0,
                system="""You are an expert at converting natural language questions into optimal PubMed search queries for Drosophila research.

Your task: Convert the user's question into 3-8 relevant keywords/phrases optimized for PubMed searching.

Guidelines:
- Extract core scientific concepts (genes, pathways, processes, phenotypes)
- Use proper biological terminology
- Remove conversational language ("tell me", "what are", "top 10", etc.)
- Include relevant synonyms when helpful
- Keep it concise - output ONLY the keywords separated by spaces
- Focus on terms that would appear in scientific papers
- Include pathway names and biological processes

Examples:
User: "What are the top 10 genes involved in aging?"
Output: aging longevity lifespan extension genes mechanisms senescence

User: "Tell me about Notch signaling in wing development"
Output: Notch signaling wing development imaginal disc patterning

User: "How does p53 work in flies?"
Output: p53 Drosophila apoptosis DNA damage tumor suppressor

User: "What genes control sleep?"
Output: sleep circadian rhythm rest behavior neurogenetics

Output ONLY the search keywords, nothing else.""",
                messages=[
                    {
                        "role": "user",
                        "content": user_query
                    }
                ]
            )
            
            reformulated = response.content[0].text.strip()
            
            # Fallback: if Claude's response is too long or seems wrong, use original query
            if len(reformulated) > 200 or '\n' in reformulated:
                print(f"  ⚠️  Claude reformulation too long, using original query")
                return user_query
            
            print(f"  ✅ Reformulated to: {reformulated}")
            return reformulated
            
        except Exception as e:
            print(f"  ⚠️  Error in query reformulation: {e}")
            # Fallback to original query if Claude call fails
            return user_query
        
    def search_pubmed(self, query, max_results=5):
        """Search PubMed for Drosophila-related papers with filtering"""
        try:
            # Reformulate query for better PubMed results
            reformulated = self.reformulate_query_for_pubmed(query)
            
            # Add Drosophila to the search query
            full_query = f"Drosophila AND ({reformulated})"
            
            print(f"  🔍 Original: {query[:60]}...")
            print(f"  🔍 PubMed query: {full_query}")
            
            # Search PubMed - get more results initially for filtering
            initial_results = max_results + 7  # Get extra papers to filter
            handle = Entrez.esearch(db="pubmed", term=full_query, retmax=initial_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            
            if not id_list:
                print("  ℹ️  No PubMed results found")
                return []
            
            print(f"  📥 Fetching {len(id_list)} papers for filtering...")
            
            # Fetch details for each paper
            handle = Entrez.efetch(db="pubmed", id=id_list, rettype="abstract", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            papers = []
            for article in records['PubmedArticle']:
                try:
                    medline = article['MedlineCitation']
                    article_data = medline['Article']
                    
                    title = article_data.get('ArticleTitle', 'No title')
                    abstract = article_data.get('Abstract', {}).get('AbstractText', ['No abstract available'])
                    if isinstance(abstract, list):
                        abstract = ' '.join(str(text) for text in abstract)
                    
                    authors = []
                    if 'AuthorList' in article_data:
                        for author in article_data['AuthorList'][:3]:  # First 3 authors
                            if 'LastName' in author:
                                authors.append(f"{author.get('LastName', '')} {author.get('Initials', '')}")
                    
                    pmid = medline['PMID']
                    
                    # Get year from multiple possible locations
                    year = 'N/A'
                    if 'ArticleDate' in article_data and article_data['ArticleDate']:
                        year = article_data['ArticleDate'][0].get('Year', 'N/A')
                    elif 'Journal' in article_data and 'JournalIssue' in article_data['Journal']:
                        pub_date = article_data['Journal']['JournalIssue'].get('PubDate', {})
                        year = pub_date.get('Year', 'N/A')
                    
                    papers.append({
                        'title': title,
                        'authors': ', '.join(authors) + ' et al.' if authors else 'Unknown authors',
                        'abstract': abstract[:800] + '...' if len(abstract) > 800 else abstract,
                        'pmid': str(pmid),
                        'year': year,
                        'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    })
                except Exception as e:
                    print(f"  ⚠️  Error parsing article: {e}")
                    continue
            
            # Filter and rank papers by relevance
            filtered_papers = self.filter_and_rank_papers(papers, query, top_n=max_results)
            
            print(f"  ✅ Successfully retrieved {len(filtered_papers)} filtered papers")
            return filtered_papers
            
        except Exception as e:
            print(f"  ❌ PubMed search error: {str(e)}")
            return []
    
    def format_papers_for_claude(self, papers):
        """Format papers into a readable string for Claude"""
        if not papers:
            return ""
        
        formatted = "\n" + "="*70 + "\n"
        formatted += "RELEVANT PUBLICATIONS FROM PUBMED\n"
        formatted += "="*70 + "\n\n"
        
        for i, paper in enumerate(papers, 1):
            formatted += f"📄 Paper {i}"
            # Add relevance score if available
            if 'relevance_score' in paper:
                formatted += f" (Relevance: {paper['relevance_score']:.2f})"
            formatted += ":\n"
            formatted += f"Title: {paper['title']}\n"
            formatted += f"Authors: {paper['authors']}\n"
            formatted += f"Year: {paper['year']}\n"
            formatted += f"PMID: {paper['pmid']}\n"
            formatted += f"URL: {paper['url']}\n"
            formatted += f"Abstract: {paper['abstract']}\n"
            formatted += "\n" + "-"*70 + "\n\n"
        
        formatted += "="*70 + "\n"
        formatted += "END OF PUBLICATIONS\n"
        formatted += "="*70 + "\n\n"
        return formatted
    
    def search_flybase_with_variants(self, gene_name: str) -> Optional[Dict]:
        """
        Search FlyBase trying multiple name variants.
        Handles common gene name formats and synonyms.
        """
        # Generate variants
        variants = self.generate_gene_name_variants(gene_name)
        
        print(f"  🔍 FlyBase lookup: {gene_name}")
        if len(variants) > 1:
            print(f"  🔄 Will try variants: {', '.join(list(variants)[:5])}")
        
        # Try each variant
        for variant in variants:
            result = self.search_flybase(variant)
            if result:
                return result
        
        print(f"  ℹ️  '{gene_name}' not found in FlyBase (tried {len(variants)} variants)")
        return None
    
    def generate_gene_name_variants(self, gene_name: str) -> List[str]:
        """
        Generate common variants of a gene name.
        Examples: FOXO -> foxo, Foxo, dFOXO, dfoxo
                  p53 -> p53, dmp53, Dmp53
        """
        gene_name = gene_name.strip()
        variants = [gene_name]  # Start with original
        
        # Lowercase version
        variants.append(gene_name.lower())
        
        # Capitalized version
        variants.append(gene_name.capitalize())
        
        # Uppercase version
        variants.append(gene_name.upper())
        
        # Add 'd' prefix (Drosophila prefix)
        if not gene_name.lower().startswith('d'):
            variants.append(f"d{gene_name}")
            variants.append(f"d{gene_name.lower()}")
            variants.append(f"D{gene_name.lower()}")
        
        # Remove 'd' prefix if present
        if gene_name.lower().startswith('d') and len(gene_name) > 1:
            variants.append(gene_name[1:])
            variants.append(gene_name[1:].lower())
        
        # Handle dashes and underscores
        if '-' in gene_name:
            variants.append(gene_name.replace('-', ''))
        if '_' in gene_name:
            variants.append(gene_name.replace('_', ''))
        
        # Return unique variants, prioritizing original and lowercase
        seen = set()
        unique_variants = []
        for v in variants:
            if v.lower() not in seen:
                unique_variants.append(v)
                seen.add(v.lower())
        
        return unique_variants
    
    def search_flybase(self, gene_name: str) -> Optional[Dict]:
        """Search FlyBase for gene information"""
        try:
            # Clean up gene name
            gene_name = gene_name.strip()
            
            # Try FlyBase lookup service
            lookup_url = f"https://flybase.org/api/v1.0/genes/autocomplete?query={gene_name}"
            
            response = requests.get(lookup_url, timeout=10)
            
            # Check if response is valid
            if response.status_code != 200:
                return None
            
            # Check if response has content
            if not response.text or response.text.strip() == '':
                return None
            
            # Try to parse JSON
            try:
                data = response.json()
            except ValueError:
                return None
            
            if data and len(data) > 0:
                # Get first match
                first_match = data[0]
                
                fbgn = first_match.get('id', '')
                symbol = first_match.get('symbol', gene_name)
                name = first_match.get('name', 'Unknown')
                
                # Track this gene in session
                self.mentioned_genes.add(symbol.lower())
                
                gene_info = {
                    'symbol': symbol,
                    'name': name,
                    'fbgn': fbgn,
                    'summary': f"Gene symbol: {symbol}. This is a Drosophila gene.",
                    'synonyms': first_match.get('synonyms', []),
                    'url': f"https://flybase.org/reports/{fbgn}"
                }
                
                print(f"  ✅ Found: {symbol} ({fbgn})")
                return gene_info
            
            return None
                
        except requests.exceptions.Timeout:
            return None
        except requests.exceptions.RequestException as e:
            return None
        except Exception as e:
            return None
    
    def format_flybase_info(self, gene_info: Dict) -> str:
        """Format FlyBase gene information for Claude"""
        if not gene_info:
            return ""
        
        formatted = "\n" + "="*70 + "\n"
        formatted += "FLYBASE GENE INFORMATION\n"
        formatted += "="*70 + "\n\n"
        
        formatted += f"🧬 Gene Symbol: {gene_info['symbol']}\n"
        formatted += f"Full Name: {gene_info['name']}\n"
        formatted += f"FlyBase ID: {gene_info['fbgn']}\n"
        
        if gene_info.get('synonyms'):
            synonyms = ', '.join(gene_info['synonyms'][:5])
            formatted += f"Synonyms: {synonyms}\n"
        
        formatted += f"\nFlyBase Link: {gene_info['url']}\n"
        
        formatted += "\n" + "="*70 + "\n"
        formatted += "END OF FLYBASE INFO\n"
        formatted += "="*70 + "\n\n"
        
        return formatted
    
    def extract_gene_names(self, text: str) -> List[str]:
        """
        Extract potential gene names from user query.
        Enhanced to detect multiple genes and comparison queries.
        """
        potential_genes = []
        
        # Words to exclude (common words that aren't genes)
        exclude_words = {
            'aging', 'gene', 'genes', 'protein', 'proteins', 'the', 'a', 'an', 
            'this', 'that', 'what', 'how', 'why', 'when', 'where', 'who',
            'is', 'are', 'was', 'were', 'be', 'been', 'being',
            'do', 'does', 'did', 'have', 'has', 'had',
            'can', 'could', 'will', 'would', 'should', 'may', 'might',
            'tell', 'find', 'show', 'give', 'make', 'take',
            'top', 'best', 'most', 'many', 'some', 'all',
            'development', 'signaling', 'pathway', 'function', 'role',
            'cancer', 'tumor', 'mutation', 'mutant', 'allele',
            'compare', 'versus', 'vs', 'between', 'and', 'or'
        }
        
        # Pattern 1: "gene X" or "X gene"
        gene_pattern = r'(?:gene\s+(\w+)|(\w+)\s+gene)'
        matches = re.findall(gene_pattern, text, re.IGNORECASE)
        for match in matches:
            gene = match[0] or match[1]
            if gene and gene.lower() not in exclude_words:
                potential_genes.append(gene)
        
        # Pattern 2: "what is X" where X might be a gene
        what_pattern = r'what\s+(?:is|does)\s+(\w+)'
        matches = re.findall(what_pattern, text, re.IGNORECASE)
        for match in matches:
            if match.lower() not in exclude_words:
                potential_genes.append(match)
        
        # Pattern 3: "tell me about X"
        tell_pattern = r'(?:tell me about|about)\s+(\w+)'
        matches = re.findall(tell_pattern, text, re.IGNORECASE)
        for match in matches:
            if match.lower() not in exclude_words:
                potential_genes.append(match)
        
        # Pattern 4: "compare X and Y" or "X vs Y"
        compare_pattern = r'(?:compare|versus|vs\.?)\s+(\w+)\s+(?:and|versus|vs\.?)\s+(\w+)'
        matches = re.findall(compare_pattern, text, re.IGNORECASE)
        for match in matches:
            if match[0].lower() not in exclude_words:
                potential_genes.append(match[0])
            if match[1].lower() not in exclude_words:
                potential_genes.append(match[1])
        
        # Pattern 5: "X and Y" (only if they look like gene names)
        and_pattern = r'\b([A-Z][a-z]+\d*)\s+and\s+([A-Z][a-z]+\d*)\b'
        matches = re.findall(and_pattern, text)
        for match in matches:
            if match[0].lower() not in exclude_words:
                potential_genes.append(match[0])
            if match[1].lower() not in exclude_words:
                potential_genes.append(match[1])
        
        # Pattern 6: Common Drosophila genes mentioned directly
        common_genes = ['notch', 'delta', 'wingless', 'hedgehog', 'decapentaplegic', 
                       'engrailed', 'even-skipped', 'fushi-tarazu', 'white', 'yellow',
                       'sevenless', 'bride', 'boss', 'torpedo', 'gurken', 'oskar',
                       'nanos', 'pumilio', 'bicoid', 'hunchback', 'kruppel', 'giant',
                       'knirps', 'tailless', 'p53', 'ras', 'myc', 'src', 'abl',
                       'inr', 'dfoxo', 'foxo', 'sir2', 'rpd3', 'mth', 'methuselah',
                       'sod', 'catalase', 'tor', 'pten', 'pi3k', 'akt', 'dmp53',
                       'jak', 'stat', 'egfr', 'insulin', 'dilp']
        
        text_lower = text.lower()
        for gene in common_genes:
            if gene in text_lower:
                potential_genes.append(gene)
        
        # Return unique genes (limit to first 5 for multi-gene queries)
        unique_genes = []
        seen = set()
        for gene in potential_genes:
            gene_lower = gene.lower()
            if gene_lower not in seen and gene_lower not in exclude_words:
                unique_genes.append(gene)
                seen.add(gene_lower)
        
        return unique_genes[:5]  # Increased from 2 to 5 for multi-gene support
    
    def chat(self, user_message):
        """Chat with Claude, searching PubMed and FlyBase as needed"""
        
        print(f"\n{'='*70}")
        print(f"User query: {user_message[:100]}...")
        print(f"{'='*70}")
        
        publication_context = ""
        flybase_context = ""
        
        # ALWAYS search PubMed for relevant papers
        print("📚 Searching PubMed for relevant research...")
        num_papers = self.determine_paper_count(user_message)
        print(f"  📊 Retrieving {num_papers} papers based on query type...")
        papers = self.search_pubmed(user_message, max_results=num_papers)
        if papers:
            publication_context = self.format_papers_for_claude(papers)
        else:
            print("  ℹ️  No recent papers found for this query")
            # Tell Claude explicitly that no papers were found
            publication_context = "\n" + "="*70 + "\n"
            publication_context += "PUBMED SEARCH RESULT\n"
            publication_context += "="*70 + "\n\n"
            publication_context += "No recent publications were found in PubMed for this specific query.\n"
            publication_context += "You may answer from your knowledge base and note that no recent papers were found.\n"
            publication_context += "\n" + "="*70 + "\n"
        
        # ALWAYS try to extract and search for genes on FlyBase (multi-gene support)
        print("🧬 Checking for gene mentions...")
        potential_genes = self.extract_gene_names(user_message)
        
        if potential_genes:
            print(f"  Potential genes detected: {potential_genes}")
            genes_found = 0
            for gene in potential_genes:
                gene_info = self.search_flybase_with_variants(gene)
                if gene_info:
                    flybase_context += self.format_flybase_info(gene_info)
                    genes_found += 1
                    # For comparison queries, get multiple genes
                    if genes_found >= 3 and 'compare' not in user_message.lower():
                        break
            
            if genes_found > 0:
                print(f"  ✅ Found {genes_found} gene(s) in FlyBase")
        else:
            print("  ℹ️  No specific gene names detected")
        
        # Build enhanced system prompt with genetic pathway emphasis
        system_prompt = """You are a specialized AI assistant for Drosophila melanogaster (fruit fly) research.

Your expertise includes:
- Drosophila genetics, development, and molecular biology
- Gene nomenclature and function
- Experimental techniques in fly research
- Classic and modern Drosophila studies
- Developmental biology and signaling pathways

IMPORTANT: For EVERY user query, the system automatically searches:
1. PubMed for recent publications (filtered for relevance)
2. FlyBase for gene information (including multiple genes for comparisons)

CRITICAL GENETIC PATHWAY EMPHASIS:

When discussing genes and their functions, you MUST emphasize:

**1. PATHWAY CONTEXT:**
   - Identify which signaling pathway(s) the gene belongs to
   - Describe upstream regulators and downstream targets
   - Explain how the gene fits into the broader regulatory network
   - Mention key pathway interactions (activators, inhibitors, feedback loops)

**2. ORTHOLOGUES AND CONSERVATION:**
   - Always mention mammalian/human orthologues when known
   - Discuss conservation of function across species
   - Highlight what fly research has taught us about human biology
   - Note any differences between fly and mammalian systems

**3. DEVELOPMENTAL AND TISSUE CONTEXT:**
   - Specify when and where the gene is expressed
   - Describe developmental stages or tissues affected
   - Explain temporal regulation if relevant
   - Connect to specific biological processes

**4. FUNCTIONAL RELATIONSHIPS:**
   - Describe genetic interactions (synthetic lethality, suppression, enhancement)
   - Identify parallel pathways or redundant functions
   - Note epistatic relationships
   - Mention compensatory mechanisms

**5. CROSS-SPECIES KNOWLEDGE TRANSFER:**
   - Explicitly connect fly findings to broader biology
   - Discuss how fly models inform disease research
   - Highlight conserved mechanisms discovered in Drosophila
   - Note clinical relevance where applicable

CRITICAL FORMATTING INSTRUCTIONS:

**NEVER use HTML tags in your response.** Always use Markdown formatting:
- For links: Use [Link Text](URL) format
- For emphasis: Use **bold** for important terms
- For pathway names: Use **bold** to highlight them

CITATION INSTRUCTIONS:

**IF you see "RELEVANT PUBLICATIONS FROM PUBMED" with papers listed:**
   ✅ This means papers WERE found - you MUST cite them!
   - Reference papers throughout your answer
   - Use markdown links: [Author et al., Year](URL)  
   - Create a "References" section at the end
   - Format: "1. Author et al. (Year). Title. PMID: 12345. [View on PubMed](URL)"
   - Synthesize information from the papers
   - Note: Papers are pre-filtered for relevance - higher scores indicate better matches
   - DO NOT say you don't have access to papers!

**IF you see "No publications found" or no PUBLICATIONS section:**
   ❌ Papers were not found
   - Answer from your knowledge base
   - Briefly note: "No recent papers were found on this specific query in PubMed"
   - Continue with your expert knowledge
   - Suggest alternative search terms if relevant

**IF you see "FLYBASE GENE INFORMATION":**
   - Use the official gene symbol and FlyBase ID
   - Include link: [View on FlyBase](URL)
   - Cite the official information provided
   - For multiple genes: compare and contrast their functions and pathways

**FOR COMPARISON QUERIES (multiple genes):**
   - Structure your answer to compare pathway roles
   - Highlight functional similarities and differences
   - Discuss potential interactions between the genes
   - Note if they're in parallel or antagonistic pathways

**General Guidelines:**
- Lead with pathway and mechanistic context
- Always mention human relevance and orthologues
- Be confident when citing papers that ARE provided
- Don't claim lack of access if papers ARE in the context
- If papers are provided, they are recent and relevant - use them!
- Synthesize across multiple papers when available
- Always be accurate and acknowledge true uncertainty
- Structure responses to emphasize biological pathways and mechanisms

Remember: Your goal is not just to describe what a gene does, but to place it in the context of cellular pathways, evolutionary conservation, and broader biological significance!"""

        # Build enhanced message
        enhanced_message = user_message
        
        if publication_context or flybase_context:
            enhanced_message = f"{user_message}\n\n"
            if flybase_context:
                enhanced_message += flybase_context
            if publication_context:
                enhanced_message += publication_context
        
        # Add to conversation history
        self.conversation_history.append({
            "role": "user",
            "content": enhanced_message
        })
        
        print("🤖 Calling Claude API...")
        
        # Call Claude API
        response = self.client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=4096,
            system=system_prompt,
            messages=self.conversation_history
        )
        
        assistant_message = response.content[0].text
        
        # Add response to history (without the extra context)
        self.conversation_history.append({
            "role": "assistant",
            "content": assistant_message
        })
        
        print("✅ Response generated")
        print(f"{'='*70}\n")
        
        return assistant_message
    
    def reset_conversation(self):
        """Clear conversation history"""
        self.conversation_history = []
        self.mentioned_genes = set()
        self.mentioned_pathways = set()
        print("🔄 Conversation history cleared")