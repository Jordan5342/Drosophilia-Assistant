import anthropic
import os
from Bio import Entrez
import json
import requests
from typing import Optional, Dict, List
import re

# Configure Entrez (PubMed API)
Entrez.email = "jordanszt@icloud.com"  # Required by NCBI

class DrosophilaAssistant:
    def __init__(self, api_key):
        self.client = anthropic.Anthropic(api_key=api_key)
        self.conversation_history = []
        
    def reformulate_query_for_pubmed(self, user_query: str) -> str:
        """Use Claude to reformulate user queries into optimal PubMed search terms"""
        try:
            # Use Claude to intelligently reformulate the query
            response = self.client.messages.create(
                model="claude-sonnet-4-5-20250929",
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
            
            return reformulated
            
        except Exception as e:
            print(f"  ⚠️  Error in query reformulation: {e}")
            # Fallback to original query if Claude call fails
            return user_query
        
    def search_pubmed(self, query, max_results=5):
        """Search PubMed for Drosophila-related papers"""
        try:
            # Reformulate query for better PubMed results
            reformulated = self.reformulate_query_for_pubmed(query)
            
            # Add Drosophila to the search query
            full_query = f"Drosophila AND ({reformulated})"
            
            print(f"  🔍 Original: {query[:60]}...")
            print(f"  🔍 PubMed query: {full_query}")
            
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=full_query, retmax=max_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            
            if not id_list:
                print("  ℹ️  No PubMed results found")
                return []
            
            print(f"  📥 Fetching {len(id_list)} papers...")
            
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
                        'abstract': abstract[:500] + '...' if len(abstract) > 500 else abstract,
                        'pmid': str(pmid),
                        'year': year,
                        'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    })
                except Exception as e:
                    print(f"  ⚠️  Error parsing article: {e}")
                    continue
            
            print(f"  ✅ Successfully retrieved {len(papers)} papers")
            return papers
            
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
            formatted += f"📄 Paper {i}:\n"
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
    
    def search_flybase(self, gene_name: str) -> Optional[Dict]:
        """Search FlyBase for gene information"""
        try:
            # Clean up gene name
            gene_name = gene_name.strip().lower()
            
            print(f"  🔍 FlyBase lookup: {gene_name}")
            
            # Try FlyBase lookup service
            lookup_url = f"https://flybase.org/api/v1.0/genes/autocomplete?query={gene_name}"
            
            response = requests.get(lookup_url, timeout=10)
            
            # Check if response is valid
            if response.status_code != 200:
                print(f"  ℹ️  FlyBase returned status {response.status_code} for '{gene_name}'")
                return None
            
            # Check if response has content
            if not response.text or response.text.strip() == '':
                print(f"  ℹ️  '{gene_name}' not found in FlyBase (empty response)")
                return None
            
            # Try to parse JSON
            try:
                data = response.json()
            except ValueError:
                print(f"  ℹ️  '{gene_name}' not found in FlyBase (invalid response)")
                return None
            
            if data and len(data) > 0:
                # Get first match
                first_match = data[0]
                
                fbgn = first_match.get('id', '')
                symbol = first_match.get('symbol', gene_name)
                name = first_match.get('name', 'Unknown')
                
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
            else:
                print(f"  ℹ️  '{gene_name}' not found in FlyBase")
            
            return None
                
        except requests.exceptions.Timeout:
            print(f"  ⚠️  FlyBase request timed out")
            return None
        except requests.exceptions.RequestException as e:
            print(f"  ⚠️  FlyBase connection error")
            return None
        except Exception as e:
            print(f"  ℹ️  Could not look up '{gene_name}' in FlyBase")
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
        """Extract potential gene names from user query"""
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
            'cancer', 'tumor', 'mutation', 'mutant', 'allele'
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
        
        # Pattern 4: Common Drosophila genes mentioned directly
        common_genes = ['notch', 'delta', 'wingless', 'hedgehog', 'decapentaplegic', 
                       'engrailed', 'even-skipped', 'fushi-tarazu', 'white', 'yellow',
                       'sevenless', 'bride', 'boss', 'torpedo', 'gurken', 'oskar',
                       'nanos', 'pumilio', 'bicoid', 'hunchback', 'kruppel', 'giant',
                       'knirps', 'tailless', 'p53', 'ras', 'myc', 'src', 'abl',
                       'inr', 'dfoxo', 'foxo', 'sir2', 'rpd3', 'mth', 'methuselah',
                       'sod', 'catalase', 'tor', 'pten', 'pi3k', 'akt']
        
        text_lower = text.lower()
        for gene in common_genes:
            if gene in text_lower:
                potential_genes.append(gene)
        
        # Return unique genes (limit to first 2 to avoid too many API calls)
        unique_genes = []
        seen = set()
        for gene in potential_genes:
            gene_lower = gene.lower()
            if gene_lower not in seen and gene_lower not in exclude_words:
                unique_genes.append(gene)
                seen.add(gene_lower)
        
        return unique_genes[:2]
    
    def should_search_pubmed(self, message: str) -> bool:
        """Determine if we should search PubMed - Always return True"""
        # Always search PubMed for every query to find relevant research
        return True
    
    def should_search_flybase(self, message: str) -> bool:
        """Determine if we should search FlyBase - Search if message contains potential gene names"""
        # Always try to extract and search for genes
        # FlyBase will return nothing if no genes found, which is fine
        return True
    
    def chat(self, user_message):
        """Chat with Claude, searching PubMed and FlyBase as needed"""
        
        print(f"\n{'='*70}")
        print(f"User query: {user_message[:100]}...")
        print(f"{'='*70}")
        
        publication_context = ""
        flybase_context = ""
        
        # ALWAYS search PubMed for relevant papers
        print("📚 Searching PubMed for relevant research...")
        papers = self.search_pubmed(user_message, max_results=8)
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
        
        # ALWAYS try to extract and search for genes on FlyBase
        print("🧬 Checking for gene mentions...")
        potential_genes = self.extract_gene_names(user_message)
        
        if potential_genes:
            print(f"  Potential genes detected: {potential_genes}")
            for gene in potential_genes:
                gene_info = self.search_flybase(gene)
                if gene_info:
                    flybase_context += self.format_flybase_info(gene_info)
                    break  # Only use first successful match
        else:
            print("  ℹ️  No specific gene names detected")
        
        # Build system prompt
        system_prompt = """You are a specialized AI assistant for Drosophila melanogaster (fruit fly) research.

Your expertise includes:
- Drosophila genetics, development, and molecular biology
- Gene nomenclature and function
- Experimental techniques in fly research
- Classic and modern Drosophila studies
- Developmental biology and signaling pathways

IMPORTANT: For EVERY user query, the system automatically searches:
1. PubMed for recent publications
2. FlyBase for gene information

CRITICAL FORMATTING INSTRUCTIONS:

**NEVER use HTML tags in your response.** Always use Markdown formatting:
- For links: Use [Link Text](URL) format
- For emphasis: Use **bold** for important terms

CITATION INSTRUCTIONS:

**IF you see "RELEVANT PUBLICATIONS FROM PUBMED" with papers listed:**
   ✅ This means papers WERE found - you MUST cite them!
   - Reference papers throughout your answer
   - Use markdown links: [Author et al., Year](URL)  
   - Create a "References" section at the end
   - Format: "1. Author et al. (Year). Title. PMID: 12345. [View on PubMed](URL)"
   - Synthesize information from the papers
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

**General Guidelines:**
- Be confident when citing papers that ARE provided
- Don't claim lack of access if papers ARE in the context
- If papers are provided, they are recent and relevant - use them!
- Synthesize across multiple papers when available
- Always be accurate and acknowledge true uncertainty

Remember: If papers are in your context window between the PUBLICATIONS markers, you DO have access to them - cite them confidently!"""

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
            model="claude-sonnet-4-5-20250929",
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
        print("🔄 Conversation history cleared")