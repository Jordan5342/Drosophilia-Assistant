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
        
    def search_pubmed(self, query, max_results=5):
        """Search PubMed for Drosophila-related papers"""
        try:
            # Add Drosophila to the search query
            full_query = f"Drosophila AND ({query})"
            
            print(f"  🔍 PubMed query: {full_query}")
            
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=full_query, retmax=max_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            
            if not id_list:
                print("  ℹ️  No PubMed results found")
                return []
            
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
            
            print(f"  ✅ Found {len(papers)} papers")
            return papers
            
        except Exception as e:
            print(f"  ❌ Error searching PubMed: {str(e)}")
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
            
            if response.status_code == 200:
                data = response.json()
                
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
                    print(f"  ℹ️  No FlyBase results for '{gene_name}'")
            
            return None
                
        except Exception as e:
            print(f"  ❌ Error searching FlyBase: {str(e)}")
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
        
        # Pattern 1: "gene X" or "X gene"
        gene_pattern = r'(?:gene\s+(\w+)|(\w+)\s+gene)'
        matches = re.findall(gene_pattern, text, re.IGNORECASE)
        for match in matches:
            gene = match[0] or match[1]
            if gene and gene.lower() not in ['the', 'a', 'this', 'that', 'my', 'your']:
                potential_genes.append(gene)
        
        # Pattern 2: "what is X" where X might be a gene
        what_pattern = r'what\s+(?:is|does)\s+(\w+)'
        matches = re.findall(what_pattern, text, re.IGNORECASE)
        for match in matches:
            if match.lower() not in ['the', 'a', 'this', 'that', 'it']:
                potential_genes.append(match)
        
        # Pattern 3: "tell me about X"
        tell_pattern = r'(?:tell me about|about)\s+(\w+)'
        matches = re.findall(tell_pattern, text, re.IGNORECASE)
        for match in matches:
            if match.lower() not in ['the', 'a', 'this', 'that']:
                potential_genes.append(match)
        
        # Pattern 4: Common Drosophila genes mentioned directly
        common_genes = ['notch', 'delta', 'wingless', 'hedgehog', 'decapentaplegic', 
                       'engrailed', 'even-skipped', 'fushi-tarazu', 'white', 'yellow',
                       'sevenless', 'bride', 'boss', 'torpedo', 'gurken', 'oskar',
                       'nanos', 'pumilio', 'bicoid', 'hunchback', 'kruppel', 'giant',
                       'knirps', 'tailless', 'gap', 'pair-rule', 'segment', 'polarity']
        
        text_lower = text.lower()
        for gene in common_genes:
            if gene in text_lower:
                potential_genes.append(gene)
        
        # Return unique genes (limit to first 2 to avoid too many API calls)
        unique_genes = list(dict.fromkeys(potential_genes))
        return unique_genes[:2]
    
    def should_search_pubmed(self, message: str) -> bool:
        """Determine if we should search PubMed"""
        search_keywords = [
            'paper', 'papers', 'study', 'studies', 'research', 'publication', 
            'publications', 'recent', 'find papers', 'literature', 'review',
            'what research', 'studies on', 'published', 'article', 'articles',
            'evidence', 'data', 'findings', 'experiments', 'experimental'
        ]
        
        message_lower = message.lower()
        return any(keyword in message_lower for keyword in search_keywords)
    
    def should_search_flybase(self, message: str) -> bool:
        """Determine if we should search FlyBase"""
        gene_keywords = [
            'gene', 'genes', 'what is', 'tell me about', 'function of', 
            'protein', 'mutant', 'mutation', 'allele', 'locus', 'chromosome'
        ]
        
        message_lower = message.lower()
        return any(keyword in message_lower for keyword in gene_keywords)
    
    def chat(self, user_message):
        """Chat with Claude, searching PubMed and FlyBase as needed"""
        
        print(f"\n{'='*70}")
        print(f"User query: {user_message[:100]}...")
        print(f"{'='*70}")
        
        publication_context = ""
        flybase_context = ""
        
        # Search PubMed if needed
        if self.should_search_pubmed(user_message):
            print("📚 Searching PubMed...")
            papers = self.search_pubmed(user_message, max_results=5)
            if papers:
                publication_context = self.format_papers_for_claude(papers)
        
        # Search FlyBase if needed
        if self.should_search_flybase(user_message):
            print("🧬 Searching FlyBase...")
            potential_genes = self.extract_gene_names(user_message)
            
            if potential_genes:
                print(f"  Potential genes: {potential_genes}")
                for gene in potential_genes:
                    gene_info = self.search_flybase(gene)
                    if gene_info:
                        flybase_context += self.format_flybase_info(gene_info)
                        break  # Only use first successful match
        
        # Build system prompt
        system_prompt = """You are a specialized AI assistant for Drosophila melanogaster (fruit fly) research.

Your expertise includes:
- Drosophila genetics, development, and molecular biology
- Gene nomenclature and function
- Experimental techniques in fly research
- Classic and modern Drosophila studies
- Developmental biology and signaling pathways

CRITICAL INSTRUCTIONS FOR CITING SOURCES:

1. When PubMed publications are provided (marked by "RELEVANT PUBLICATIONS"):
   - Reference papers naturally in your answer
   - ALWAYS include clickable links using markdown: [Author et al., Year](URL)
   - Create a "References" section at the end
   - Format each reference as: "Author et al. (Year). Title. PMID: 12345. [View on PubMed](URL)"

2. When FlyBase gene information is provided (marked by "FLYBASE GENE INFORMATION"):
   - Use the official gene symbol and FlyBase ID
   - Include the FlyBase link as: [View on FlyBase](URL)
   - Mention synonyms if they help clarify

3. If NO external data is provided:
   - Answer from your knowledge base
   - Be clear that information is from your training data
   - Suggest they could search for recent papers if relevant

Always be accurate, acknowledge uncertainty, and provide clear scientific explanations."""

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