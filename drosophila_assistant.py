import anthropic
import os
from Bio import Entrez
import json
import requests
from typing import Optional, Dict, List

# Configure Entrez (PubMed API)
Entrez.email = "jordanszt@icloud.com"  # Required by NCBI

class DrosophilaAssistant:
    def __init__(self, api_key):
        self.client = anthropic.Anthropic(api_key=api_key)
        self.conversation_history = []
        
    def search_pubmed(self, query, max_results=10):
        """Search PubMed for Drosophila-related papers"""
        try:
            # Add Drosophila to the search query
            full_query = f"Drosophila AND ({query})"
            
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=full_query, retmax=max_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            
            if not id_list:
                return "No publications found."
            
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
                    year = article_data.get('ArticleDate', [{}])[0].get('Year', 'N/A') if article_data.get('ArticleDate') else \
                           medline.get('DateCompleted', {}).get('Year', 'N/A')
                    
                    papers.append({
                        'title': title,
                        'authors': ', '.join(authors) + ' et al.' if authors else 'Unknown authors',
                        'abstract': abstract,
                        'pmid': str(pmid),
                        'year': year,
                        'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    })
                except Exception as e:
                    print(f"Error parsing article: {e}")
                    continue
            
            return papers
            
        except Exception as e:
            return f"Error searching PubMed: {str(e)}"
    
    def format_papers_for_claude(self, papers):
        """Format papers into a readable string for Claude"""
        if isinstance(papers, str):
            return papers
        
        formatted = "Here are relevant publications from PubMed:\n\n"
        for i, paper in enumerate(papers, 1):
            formatted += f"Paper {i}:\n"
            formatted += f"Title: {paper['title']}\n"
            formatted += f"Authors: {paper['authors']}\n"
            formatted += f"Year: {paper['year']}\n"
            formatted += f"PMID: {paper['pmid']}\n"
            formatted += f"URL: {paper['url']}\n"
            formatted += f"Abstract: {paper['abstract']}\n\n"
        
        return formatted
    
    def chat(self, user_message, search_publications=True):
        """Chat with Claude, optionally searching publications first"""
        
        # Check if the query seems to need publication search
        search_keywords = ['paper', 'study', 'research', 'publication', 'gene', 'protein', 'pathway', 'mutation', 
                          'find', 'recent', 'what', 'how', 'development', 'signaling', 'regulation']
        should_search = search_publications and any(keyword in user_message.lower() for keyword in search_keywords)
        
        publication_context = ""
        if should_search:
            print("Searching PubMed for relevant publications...")
            papers = self.search_pubmed(user_message)
            publication_context = self.format_papers_for_claude(papers)
        
        # Check if query mentions specific genes - search FlyBase
        flybase_context = ""
        gene_keywords = ['gene', 'protein', 'mutant', 'allele', 'what is', 'tell me about', 'function']
        if any(keyword in user_message.lower() for keyword in gene_keywords):
            print("Searching FlyBase for gene information...")
            potential_genes = self.extract_gene_names(user_message)
            
            for gene in potential_genes:
                gene_info = self.search_flybase(gene)
                if gene_info:
                    flybase_context += self.format_flybase_info(gene_info)
                    break  # Only get info for first valid gene found
        
        # Build the system prompt
        system_prompt = """You are a specialized AI assistant for Drosophila melanogaster (fruit fly) research. 

Your expertise includes:
- Drosophila genetics, development, and molecular biology
- Gene nomenclature and function
- Experimental techniques in fly research
- Classic and modern Drosophila studies
- Comparative biology and model organism research

IMPORTANT: When publications are provided, you MUST:
1. Cite each paper with its PMID number in your answer
2. Include the full PubMed URL for each paper mentioned: https://pubmed.ncbi.nlm.nih.gov/PMID/
3. List the paper titles and authors when referencing them
4. Create a "References" section at the end with all cited papers

When FlyBase gene information is provided:
1. Use the official gene symbol and FlyBase ID (FBgn)
2. Include the FlyBase link for reference
3. Mention synonyms if relevant to help users find the gene

Always provide accurate scientific information and acknowledge uncertainty when appropriate.
Use clear scientific language but explain complex concepts when needed."""

        # Add publication context and FlyBase context to the user message if available
        enhanced_message = user_message
        if flybase_context:
            enhanced_message = f"{user_message}\n\n{flybase_context}"
        if publication_context:
            enhanced_message = f"{enhanced_message}\n\n{publication_context}"
        
        # Add to conversation history
        self.conversation_history.append({
            "role": "user",
            "content": enhanced_message
        })
        
        # Call Claude API
        response = self.client.messages.create(
            model="claude-sonnet-4-5-20250929",
            max_tokens=4096,
            system=system_prompt,
            messages=self.conversation_history
        )
        
        assistant_message = response.content[0].text
        
        # Add response to history
        self.conversation_history.append({
            "role": "assistant",
            "content": assistant_message
        })
        
        return assistant_message
    
    def reset_conversation(self):
        """Clear conversation history"""
        self.conversation_history = []
    
    def search_flybase(self, gene_name: str) -> Optional[Dict]:
        """Search FlyBase for gene information"""
        try:
            # FlyBase API endpoint
            base_url = "https://api.flybase.org/api/v1.0"
            
            # Search for gene by symbol or name
            search_url = f"{base_url}/gene/symbol/{gene_name}"
            
            response = requests.get(search_url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                
                # Extract key information
                gene_info = {
                    'symbol': data.get('symbol', gene_name),
                    'name': data.get('name', 'Unknown'),
                    'fbgn': data.get('FBgn', 'Unknown'),
                    'summary': data.get('automated_gene_summary', 'No summary available'),
                    'synonyms': data.get('synonyms', []),
                    'url': f"https://flybase.org/reports/{data.get('FBgn', '')}"
                }
                
                return gene_info
            else:
                # Try autocomplete search as fallback
                autocomplete_url = f"{base_url}/autocomplete?query={gene_name}"
                response = requests.get(autocomplete_url, timeout=10)
                
                if response.status_code == 200:
                    results = response.json()
                    if results and len(results) > 0:
                        # Get first match
                        first_match = results[0]
                        fbgn = first_match.get('id', '')
                        
                        if fbgn:
                            # Fetch full gene details
                            gene_url = f"{base_url}/gene/{fbgn}"
                            gene_response = requests.get(gene_url, timeout=10)
                            
                            if gene_response.status_code == 200:
                                data = gene_response.json()
                                gene_info = {
                                    'symbol': data.get('symbol', gene_name),
                                    'name': data.get('name', 'Unknown'),
                                    'fbgn': fbgn,
                                    'summary': data.get('automated_gene_summary', 'No summary available'),
                                    'synonyms': data.get('synonyms', []),
                                    'url': f"https://flybase.org/reports/{fbgn}"
                                }
                                return gene_info
                
                return None
                
        except Exception as e:
            print(f"Error searching FlyBase: {str(e)}")
            return None
    
    def format_flybase_info(self, gene_info: Dict) -> str:
        """Format FlyBase gene information for Claude"""
        if not gene_info:
            return ""
        
        formatted = f"\n--- FlyBase Information ---\n"
        formatted += f"Gene Symbol: {gene_info['symbol']}\n"
        formatted += f"Full Name: {gene_info['name']}\n"
        formatted += f"FlyBase ID: {gene_info['fbgn']}\n"
        
        if gene_info.get('synonyms'):
            synonyms = ', '.join(gene_info['synonyms'][:5])  # Show first 5
            formatted += f"Synonyms: {synonyms}\n"
        
        formatted += f"Summary: {gene_info['summary']}\n"
        formatted += f"FlyBase Link: {gene_info['url']}\n"
        formatted += "---\n\n"
        
        return formatted
    
    def extract_gene_names(self, text: str) -> List[str]:
        """Extract potential gene names from user query"""
        # Common Drosophila gene patterns
        import re
        
        # Look for common gene name patterns
        # Italicized or capitalized gene names
        potential_genes = []
        
        # Pattern 1: Words in context that might be genes
        words = text.split()
        for word in words:
            # Remove punctuation
            clean_word = word.strip('.,!?;:()[]{}"\'-')
            
            # Drosophila genes often:
            # - Start with lowercase (wingless, notch)
            # - Have mixed case (Notch, Delta)
            # - Are short (2-5 letters common)
            if len(clean_word) >= 2:
                potential_genes.append(clean_word)
        
        # Pattern 2: Look for "gene X" patterns
        gene_pattern = r'gene\s+(\w+)|(\w+)\s+gene'
        matches = re.findall(gene_pattern, text, re.IGNORECASE)
        for match in matches:
            gene = match[0] or match[1]
            if gene:
                potential_genes.append(gene)
        
        # Return unique genes (limit to first 3 to avoid too many API calls)
        unique_genes = list(dict.fromkeys(potential_genes))
        return unique_genes[:3]


# Example usage
if __name__ == "__main__":
    # Get API key from environment variable
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    
    if not api_key:
        print("Please set ANTHROPIC_API_KEY environment variable")
        exit(1)
    
    assistant = DrosophilaAssistant(api_key)
    
    print("Drosophila Research Assistant")
    print("=" * 50)
    print("Ask questions about Drosophila research!")
    print("Type 'quit' to exit, 'reset' to clear conversation\n")
    
    while True:
        user_input = input("You: ").strip()
        
        if user_input.lower() == 'quit':
            break
        elif user_input.lower() == 'reset':
            assistant.reset_conversation()
            print("Conversation reset!\n")
            continue
        elif not user_input:
            continue
        
        try:
            response = assistant.chat(user_input)
            print(f"\nAssistant: {response}\n")
        except Exception as e:
            print(f"Error: {str(e)}\n")