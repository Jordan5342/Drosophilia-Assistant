import anthropic
import os
from Bio import Entrez
import json

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
            formatted += f"Abstract: {paper['abstract'][:500]}...\n\n"
        
        return formatted
    
    def chat(self, user_message, search_publications=True):
        """Chat with Claude, optionally searching publications first"""
        
        # Check if the query seems to need publication search
        search_keywords = ['paper', 'study', 'research', 'publication', 'gene', 'protein', 'pathway', 'mutation']
        should_search = search_publications and any(keyword in user_message.lower() for keyword in search_keywords)
        
        publication_context = ""
        if should_search:
            print("Searching PubMed for relevant publications...")
            papers = self.search_pubmed(user_message)
            publication_context = self.format_papers_for_claude(papers)
        
        # Build the system prompt
        system_prompt = """You are a specialized AI assistant for Drosophila melanogaster (fruit fly) research. 

Your expertise includes:
- Drosophila genetics, development, and molecular biology
- Gene nomenclature and function
- Experimental techniques in fly research
- Classic and modern Drosophila studies
- Comparative biology and model organism research

When publications are provided, cite them appropriately using PMID numbers. 
Always provide accurate scientific information and acknowledge uncertainty when appropriate.
Use clear scientific language but explain complex concepts when needed."""

        # Add publication context to the user message if available
        enhanced_message = user_message
        if publication_context:
            enhanced_message = f"{user_message}\n\n{publication_context}"
        
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