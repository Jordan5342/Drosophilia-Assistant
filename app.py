from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import os
from datetime import datetime
import json
import signal
from drosophila_assistant_enhanced import DrosophilaAssistant

app = Flask(__name__, static_folder='static')
CORS(app)

# Initialize assistant
api_key = os.environ.get("ANTHROPIC_API_KEY")
assistant = DrosophilaAssistant(api_key) if api_key else None

# Store conversation history (in production, use a database)
conversations = {}

# Timeout handling
class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException()

@app.route('/')
def home():
    return send_from_directory('static', 'index.html')

@app.route('/api/chat', methods=['POST'])
def chat():
    if not assistant:
        return jsonify({'error': 'API key not configured. Please set ANTHROPIC_API_KEY environment variable.'}), 500
    
    try:
        data = request.get_json()
        user_message = data.get('message', '')
        session_id = data.get('session_id', 'default')
        
        if not user_message:
            return jsonify({'error': 'No message provided'}), 400
        
        # Get or create conversation for this session
        if session_id not in conversations:
            conversations[session_id] = DrosophilaAssistant(api_key)
        
        session_assistant = conversations[session_id]
        
        # Set a timeout alarm (110 seconds - less than gunicorn's 120)
        # This only works on Unix-like systems
        try:
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(110)
        except (AttributeError, ValueError):
            # signal.SIGALRM not available on Windows
            pass
        
        try:
            response = session_assistant.chat(user_message)
            
            # Cancel the alarm
            try:
                signal.alarm(0)
            except (AttributeError, ValueError):
                pass
            
            return jsonify({
                'response': response,
                'timestamp': datetime.now().isoformat()
            })
        
        except TimeoutException:
            # Cancel the alarm
            try:
                signal.alarm(0)
            except (AttributeError, ValueError):
                pass
            
            return jsonify({
                'error': 'Request took too long to process. This can happen with complex queries. Try:\n• A more specific question\n• Fewer genes to compare\n• A simpler topic',
                'timeout': True
            }), 408
    
    except Exception as e:
        print(f"Error in chat endpoint: {str(e)}")
        import traceback
        traceback.print_exc()
        
        # Cancel alarm on any error
        try:
            signal.alarm(0)
        except (AttributeError, ValueError):
            pass
        
        return jsonify({'error': str(e)}), 500

@app.route('/api/reset', methods=['POST'])
def reset():
    try:
        data = request.get_json()
        session_id = data.get('session_id', 'default')
        
        if session_id in conversations:
            conversations[session_id].reset_conversation()
        
        return jsonify({'status': 'success'})
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/export', methods=['POST'])
def export_conversation():
    try:
        data = request.get_json()
        session_id = data.get('session_id', 'default')
        
        if session_id not in conversations:
            return jsonify({'error': 'No conversation found'}), 404
        
        session_assistant = conversations[session_id]
        history = session_assistant.conversation_history
        
        # Format as markdown
        markdown = "# Drosophila Research Conversation\n\n"
        markdown += f"Exported: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
        markdown += "---\n\n"
        
        for msg in history:
            role = "**You**" if msg['role'] == 'user' else "**Assistant**"
            # Strip the context additions from user messages
            content = msg['content']
            # If content has FlyBase or PubMed sections, only show the original query
            if '======' in content:
                content = content.split('\n\n')[0]
            markdown += f"{role}:\n\n{content}\n\n---\n\n"
        
        return jsonify({
            'markdown': markdown,
            'filename': f"drosophila_conversation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/health', methods=['GET'])
def health():
    return jsonify({
        'status': 'healthy',
        'api_configured': assistant is not None,
        'features': {
            'pubmed_search': True,
            'flybase_lookup': True,
            'multi_gene_support': True,
            'result_filtering': True,
            'pathway_emphasis': True
        }
    })

@app.route('/api/stats', methods=['GET'])
def stats():
    """Get statistics about active conversations"""
    return jsonify({
        'active_sessions': len(conversations),
        'total_conversations': sum(len(conv.conversation_history) for conv in conversations.values()) // 2
    })

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    
    if not api_key:
        print("\n" + "="*60)
        print("⚠️  WARNING: ANTHROPIC_API_KEY not set!")
        print("="*60)
        print("The app will start but won't work until you set your API key.")
        print("Set it as an environment variable:")
        print("  export ANTHROPIC_API_KEY=your_key_here")
        print("="*60 + "\n")
    
    print("\n" + "="*60)
    print("🪰 Drosophila Research Assistant - ENHANCED")
    print("="*60)
    print(f"\n✓ Server starting on http://localhost:{port}")
    print("\n✨ Enhanced Features:")
    print("  📚 PubMed literature search with relevance filtering")
    print("  🧬 FlyBase gene lookup with variant name support")
    print("  🔬 Genetic pathway emphasis in responses")
    print("  🧬 Multi-gene comparison support")
    print("  🤖 Claude AI integration (Sonnet 4)")
    print("  💾 Memory-optimized (50% reduction)")
    print("  ⏱️  Timeout protection (110s)")
    print("\n💡 Tips:")
    print("  • Open the URL in your browser")
    print("  • Ask about genes, papers, or research topics")
    print("  • Try comparisons: 'Compare Notch and Delta'")
    print("  • Use gene names flexibly: 'dFOXO', 'foxo', 'FOXO' all work")
    print("  • Press Ctrl+C to stop\n")
    print("⚙️  Configuration:")
    print(f"  • API Key: {'✓ Configured' if api_key else '✗ Missing'}")
    print(f"  • Port: {port}")
    print(f"  • Debug: False")
    print("="*60 + "\n")
    
    app.run(host='0.0.0.0', port=port, debug=False)