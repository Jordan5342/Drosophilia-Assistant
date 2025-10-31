from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import os
from datetime import datetime
import json
from drosophila_assistant import DrosophilaAssistant

app = Flask(__name__, static_folder='static')
CORS(app)

# Initialize assistant
api_key = os.environ.get("ANTHROPIC_API_KEY")
assistant = DrosophilaAssistant(api_key) if api_key else None

# Store conversation history (in production, use a database)
conversations = {}

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
        response = session_assistant.chat(user_message)
        
        return jsonify({
            'response': response,
            'timestamp': datetime.now().isoformat()
        })
    
    except Exception as e:
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
            markdown += f"{role}:\n\n{msg['content']}\n\n---\n\n"
        
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
        'api_configured': assistant is not None
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
    print("🪰 Drosophila Research Assistant")
    print("="*60)
    print(f"\n✓ Server starting on http://localhost:{port}")
    print("\n💡 Tips:")
    print("  • Open the URL in your browser")
    print("  • Share with colleagues - they just need the link!")
    print("  • Press Ctrl+C to stop\n")
    
    app.run(host='0.0.0.0', port=port, debug=False)