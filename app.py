from flask import Flask, request, jsonify, send_from_directory, Response, stream_with_context
from flask_cors import CORS
import os
from datetime import datetime
import json
import signal
from drosophila_assistant import DrosophilaAssistant
from agents.pipeline import AgentPipeline
import database

DEBUG_AGENTS = os.environ.get("DEBUG_AGENTS", "false").lower() == "true"

app = Flask(__name__, static_folder='static')
CORS(app)

# Runs under both gunicorn and direct python execution
if os.environ.get("DATABASE_URL"):
    try:
        database.init_db()
    except Exception as _db_init_exc:
        print(f"⚠️  Database init failed at startup: {_db_init_exc}")

api_key = os.environ.get("ANTHROPIC_API_KEY")
assistant = DrosophilaAssistant(api_key) if api_key else None

# In-memory cache per worker — DB is the source of truth.
# On a cache miss, get_session() hydrates from DB before returning.
conversations = {}


class TimeoutException(Exception):
    pass


def timeout_handler(signum, frame):
    raise TimeoutException()


def get_session(session_id: str) -> DrosophilaAssistant:
    if session_id in conversations:
        return conversations[session_id]

    sess = DrosophilaAssistant(api_key)

    try:
        state = database.load_session(session_id)
        if state:
            sess.last_topic = state["last_topic"]
            sess.last_genes = state["last_genes"]
            sess.planning_mode = state["planning_mode"]
            sess.awaiting_clarification = state["awaiting_clarification"]

        history = database.load_conversation(session_id)
        if history:
            sess.conversation_history = history

        proposal_data = database.load_proposal(session_id)
        if proposal_data:
            sess.planner.current_proposal = proposal_data["proposal"]
            ctx = proposal_data["proposal_context"]
            if ctx:
                sess.planner.proposal_context = ctx

        papers = database.load_papers(session_id)
        if papers:
            sess.last_papers = papers

    except Exception as exc:
        print(f"  ⚠️  DB load error for {session_id} (starting fresh): {exc}")

    conversations[session_id] = sess
    return sess


def _persist_session(session_id: str, sess: DrosophilaAssistant) -> None:
    """Save all mutable session state back to DB. Non-fatal on failure."""
    try:
        database.save_session(
            session_id,
            last_topic=sess.last_topic or "",
            last_genes=list(sess.last_genes) if sess.last_genes else [],
            planning_mode=bool(sess.planning_mode),
            awaiting_clarification=bool(sess.awaiting_clarification),
        )
        database.save_conversation(session_id, sess.conversation_history or [])
        if sess.planner.current_proposal is not None:
            database.save_proposal(
                session_id,
                sess.planner.current_proposal,
                sess.planner.proposal_context or {},
            )
        if sess.last_papers:
            database.save_papers(session_id, sess.last_papers)
    except Exception as exc:
        print(f"  ⚠️  DB save error for {session_id} (non-fatal): {exc}")


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
        force_planning = data.get('force_planning', False)  # from "Plan Research" button
        client_awaiting = data.get('awaiting_clarification', False)  # client-side state

        if not user_message:
            return jsonify({'error': 'No message provided'}), 400

        session_assistant = get_session(session_id)

        # Restore client-side planning state if server session was reset (e.g. after restart)
        print(f"  🔍 State: client_awaiting={client_awaiting}, server_awaiting={session_assistant.awaiting_clarification}, force={force_planning}")
        if client_awaiting and not session_assistant.awaiting_clarification:
            print(f"  ♻️  Restoring awaiting_clarification from client state")
            session_assistant.awaiting_clarification = True
            session_assistant.planning_mode = True

        try:
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(110)
        except (AttributeError, ValueError):
            pass

        try:
            result = session_assistant.chat(user_message, force_planning=force_planning)

            try:
                signal.alarm(0)
            except (AttributeError, ValueError):
                pass

            _persist_session(session_id, session_assistant)

            return jsonify({
                'response': result['response'],
                'is_planning': result['is_planning'],
                'ready_for_export': result['ready_for_export'],
                'has_proposal': result['proposal'] is not None,
                'proposal_json': result['proposal'] if result['proposal'] else None,
                'timestamp': datetime.now().isoformat()
            })

        except TimeoutException:
            try:
                signal.alarm(0)
            except (AttributeError, ValueError):
                pass
            return jsonify({
                'error': 'Request took too long. Try a more specific question or fewer genes.',
                'timeout': True
            }), 408

    except Exception as e:
        print(f"Error in chat endpoint: {str(e)}")
        import traceback
        traceback.print_exc()
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
            del conversations[session_id]
        try:
            database.delete_session(session_id)
        except Exception as exc:
            print(f"  ⚠️  DB delete error on reset (non-fatal): {exc}")
        return jsonify({'status': 'success'})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/chat_history', methods=['POST'])
def chat_history():
    try:
        data = request.get_json()
        session_id = data.get('session_id', 'default')
        sess = get_session(session_id)
        return jsonify({'history': sess.conversation_history or []})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/export', methods=['POST'])
def export_conversation():
    try:
        data = request.get_json()
        session_id = data.get('session_id', 'default')

        session_assistant = get_session(session_id)
        history = session_assistant.conversation_history

        if not history:
            return jsonify({'error': 'No conversation found'}), 404

        markdown = "# Drosophila Research Conversation\n\n"
        markdown += f"Exported: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n---\n\n"

        for msg in history:
            role = "**You**" if msg['role'] == 'user' else "**Assistant**"
            content = msg['content']
            if '======' in content:
                content = content.split('\n\n')[0]
            markdown += f"{role}:\n\n{content}\n\n---\n\n"

        return jsonify({
            'markdown': markdown,
            'filename': f"drosophila_conversation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/export_proposal', methods=['POST'])
def export_proposal():
    """Export the current proposal as a .docx file using docx-js via Node."""
    try:
        data = request.get_json()
        session_id = data.get('session_id', 'default')

        session_assistant = get_session(session_id)
        proposal = session_assistant.planner.get_proposal_for_export()

        if not proposal:
            return jsonify({'error': 'No proposal to export. Generate a proposal first.'}), 404

        import base64
        from docx_export import build_proposal_docx

        docx_bytes = build_proposal_docx(proposal)
        encoded = base64.b64encode(docx_bytes).decode('utf-8')
        filename = f"research_proposal_{datetime.now().strftime('%Y%m%d_%H%M%S')}.docx"

        return jsonify({
            'docx_base64': encoded,
            'filename': filename,
            'status': 'success'
        })

    except Exception as e:
        print(f"Export proposal error: {str(e)}")
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500



@app.route('/api/design_experiment', methods=['POST'])
def design_experiment():
    """Generate a detailed experiment design for a specific aim."""
    try:
        data = request.get_json()
        session_id = data.get('session_id', 'default')
        aim_number = data.get('aim_number', 1)
        client_proposal = data.get('proposal', None)  # fallback from browser cache

        session_assistant = get_session(session_id)
        proposal = session_assistant.planner.get_proposal_for_export()

        # If server lost the proposal (restart), restore from client-side cache
        if not proposal and client_proposal:
            print(f"  ♻️  Restoring proposal from client cache")
            session_assistant.planner.current_proposal = client_proposal
            proposal = client_proposal

        if not proposal:
            return jsonify({'error': 'No proposal found. Generate a proposal first.'}), 404

        aims = proposal.get('specific_aims', [])
        aim = next((a for a in aims if a.get('aim_number') == aim_number), None)

        if not aim:
            return jsonify({'error': f'Aim {aim_number} not found in proposal.'}), 404

        print(f"  🧪 Generating experiment design for Aim {aim_number}...")
        design = session_assistant.planner.generate_experiment_design(aim, proposal)
        formatted = session_assistant.planner.format_experiment_design(design)

        return jsonify({
            'design': design,
            'formatted': formatted,
            'aim_number': aim_number,
            'aim_title': aim.get('title', ''),
            'status': 'success'
        })

    except Exception as e:
        print(f"Design experiment error: {str(e)}")
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/export_experiment_design', methods=['POST'])
def export_experiment_design():
    """Export a single experiment design aim as .docx"""
    try:
        data = request.get_json()
        session_id = data.get('session_id', 'default')
        design = data.get('design', None)  # full design JSON from client

        if not design:
            return jsonify({'error': 'No experiment design provided.'}), 400

        import base64
        from docx_export import build_experiment_design_docx

        docx_bytes = build_experiment_design_docx(design)
        encoded = base64.b64encode(docx_bytes).decode('utf-8')
        aim_num = design.get('aim_number', 1)
        filename = f"aim{aim_num}_experiment_design_{datetime.now().strftime('%Y%m%d_%H%M%S')}.docx"

        return jsonify({'docx_base64': encoded, 'filename': filename, 'status': 'success'})

    except Exception as e:
        print(f"Export design error: {str(e)}")
        import traceback; traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/api/agent_pipeline', methods=['POST'])
def agent_pipeline():
    """
    Run the three-agent pipeline (Literature → Hypothesis → Critic) with SSE streaming.
    Set DEBUG_AGENTS=true in env to receive full per-iteration JSON in the stream.
    """
    if not assistant:
        return jsonify({'error': 'API key not configured. Set ANTHROPIC_API_KEY.'}), 500

    try:
        data = request.get_json()
        query = data.get('query', '').strip()
        session_id = data.get('session_id', 'default')

        if not query:
            return jsonify({'error': 'No query provided'}), 400

        session_assistant = get_session(session_id)

        # Extract topic and gene list using existing helpers
        topic = session_assistant.planner.extract_topic_from_message(query) or query[:120]
        genes = session_assistant.extract_gene_names(query)
        if not genes and session_assistant.last_genes:
            genes = session_assistant.last_genes

        print(f"\n{'='*70}")
        print(f"[Agent Pipeline] Query: {query[:80]}")
        print(f"[Agent Pipeline] Topic: {topic[:60]} | Genes: {genes}")
        print(f"[Agent Pipeline] DEBUG_AGENTS={DEBUG_AGENTS}")
        print(f"{'='*70}")

        pipeline = AgentPipeline(
            client=session_assistant.client,
            fetch_literature_fn=session_assistant.fetch_literature_for_proposal
        )

        def generate():
            try:
                for event in pipeline.run_streaming(
                    query=query,
                    topic=topic,
                    genes=genes,
                    debug=DEBUG_AGENTS
                ):
                    yield f"event: {event['type']}\ndata: {json.dumps(event['data'])}\n\n"
            except Exception as exc:
                import traceback
                traceback.print_exc()
                yield f"event: error\ndata: {json.dumps({'message': str(exc)})}\n\n"

        return Response(
            stream_with_context(generate()),
            content_type='text/event-stream',
            headers={
                'Cache-Control': 'no-cache',
                'X-Accel-Buffering': 'no',
                'Connection': 'keep-alive'
            }
        )

    except Exception as exc:
        print(f"Agent pipeline error: {exc}")
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(exc)}), 500


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
            'pathway_emphasis': True,
            'research_planner': True,
            'proposal_export_docx': True
        }
    })


@app.route('/api/stats', methods=['GET'])
def stats():
    return jsonify({
        'active_sessions': len(conversations),
        'total_conversations': sum(len(conv.conversation_history) for conv in conversations.values()) // 2
    })


if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))

    if os.environ.get("DATABASE_URL"):
        try:
            database.init_db()
        except Exception as _db_exc:
            print(f"⚠️  Database init failed: {_db_exc}")
    else:
        print("⚠️  DATABASE_URL not set — running without persistence")

    if not api_key:
        print("\n" + "="*60)
        print("⚠️  WARNING: ANTHROPIC_API_KEY not set!")
        print("="*60)

    print("\n" + "="*60)
    print("🪰 Drosophila Research Assistant - WITH RESEARCH PLANNER")
    print("="*60)
    print(f"\n✓ Server starting on http://localhost:{port}")
    print("\n✨ Features:")
    print("  📚 PubMed + FlyBase literature search")
    print("  🔬 Research Proposal Generator")
    print("  📄 .docx export for proposals")
    print("  🤖 Claude AI (Sonnet 4)")
    print("="*60 + "\n")

    app.run(host='0.0.0.0', port=port, debug=False)