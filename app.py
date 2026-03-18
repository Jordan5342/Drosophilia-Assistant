from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import os
from datetime import datetime
import json
import signal
from drosophila_assistant import DrosophilaAssistant

app = Flask(__name__, static_folder='static')
CORS(app)

api_key = os.environ.get("ANTHROPIC_API_KEY")
assistant = DrosophilaAssistant(api_key) if api_key else None
conversations = {}


class TimeoutException(Exception):
    pass


def timeout_handler(signum, frame):
    raise TimeoutException()


def get_session(session_id: str) -> DrosophilaAssistant:
    if session_id not in conversations:
        conversations[session_id] = DrosophilaAssistant(api_key)
    return conversations[session_id]


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

            return jsonify({
                'response': result['response'],
                'is_planning': result['is_planning'],
                'ready_for_export': result['ready_for_export'],
                'has_proposal': result['proposal'] is not None,
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

        if session_id not in conversations:
            return jsonify({'error': 'No active session found'}), 404

        session_assistant = conversations[session_id]
        proposal = session_assistant.planner.get_proposal_for_export()

        if not proposal:
            return jsonify({'error': 'No proposal to export. Generate a proposal first.'}), 404

        # Write the proposal JSON to a temp file and run the Node docx generator
        import subprocess
        import tempfile

        proposal_json = json.dumps(proposal)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_path = f"/tmp/proposal_{timestamp}.docx"

        # Write Node script inline
        node_script = f"""
const {{ Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
        LevelFormat, ExternalHyperlink, BorderStyle }} = require('docx');
const fs = require('fs');

const proposal = {proposal_json};

const BLUE = "2E5FA3";
const DARK = "1A1A2E";
const GRAY = "555555";

function makeHeading1(text) {{
  return new Paragraph({{
    heading: HeadingLevel.HEADING_1,
    spacing: {{ before: 360, after: 120 }},
    border: {{ bottom: {{ style: BorderStyle.SINGLE, size: 6, color: BLUE, space: 4 }} }},
    children: [new TextRun({{ text, bold: true, color: BLUE, size: 28, font: "Arial" }})]
  }});
}}

function makeHeading2(text) {{
  return new Paragraph({{
    heading: HeadingLevel.HEADING_2,
    spacing: {{ before: 240, after: 80 }},
    children: [new TextRun({{ text, bold: true, color: DARK, size: 24, font: "Arial" }})]
  }});
}}

function makeBody(text, options = {{}}) {{
  const runs = [];
  // Handle bold markers **text**
  const parts = (text || '').split(/\\*\\*([^*]+)\\*\\*/g);
  parts.forEach((part, i) => {{
    runs.push(new TextRun({{
      text: part,
      bold: i % 2 === 1,
      size: 22,
      font: "Arial",
      color: GRAY,
      ...options
    }}));
  }});
  return new Paragraph({{
    spacing: {{ before: 60, after: 120 }},
    children: runs
  }});
}}

function makeBullet(text, ref) {{
  return new Paragraph({{
    numbering: {{ reference: ref, level: 0 }},
    spacing: {{ before: 40, after: 40 }},
    children: [new TextRun({{ text: text || '', size: 22, font: "Arial", color: GRAY }})]
  }});
}}

function spacer() {{
  return new Paragraph({{ spacing: {{ before: 60, after: 60 }}, children: [new TextRun('')] }});
}}

const children = [];

// Title page
children.push(new Paragraph({{
  alignment: AlignmentType.CENTER,
  spacing: {{ before: 480, after: 240 }},
  children: [new TextRun({{ text: proposal.title || 'Research Proposal', bold: true, size: 40, font: "Arial", color: DARK }})]
}}));

children.push(new Paragraph({{
  alignment: AlignmentType.CENTER,
  spacing: {{ before: 0, after: 120 }},
  children: [new TextRun({{ text: 'Drosophila Research Proposal', size: 24, font: "Arial", color: BLUE }})]
}}));

children.push(new Paragraph({{
  alignment: AlignmentType.CENTER,
  spacing: {{ before: 0, after: 480 }},
  children: [new TextRun({{ text: new Date().toLocaleDateString('en-US', {{ year: 'numeric', month: 'long', day: 'numeric' }}), size: 22, font: "Arial", color: GRAY }})]
}}));

children.push(spacer());

// Background
children.push(makeHeading1('Background & Significance'));
children.push(makeBody(proposal.background || ''));
children.push(spacer());

// Hypothesis
children.push(makeHeading1('Central Hypothesis'));
children.push(new Paragraph({{
  spacing: {{ before: 80, after: 80 }},
  children: [
    new TextRun({{ text: 'H\u2081: ', bold: true, size: 22, font: "Arial", color: BLUE }}),
    new TextRun({{ text: proposal.central_hypothesis || '', size: 22, font: "Arial", color: DARK }})
  ]
}}));

if (proposal.null_hypothesis) {{
  children.push(new Paragraph({{
    spacing: {{ before: 80, after: 80 }},
    children: [
      new TextRun({{ text: 'H\u2080: ', bold: true, size: 22, font: "Arial", color: GRAY }}),
      new TextRun({{ text: proposal.null_hypothesis, size: 22, font: "Arial", color: GRAY }})
    ]
  }}));
}}

if (proposal.rationale) {{
  children.push(spacer());
  children.push(makeHeading2('Rationale'));
  children.push(makeBody(proposal.rationale));
}}
children.push(spacer());

// Specific Aims
if (proposal.specific_aims && proposal.specific_aims.length > 0) {{
  children.push(makeHeading1('Specific Aims'));
  proposal.specific_aims.forEach(aim => {{
    children.push(makeHeading2(`Aim ${{aim.aim_number}}: ${{aim.title}}`));
    if (aim.objective) {{
      children.push(new Paragraph({{
        spacing: {{ before: 60, after: 40 }},
        children: [
          new TextRun({{ text: 'Objective: ', bold: true, size: 22, font: "Arial", color: DARK }}),
          new TextRun({{ text: aim.objective, size: 22, font: "Arial", color: GRAY }})
        ]
      }}));
    }}
    if (aim.approach) {{
      children.push(new Paragraph({{
        spacing: {{ before: 40, after: 40 }},
        children: [
          new TextRun({{ text: 'Approach: ', bold: true, size: 22, font: "Arial", color: DARK }}),
          new TextRun({{ text: aim.approach, size: 22, font: "Arial", color: GRAY }})
        ]
      }}));
    }}
    if (aim.expected_outcomes) {{
      children.push(new Paragraph({{
        spacing: {{ before: 40, after: 40 }},
        children: [
          new TextRun({{ text: 'Expected Outcomes: ', bold: true, size: 22, font: "Arial", color: DARK }}),
          new TextRun({{ text: aim.expected_outcomes, size: 22, font: "Arial", color: GRAY }})
        ]
      }}));
    }}
    if (aim.potential_pitfalls) {{
      children.push(new Paragraph({{
        spacing: {{ before: 40, after: 40 }},
        children: [
          new TextRun({{ text: 'Potential Pitfalls: ', bold: true, size: 22, font: "Arial", color: DARK }}),
          new TextRun({{ text: aim.potential_pitfalls, size: 22, font: "Arial", color: GRAY }})
        ]
      }}));
    }}
    if (aim.mitigation) {{
      children.push(new Paragraph({{
        spacing: {{ before: 40, after: 80 }},
        children: [
          new TextRun({{ text: 'Mitigation: ', bold: true, size: 22, font: "Arial", color: DARK }}),
          new TextRun({{ text: aim.mitigation, size: 22, font: "Arial", color: GRAY }})
        ]
      }}));
    }}
    children.push(spacer());
  }});
}}

// Experimental Approach
if (proposal.experimental_approach) {{
  children.push(makeHeading1('Experimental Approach'));
  children.push(makeBody(proposal.experimental_approach));
  children.push(spacer());
}}

// Controls
if (proposal.controls) {{
  children.push(makeHeading1('Controls'));
  children.push(makeBody(proposal.controls));
  children.push(spacer());
}}

// Literature Gaps
if (proposal.literature_gaps && proposal.literature_gaps.length > 0) {{
  children.push(makeHeading1('\u26a0\ufe0f Literature Gaps & Research Opportunities'));
  proposal.literature_gaps.forEach(gap => {{
    children.push(makeBullet(gap, 'bullets'));
  }});
  children.push(spacer());
}}

// Timeline
if (proposal.timeline && proposal.timeline.length > 0) {{
  children.push(makeHeading1('Timeline'));
  proposal.timeline.forEach(phase => {{
    children.push(new Paragraph({{
      spacing: {{ before: 80, after: 40 }},
      children: [
        new TextRun({{ text: phase.phase + ': ', bold: true, size: 22, font: "Arial", color: BLUE }}),
        new TextRun({{ text: phase.milestones, size: 22, font: "Arial", color: GRAY }})
      ]
    }}));
  }});
  children.push(spacer());
}}

// References
if (proposal.references && proposal.references.length > 0) {{
  children.push(makeHeading1('References'));
  proposal.references.forEach(ref => {{
    const refChildren = [
      new TextRun({{ text: '[' + ref.number + '] ', bold: true, size: 20, font: "Arial", color: BLUE }}),
      new TextRun({{ text: ref.citation, size: 20, font: "Arial", color: GRAY }})
    ];
    if (ref.url) {{
      refChildren.push(new TextRun({{ text: ' ', size: 20 }}));
      refChildren.push(new ExternalHyperlink({{
        link: ref.url,
        children: [new TextRun({{ text: '[Link]', style: 'Hyperlink', size: 20, font: "Arial" }})]
      }}));
    }}
    children.push(new Paragraph({{
      spacing: {{ before: 40, after: 40 }},
      children: refChildren
    }}));
  }});
}}

const doc = new Document({{
  numbering: {{
    config: [
      {{ reference: "bullets",
        levels: [{{ level: 0, format: LevelFormat.BULLET, text: "\u2022", alignment: AlignmentType.LEFT,
          style: {{ paragraph: {{ indent: {{ left: 720, hanging: 360 }} }} }} }}] }}
    ]
  }},
  styles: {{
    default: {{ document: {{ run: {{ font: "Arial", size: 22 }} }} }},
    paragraphStyles: [
      {{ id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: {{ size: 28, bold: true, font: "Arial", color: BLUE }},
        paragraph: {{ spacing: {{ before: 360, after: 120 }}, outlineLevel: 0 }} }},
      {{ id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: {{ size: 24, bold: true, font: "Arial", color: DARK }},
        paragraph: {{ spacing: {{ before: 240, after: 80 }}, outlineLevel: 1 }} }},
    ]
  }},
  sections: [{{
    properties: {{
      page: {{
        size: {{ width: 12240, height: 15840 }},
        margin: {{ top: 1440, right: 1440, bottom: 1440, left: 1440 }}
      }}
    }},
    children
  }}]
}});

Packer.toBuffer(doc).then(buffer => {{
  fs.writeFileSync('{output_path}', buffer);
  console.log('SUCCESS:{output_path}');
}}).catch(err => {{
  console.error('ERROR:' + err.message);
  process.exit(1);
}});
"""

        # Ensure docx npm package is available — install to /tmp if needed
        node_modules_path = "/tmp/node_modules"
        docx_check = os.path.join(node_modules_path, "docx")
        if not os.path.exists(docx_check):
            print("  Installing docx npm package to /tmp...")
            install_result = subprocess.run(
                ['npm', 'install', '--prefix', '/tmp', 'docx'],
                capture_output=True, text=True, timeout=60
            )
            if install_result.returncode != 0:
                print(f"npm install error: {install_result.stderr}")
                return jsonify({'error': f'Failed to install docx package: {install_result.stderr}'}), 500
            print("  docx installed")

        # Write script
        script_path = f"/tmp/gen_proposal_{timestamp}.js"
        with open(script_path, 'w') as f:
            f.write(node_script)

        # Run with NODE_PATH pointing to /tmp/node_modules
        env = os.environ.copy()
        env['NODE_PATH'] = node_modules_path

        result = subprocess.run(
            ['node', script_path],
            capture_output=True, text=True, timeout=30,
            env=env
        )

        os.unlink(script_path)

        if result.returncode != 0 or not os.path.exists(output_path):
            print(f"Node error: {result.stderr}")
            return jsonify({'error': f'Failed to generate .docx: {result.stderr}'}), 500

        # Read and return the file as base64
        import base64
        with open(output_path, 'rb') as f:
            docx_bytes = f.read()
        os.unlink(output_path)

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

        if session_id not in conversations:
            return jsonify({'error': 'No active session found'}), 404

        session_assistant = conversations[session_id]
        proposal = session_assistant.planner.get_proposal_for_export()

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