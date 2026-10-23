import json
import os
import re
from typing import Dict

# Default Claude model for chat, planning, the pipeline agents and the evidence check.
# Override with the CLAUDE_MODEL env var; the app's model menu overrides it per request.
CLAUDE_MODEL = os.environ.get("CLAUDE_MODEL", "claude-sonnet-4-6")

# Models the app's model menu may request, so a request can't name an arbitrary model.
CLAUDE_MODELS = {
    "claude-sonnet-4-6": "Sonnet 4.6",
    "claude-sonnet-5-5": "Sonnet 5.5",
    "claude-opus-5-5": "Opus 5.5",
}

# These models reject a forced tool_choice ("tool"/"any") with a 400.
NO_FORCED_TOOL_CHOICE = {"claude-sonnet-5-5", "claude-opus-5-5"}


def resolve_model(requested) -> str:
    """The requested model if it is allowed, otherwise the default."""
    if requested in CLAUDE_MODELS or requested == CLAUDE_MODEL:
        return requested
    return CLAUDE_MODEL


def response_text(response) -> str:
    """Join a response's text blocks; newer models can return a thinking block first."""
    return "".join(
        getattr(block, "text", "")
        for block in response.content
        if getattr(block, "type", "text") == "text"
    )


def parse_json_response(raw: str, agent_name: str) -> Dict:
    """Extract and parse JSON from an agent response, with fallback on malformed output."""
    text = re.sub(r'^```(?:json)?\s*', '', raw, flags=re.MULTILINE)
    text = re.sub(r'\s*```$', '', text, flags=re.MULTILINE)
    text = text.strip()

    try:
        return json.loads(text)
    except json.JSONDecodeError:
        match = re.search(r'\{[\s\S]*\}', text)
        if match:
            try:
                return json.loads(match.group(0))
            except json.JSONDecodeError:
                pass

    return {"_parse_error": True, "raw": raw[:500], "agent": agent_name}
