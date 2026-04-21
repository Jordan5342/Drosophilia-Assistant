import json
import re
from typing import Dict


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
