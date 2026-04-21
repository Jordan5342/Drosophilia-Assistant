"""
Persistent storage layer for the Drosophila Research Assistant.
Backed by Supabase PostgreSQL via psycopg2. DATABASE_URL must be set in env.

All JSON fields are stored as JSONB in Postgres. psycopg2 returns JSONB columns as
already-deserialized Python objects (dict/list), so callers never need to call
json.loads() on the results.
"""

import json
import os
from typing import Dict, List, Optional

import psycopg2
import psycopg2.pool
import psycopg2.extras  # registers UUID/JSON adapters

# Module-level connection pool — created once, shared across threads.
_pool: Optional[psycopg2.pool.ThreadedConnectionPool] = None


def _get_pool() -> psycopg2.pool.ThreadedConnectionPool:
    global _pool
    if _pool is None:
        database_url = os.environ.get("DATABASE_URL", "")
        if not database_url:
            raise RuntimeError("DATABASE_URL environment variable is not set.")
        # Supabase / Heroku often use postgres:// but psycopg2 requires postgresql://
        if database_url.startswith("postgres://"):
            database_url = database_url.replace("postgres://", "postgresql://", 1)
        _pool = psycopg2.pool.ThreadedConnectionPool(1, 5, dsn=database_url)
    return _pool


def _acquire():
    return _get_pool().getconn()


def _release(conn):
    _get_pool().putconn(conn)


# ── Schema ────────────────────────────────────────────────────────────────────

def init_db() -> None:
    """Create all tables if they don't already exist."""
    conn = _acquire()
    try:
        with conn:
            with conn.cursor() as cur:
                cur.execute("""
                    CREATE TABLE IF NOT EXISTS sessions (
                        session_id            TEXT PRIMARY KEY,
                        last_topic            TEXT    NOT NULL DEFAULT '',
                        last_genes            JSONB   NOT NULL DEFAULT '[]',
                        planning_mode         BOOLEAN NOT NULL DEFAULT FALSE,
                        awaiting_clarification BOOLEAN NOT NULL DEFAULT FALSE,
                        updated_at            TIMESTAMPTZ NOT NULL DEFAULT NOW()
                    )
                """)
                cur.execute("""
                    CREATE TABLE IF NOT EXISTS conversations (
                        session_id  TEXT PRIMARY KEY,
                        history     JSONB NOT NULL DEFAULT '[]',
                        updated_at  TIMESTAMPTZ NOT NULL DEFAULT NOW()
                    )
                """)
                cur.execute("""
                    CREATE TABLE IF NOT EXISTS proposals (
                        session_id        TEXT PRIMARY KEY,
                        proposal          JSONB,
                        proposal_context  JSONB NOT NULL DEFAULT '{}',
                        updated_at        TIMESTAMPTZ NOT NULL DEFAULT NOW()
                    )
                """)
                cur.execute("""
                    CREATE TABLE IF NOT EXISTS papers (
                        session_id  TEXT PRIMARY KEY,
                        papers      JSONB NOT NULL DEFAULT '[]',
                        updated_at  TIMESTAMPTZ NOT NULL DEFAULT NOW()
                    )
                """)
        print("✅ Database tables ready")
    except Exception as exc:
        print(f"⚠️  init_db error: {exc}")
        raise
    finally:
        _release(conn)


# ── Session state ─────────────────────────────────────────────────────────────

def load_session(session_id: str) -> Optional[Dict]:
    """
    Returns {'last_topic', 'last_genes', 'planning_mode', 'awaiting_clarification'}
    or None if this session has never been saved.
    """
    conn = _acquire()
    try:
        with conn.cursor() as cur:
            cur.execute(
                """SELECT last_topic, last_genes, planning_mode, awaiting_clarification
                   FROM sessions WHERE session_id = %s""",
                (session_id,)
            )
            row = cur.fetchone()
        if row is None:
            return None
        last_topic, last_genes, planning_mode, awaiting_clarification = row
        return {
            "last_topic": last_topic or "",
            # psycopg2 deserializes JSONB automatically; guard against None just in case
            "last_genes": (last_genes if isinstance(last_genes, list) else []),
            "planning_mode": bool(planning_mode),
            "awaiting_clarification": bool(awaiting_clarification),
        }
    finally:
        _release(conn)


def save_session(
    session_id: str,
    last_topic: str,
    last_genes: List[str],
    planning_mode: bool,
    awaiting_clarification: bool,
) -> None:
    # Coerce to list — last_genes could theoretically be a set in edge cases
    genes_list = list(last_genes) if last_genes else []
    conn = _acquire()
    try:
        with conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    INSERT INTO sessions
                        (session_id, last_topic, last_genes, planning_mode,
                         awaiting_clarification, updated_at)
                    VALUES (%s, %s, %s::jsonb, %s, %s, NOW())
                    ON CONFLICT (session_id) DO UPDATE SET
                        last_topic             = EXCLUDED.last_topic,
                        last_genes             = EXCLUDED.last_genes,
                        planning_mode          = EXCLUDED.planning_mode,
                        awaiting_clarification = EXCLUDED.awaiting_clarification,
                        updated_at             = NOW()
                    """,
                    (session_id, last_topic, json.dumps(genes_list),
                     planning_mode, awaiting_clarification)
                )
    finally:
        _release(conn)


def delete_session(session_id: str) -> None:
    """Remove all persisted data for a session (used by /api/reset)."""
    conn = _acquire()
    try:
        with conn:
            with conn.cursor() as cur:
                for table in ("sessions", "conversations", "proposals", "papers"):
                    cur.execute(
                        f"DELETE FROM {table} WHERE session_id = %s",
                        (session_id,)
                    )
    finally:
        _release(conn)


# ── Conversation history ──────────────────────────────────────────────────────

def load_conversation(session_id: str) -> List[Dict]:
    """Returns the conversation history list, or [] if none stored."""
    conn = _acquire()
    try:
        with conn.cursor() as cur:
            cur.execute(
                "SELECT history FROM conversations WHERE session_id = %s",
                (session_id,)
            )
            row = cur.fetchone()
        if row is None or row[0] is None:
            return []
        history = row[0]
        return history if isinstance(history, list) else []
    finally:
        _release(conn)


def save_conversation(session_id: str, history: List[Dict]) -> None:
    """Persist the full conversation history list (replaces any existing)."""
    conn = _acquire()
    try:
        with conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    INSERT INTO conversations (session_id, history, updated_at)
                    VALUES (%s, %s::jsonb, NOW())
                    ON CONFLICT (session_id) DO UPDATE SET
                        history    = EXCLUDED.history,
                        updated_at = NOW()
                    """,
                    (session_id, json.dumps(history))
                )
    finally:
        _release(conn)


def append_message(session_id: str, role: str, content: str) -> None:
    """Append a single message to the stored conversation history."""
    history = load_conversation(session_id)
    history.append({"role": role, "content": content})
    save_conversation(session_id, history)


def clear_conversation(session_id: str) -> None:
    """Reset the stored conversation history to an empty list."""
    conn = _acquire()
    try:
        with conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    INSERT INTO conversations (session_id, history, updated_at)
                    VALUES (%s, '[]'::jsonb, NOW())
                    ON CONFLICT (session_id) DO UPDATE SET
                        history    = '[]'::jsonb,
                        updated_at = NOW()
                    """,
                    (session_id,)
                )
    finally:
        _release(conn)


# ── Proposals ─────────────────────────────────────────────────────────────────

def save_proposal(
    session_id: str,
    proposal: Optional[Dict],
    proposal_context: Optional[Dict] = None,
) -> None:
    conn = _acquire()
    try:
        with conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    INSERT INTO proposals
                        (session_id, proposal, proposal_context, updated_at)
                    VALUES (%s, %s::jsonb, %s::jsonb, NOW())
                    ON CONFLICT (session_id) DO UPDATE SET
                        proposal         = EXCLUDED.proposal,
                        proposal_context = EXCLUDED.proposal_context,
                        updated_at       = NOW()
                    """,
                    (
                        session_id,
                        json.dumps(proposal) if proposal is not None else None,
                        json.dumps(proposal_context or {}),
                    )
                )
    finally:
        _release(conn)


def load_proposal(session_id: str) -> Optional[Dict]:
    """
    Returns {'proposal': dict|None, 'proposal_context': dict}
    or None if nothing has been saved for this session.
    """
    conn = _acquire()
    try:
        with conn.cursor() as cur:
            cur.execute(
                "SELECT proposal, proposal_context FROM proposals WHERE session_id = %s",
                (session_id,)
            )
            row = cur.fetchone()
        if row is None:
            return None
        proposal, proposal_context = row
        return {
            "proposal": proposal,  # already a dict (psycopg2 JSONB auto-deserialization)
            "proposal_context": (proposal_context if isinstance(proposal_context, dict) else {}),
        }
    finally:
        _release(conn)


# ── Papers ────────────────────────────────────────────────────────────────────

def save_papers(session_id: str, papers: List[Dict]) -> None:
    conn = _acquire()
    try:
        with conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    INSERT INTO papers (session_id, papers, updated_at)
                    VALUES (%s, %s::jsonb, NOW())
                    ON CONFLICT (session_id) DO UPDATE SET
                        papers     = EXCLUDED.papers,
                        updated_at = NOW()
                    """,
                    (session_id, json.dumps(papers))
                )
    finally:
        _release(conn)


def load_papers(session_id: str) -> List[Dict]:
    """Returns the last_papers list for a session, or [] if none stored."""
    conn = _acquire()
    try:
        with conn.cursor() as cur:
            cur.execute(
                "SELECT papers FROM papers WHERE session_id = %s",
                (session_id,)
            )
            row = cur.fetchone()
        if row is None or row[0] is None:
            return []
        papers = row[0]
        return papers if isinstance(papers, list) else []
    finally:
        _release(conn)
