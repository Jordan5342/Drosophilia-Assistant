"""
Batch-run the lab's question list through the literature -> hypothesis agents on
several Claude models, so nobody has to paste 50 questions into the app one by one.

Questions live in lab_questions.md. Each question gets one PubMed search, and every
model sees the same papers, so the models can be compared fairly. Results are saved
after every finished run, and a rerun only does the question/model pairs that are new,
edited, or failed last time.

Usage:
    python run_questions.py                          # all questions, all models
    python run_questions.py --models claude-sonnet-4-6,claude-opus-5-5
    python run_questions.py --only 26,40-42          # just these questions
    python run_questions.py --only 26 --redo         # rerun even if unchanged

Writes question_runs/results.json (full output) and question_runs/results.md
(readable). Requires ANTHROPIC_API_KEY; NCBI_API_KEY is optional.
"""

import argparse
import hashlib
import json
import os
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

import anthropic

from agents import CLAUDE_MODELS
from agents.hypothesis_agent import HypothesisAgent
from agents.literature_agent import LiteratureAgent
from drosophila_assistant import ANTHROPIC_TIMEOUT_SECONDS, DrosophilaAssistant

QUESTION_RE = re.compile(r"^(\d+)\.\s+(.+)$")
OVERRIDE_RE = re.compile(r"^\s+(search|genes):\s*(.*)$", re.IGNORECASE)


def load_questions(path: Path) -> dict:
    """Parse numbered questions plus optional indented search:/genes: overrides."""
    text = re.sub(r"<!--.*?-->", "", path.read_text(), flags=re.DOTALL)
    questions, current = {}, None
    for line in text.splitlines():
        match = QUESTION_RE.match(line.strip()) if not line.startswith((" ", "\t")) else None
        if match:
            current = {"question": match.group(2).strip(), "search": None, "genes": None}
            questions[int(match.group(1))] = current
            continue
        override = OVERRIDE_RE.match(line)
        if override and current:
            key, value = override.group(1).lower(), override.group(2).strip()
            current[key] = [g.strip() for g in value.split(",") if g.strip()] if key == "genes" else value
    return questions


def parse_numbers(spec: str) -> set:
    numbers = set()
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            lo, hi = part.split("-")
            numbers.update(range(int(lo), int(hi) + 1))
        elif part:
            numbers.add(int(part))
    return numbers


def fingerprint(q: dict) -> str:
    """Changes when the question or its search overrides are edited."""
    return hashlib.sha1(json.dumps([q["question"], q["search"], q["genes"]]).encode()).hexdigest()[:12]


def is_good(entry: dict) -> bool:
    if not entry or "_error" in entry or entry.get("refused"):
        return False
    return not any(entry.get(k, {}).get("_parse_error") for k in ("literature", "hypothesis"))


class Runner:
    def __init__(self, out_dir: Path, models: list):
        self.out_dir = out_dir
        self.models = models
        self.results_path = out_dir / "results.json"
        self.papers_path = out_dir / "papers.json"
        self.results = json.loads(self.results_path.read_text()) if self.results_path.exists() else {}
        self.paper_cache = json.loads(self.papers_path.read_text()) if self.papers_path.exists() else {}
        self.assistant = DrosophilaAssistant(api_key=os.environ["ANTHROPIC_API_KEY"])
        # PubMed allows only a few requests per second, so searches run one at a time;
        # the Claude calls are what run in parallel.
        self.pubmed_lock = threading.Lock()
        self.save_lock = threading.Lock()

    def papers_for(self, number: int, q: dict) -> dict:
        key = f"{number}:{fingerprint(q)}"
        with self.pubmed_lock:
            if key not in self.paper_cache:
                topic = q["search"] or self.assistant.build_chat_search_query(q["question"])
                genes = q["genes"] if q["genes"] is not None else self.assistant.extract_gene_names(q["question"])
                papers = self.assistant.fetch_literature_for_proposal(topic=topic, genes=genes, user_clarifications="")
                self.paper_cache[key] = {"search_topic": topic, "search_genes": genes, "papers": papers}
                with self.save_lock:
                    self.papers_path.write_text(json.dumps(self.paper_cache, indent=2))
            return self.paper_cache[key]

    def run_one(self, number: int, q: dict, model: str) -> dict:
        search = self.papers_for(number, q)
        client = anthropic.Anthropic(timeout=ANTHROPIC_TIMEOUT_SECONDS)
        stops = []
        create = client.messages.create

        def recording_create(**kwargs):
            response = create(**kwargs)
            stops.append(response.stop_reason)
            return response

        client.messages.create = recording_create
        start = time.time()
        try:
            literature = LiteratureAgent(client, model).run(q["question"], search["papers"])
            hypothesis = HypothesisAgent(client, model).run(literature, q["question"])
            entry = {"literature": literature, "hypothesis": hypothesis}
        except Exception as exc:
            entry = {"_error": f"{type(exc).__name__}: {exc}"}
        entry.update(
            question_number=number, question=q["question"], model=model, fingerprint=fingerprint(q),
            search_topic=search["search_topic"], search_genes=search["search_genes"],
            papers=[{"pmid": p.get("pmid"), "title": p.get("title"), "year": p.get("year")} for p in search["papers"]],
            stop_reasons=stops, refused="refusal" in stops,
            seconds=round(time.time() - start), finished_at=time.strftime("%Y-%m-%d %H:%M:%S"),
        )
        with self.save_lock:
            self.results.setdefault(str(number), {})[model] = entry
            self.results_path.write_text(json.dumps(self.results, indent=2))
        return entry

    def pending(self, questions: dict, redo: bool) -> list:
        jobs = []
        for number, q in sorted(questions.items()):
            for model in self.models:
                old = self.results.get(str(number), {}).get(model)
                if redo or not is_good(old) or old.get("fingerprint") != fingerprint(q):
                    jobs.append((number, q, model))
        return jobs


def status_of(entry: dict) -> str:
    if "_error" in entry:
        return f"ERROR {entry['_error'][:120]}"
    if entry.get("refused"):
        return "REFUSED"
    if not is_good(entry):
        return "unparseable output"
    return "ok"


def write_markdown(path: Path, questions: dict, results: dict, models: list):
    lines = ["# Batch hypothesis runs", "",
             "Unreviewed drafts written by AI. Each model got the same question and the same papers.", ""]
    for number, q in sorted(questions.items()):
        lines += [f"## {number}. {q['question']}", ""]
        runs = results.get(str(number), {})
        stale = [m for m in models if runs.get(m, {}).get("fingerprint") not in (None, fingerprint(q))]
        if not any(m in runs for m in models):
            lines += ["_Not run yet._", ""]
            continue
        first = next(runs[m] for m in models if m in runs)
        lines += [f"PubMed search: `{first['search_topic']}`, genes {first['search_genes'] or 'none'}, "
                  f"{len(first['papers'])} papers", ""]
        for model in models:
            entry = runs.get(model)
            if not entry:
                continue
            label = CLAUDE_MODELS.get(model, model)
            note = " (question edited since this run)" if model in stale else ""
            lines += [f"### {label}{note}", ""]
            status = status_of(entry)
            if status != "ok":
                lines += [f"_{status}_", ""]
                continue
            h = entry["hypothesis"]
            lines += [f"**Hypothesis:** {h.get('hypothesis', '')}", ""]
            if h.get("specific_aims"):
                lines += ["**Aims:**"] + [f"- {a}" for a in h["specific_aims"]] + [""]
            for field, title in (("experimental_approach", "Approach"), ("expected_outcomes", "Expected outcomes"),
                                 ("novelty_claim", "Novelty"), ("feasibility_notes", "Feasibility")):
                if h.get(field):
                    lines += [f"**{title}:** {h[field]}", ""]
    path.write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description="Batch-run lab questions through the hypothesis pipeline.")
    parser.add_argument("--questions", default="lab_questions.md", help="Question file (default lab_questions.md)")
    parser.add_argument("--models", default=",".join(CLAUDE_MODELS), help="Comma-separated Claude model IDs")
    parser.add_argument("--only", help="Question numbers to run, e.g. 26,40-42")
    parser.add_argument("--redo", action="store_true", help="Rerun the selected questions even if unchanged")
    parser.add_argument("--workers", type=int, default=4, help="Claude runs in parallel (default 4)")
    parser.add_argument("--out-dir", default="question_runs", help="Output folder (default question_runs)")
    args = parser.parse_args()

    if not os.environ.get("ANTHROPIC_API_KEY"):
        sys.exit("❌ ANTHROPIC_API_KEY not set. Add it to your .env file.")
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    unknown = [m for m in models if m not in CLAUDE_MODELS]
    if unknown:
        sys.exit(f"❌ Unknown model(s) {unknown}. Choose from {list(CLAUDE_MODELS)}")

    questions = load_questions(Path(args.questions))
    selected = {n: q for n, q in questions.items() if not args.only or n in parse_numbers(args.only)}
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)
    runner = Runner(out_dir, models)
    jobs = runner.pending(selected, args.redo)
    print(f"🧪 {len(selected)} questions × {len(models)} models: {len(jobs)} to run, "
          f"{len(selected) * len(models) - len(jobs)} already done")

    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(runner.run_one, *job): job for job in jobs}
        for future in as_completed(futures):
            number, _, model = futures[future]
            done += 1
            try:
                entry = future.result()
                print(f"[{done}/{len(jobs)}] Q{number} {model}: {status_of(entry)} ({entry['seconds']}s)", flush=True)
            except Exception as exc:
                print(f"[{done}/{len(jobs)}] Q{number} {model}: ERROR {exc}", flush=True)
            write_markdown(out_dir / "results.md", questions, runner.results, models)

    write_markdown(out_dir / "results.md", questions, runner.results, models)
    print(f"✅ Done. Results in {(out_dir / 'results.md').resolve()}")


if __name__ == "__main__":
    main()
