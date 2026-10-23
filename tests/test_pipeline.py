import importlib
import json
import os
import sys
import types
import unittest
from unittest.mock import Mock, patch


def _import_agent_modules():
    anthropic_stub = types.ModuleType("anthropic")
    anthropic_stub.Anthropic = object

    google_stub = types.ModuleType("google")
    google_stub.__path__ = []
    genai_stub = types.ModuleType("google.genai")
    genai_stub.__path__ = []
    genai_stub.Client = object
    genai_types_stub = types.ModuleType("google.genai.types")
    google_stub.genai = genai_stub

    stubs = {
        "anthropic": anthropic_stub,
        "google": google_stub,
        "google.genai": genai_stub,
        "google.genai.types": genai_types_stub,
    }
    with patch.dict(sys.modules, stubs):
        pipeline_module = importlib.import_module("agents.pipeline")
        planner_module = importlib.import_module("research_planner")
    return pipeline_module, planner_module


pipeline_module, planner_module = _import_agent_modules()
verifier_module = importlib.import_module("evidence_verifier")


class PipelineTests(unittest.TestCase):
    def run_pipeline(self, hypothesis_run=None, critic_run=None):
        retrieved = [{
            "pmid": "12345678",
            "title": "Verified paper title",
            "authors": "Author A et al.",
            "year": "2024",
            "url": "https://pubmed.ncbi.nlm.nih.gov/12345678/",
            "abstract": "Abstract text",
        }]
        literature_agent = Mock()
        literature_agent.run.return_value = {
            "genes": ["foxo"],
            "key_papers": [
                {"pmid": "12345678", "title": "Invented title", "url": "https://bad.example"},
                {"pmid": "99999999", "title": "Unretrieved paper"},
            ],
            "summary": "Summary",
            "identified_gaps": [],
        }
        hypothesis_agent = Mock()
        if hypothesis_run is None:
            hypothesis_agent.run.return_value = {
                "hypothesis": "A falsifiable hypothesis",
                "specific_aims": [],
            }
        else:
            hypothesis_agent.run.side_effect = hypothesis_run
        critic = Mock()
        critic.run.return_value = {
            "scores": {
                "hallucination_defense": 4,
                "genetic_feasibility": 3,
                "control_reasoning": 5,
            },
            "issues": [],
            "suggestions": [],
            "reasoning": "Scores are descriptive.",
        }
        if critic_run is not None:
            critic.run.side_effect = critic_run
        fetch_literature = Mock(return_value=retrieved)

        with (
            patch.object(pipeline_module, "LiteratureAgent", return_value=literature_agent),
            patch.object(pipeline_module, "HypothesisAgent", return_value=hypothesis_agent),
            patch.dict(os.environ, {"GEMINI_API_KEY": "test-gemini-key"}),
            patch.object(pipeline_module, "GeminiCriticAgent", return_value=critic) as critic_factory,
        ):
            pipeline = pipeline_module.AgentPipeline(
                client=object(),
                fetch_literature_fn=fetch_literature
            )
            events = list(pipeline.run_streaming(
                query="test query",
                topic="test topic",
                genes=["foxo"]
            ))
        critic_factory.assert_called_once_with(api_key="test-gemini-key")
        fetch_literature.assert_called_once_with(
            topic="test topic",
            genes=["foxo"],
            user_clarifications="test query"
        )
        return events, literature_agent, hypothesis_agent, critic

    def test_pipeline_runs_each_agent_once_and_verifies_literature_sources(self):
        events, literature_agent, hypothesis_agent, critic = self.run_pipeline()

        started_agents = [
            event["data"]["agent"]
            for event in events
            if event["type"] == "agent_start"
        ]
        self.assertEqual(started_agents, ["literature", "hypothesis", "critic"])
        literature_agent.run.assert_called_once()
        hypothesis_agent.run.assert_called_once()
        critic.run.assert_called_once()

        complete = next(event["data"] for event in events if event["type"] == "pipeline_complete")
        self.assertIsNone(complete["passed"])
        self.assertEqual(complete["status"], "complete")
        self.assertEqual(complete["best_score"], 12)
        self.assertEqual(
            complete["final_literature"]["key_papers"],
            [{
                "pmid": "12345678",
                "title": "Verified paper title",
                "authors": "Author A et al.",
                "year": "2024",
                "url": "https://pubmed.ncbi.nlm.nih.gov/12345678/",
                "abstract": "Abstract text",
            }]
        )
        routing = next(event["data"] for event in events if event["type"] == "routing")
        self.assertIn("no pass/fail verdict", routing["message"])

    def test_agent_exception_is_reported_as_partial_not_success(self):
        events, _, hypothesis_agent, _ = self.run_pipeline(
            hypothesis_run=RuntimeError("model unavailable")
        )

        self.assertEqual(hypothesis_agent.run.call_count, 1)
        hypothesis_event = next(
            event["data"] for event in events
            if event["type"] == "agent_complete" and event["data"]["agent"] == "hypothesis"
        )
        self.assertEqual(hypothesis_event["status"], "error")
        self.assertEqual(hypothesis_event["error"], "model unavailable")

        complete = next(event["data"] for event in events if event["type"] == "pipeline_complete")
        self.assertEqual(complete["status"], "partial")
        self.assertEqual(complete["failed_agents"], ["hypothesis"])
        self.assertIsNone(complete["passed"])

    def test_critic_exception_does_not_fabricate_a_score_or_verdict(self):
        events, _, _, critic = self.run_pipeline(
            critic_run=RuntimeError("critic unavailable")
        )

        self.assertEqual(critic.run.call_count, 1)
        critic_event = next(
            event["data"] for event in events
            if event["type"] == "agent_complete" and event["data"]["agent"] == "critic"
        )
        self.assertEqual(critic_event["status"], "error")
        self.assertEqual(critic_event["scores"], {})
        self.assertIsNone(critic_event["verdict"])

        complete = next(event["data"] for event in events if event["type"] == "pipeline_complete")
        self.assertEqual(complete["status"], "partial")
        self.assertEqual(complete["failed_agents"], ["critic"])
        self.assertIsNone(complete["passed"])
        self.assertEqual(complete["best_score"], 0)


class GeminiCriticTests(unittest.TestCase):
    def test_failed_reviewers_get_no_score_instead_of_default(self):
        critic_class = pipeline_module.GeminiCriticAgent
        fake_types = types.SimpleNamespace(
            GenerateContentConfig=lambda **kwargs: kwargs,
            ThinkingConfig=lambda **kwargs: kwargs,
        )
        responses = {
            "hallucination_defense": '{"score": 4, "issues": [], "suggestions": [], "reasoning": "ok"}',
            "genetic_feasibility": '```json\n{\n  "score": 5,\n  "reasoning": "',  # truncated
        }

        def generate_content(model, contents, config):
            name = next(n for n, p in critic_globals["REVIEWER_PROMPTS"].items() if p == config["system_instruction"])
            if name not in responses:
                raise TimeoutError("request timed out")
            return types.SimpleNamespace(text=responses[name])

        critic_globals = critic_class.run.__globals__
        critic = object.__new__(critic_class)
        critic.client = types.SimpleNamespace(models=types.SimpleNamespace(generate_content=generate_content))
        with patch.dict(critic_globals, {"types": fake_types}):
            result = critic.run({"hypothesis": "h"}, {"summary": "s"}, "topic")

        self.assertEqual(result["scores"], {"hallucination_defense": 4})
        self.assertEqual(result["failed_reviewers"], ["genetic_feasibility", "control_reasoning"])
        self.assertNotIn("_error", result)

        responses.clear()
        with patch.dict(critic_globals, {"types": fake_types}):
            result = critic.run({"hypothesis": "h"}, {"summary": "s"}, "topic")
        self.assertEqual(result["scores"], {})
        self.assertIn("_error", result)


class PipelineProposalTests(unittest.TestCase):
    def test_pipeline_output_maps_to_persistable_proposal_with_traceable_sources(self):
        proposal = pipeline_module.build_pipeline_proposal(
            "sleep regulation",
            {
                "summary": "Prior work links sleep to neural circuits.",
                "pathways": ["circadian clock"],
                "identified_gaps": ["The cell-specific mechanism is unknown."],
                "evidence_claims": [{
                    "claim": "Prior work links sleep to neural circuits.",
                    "evidence_type": "finding",
                    "source_ids": ["12345678"],
                }],
                "key_papers": [{
                    "pmid": "12345678",
                    "title": "Sleep in Drosophila",
                    "authors": "Author A et al.",
                    "year": "2024",
                    "journal": "Fly Journal",
                    "doi": "10.1234/example",
                    "url": "https://pubmed.ncbi.nlm.nih.gov/12345678/",
                    "source": "PubMed",
                    "_relevance_score": 0.42,
                }],
            },
            {
                "hypothesis": "Neural circuit X controls sleep.",
                "specific_aims": ["Map circuit X — Test activity across sleep states"],
                "experimental_approach": "Use genetic labeling and activity assays.",
                "fly_lines_needed": [{"line": "UAS-GCaMP", "purpose": "Measure activity"}],
                "novelty_claim": "Tests a previously untested circuit.",
                "expected_outcomes": "Circuit activity will correlate with sleep.",
                "evidence_claims": [{
                    "claim": "Neural circuit X controls sleep.",
                    "evidence_type": "prediction",
                    "source_ids": [],
                }],
            },
            {
                "scores": {"genetic_feasibility": 4},
                "issues": ["Validate the driver line."],
                "suggestions": ["Add an independent driver."],
                "reasoning": "Promising but needs validation.",
            },
        )

        self.assertEqual(proposal["pipeline_status"], "complete")
        self.assertEqual(proposal["specific_aims"][0]["objective"], "Test activity across sleep states")
        self.assertEqual(proposal["references"][0]["pmid"], "12345678")
        self.assertEqual(proposal["references"][0]["doi"], "10.1234/example")
        self.assertEqual(proposal["references"][0]["source"], "PubMed")
        self.assertEqual(proposal["references"][0]["relevance_score"], 0.42)
        self.assertEqual(proposal["evidence_status"], "retrieved_literature")
        self.assertEqual(proposal["literature_findings"]["pathways"], ["circadian clock"])
        self.assertEqual(proposal["fly_lines_needed"][0]["line"], "UAS-GCaMP")
        self.assertEqual(len(proposal["evidence_claims"]), 2)

        planner = planner_module.ResearchPlanner(Mock())
        markdown = planner.format_proposal_for_chat(proposal)
        self.assertIn("Independent Review", markdown)
        self.assertIn("TF-IDF text-match heuristic: 0.42", markdown)
        self.assertIn("circadian clock", markdown)
        self.assertIn("UAS-GCaMP", markdown)

    def test_pipeline_proposal_marks_agent_errors_as_partial(self):
        proposal = pipeline_module.build_pipeline_proposal(
            "test topic",
            {"_error": "literature unavailable"},
            {},
            {},
        )

        self.assertEqual(proposal["pipeline_status"], "partial")
        self.assertEqual(proposal["references"], [])


class LiteratureRelevanceTests(unittest.TestCase):
    def test_gene_terms_are_separated_and_influence_literature_ranking(self):
        anthropic_stub = types.ModuleType("anthropic")
        anthropic_stub.Anthropic = object
        bio_stub = types.ModuleType("Bio")
        bio_stub.__path__ = []
        bio_stub.Entrez = types.ModuleType("Bio.Entrez")
        requests_stub = types.ModuleType("requests")
        numpy_stub = types.ModuleType("numpy")
        sklearn_stub = types.ModuleType("sklearn")
        sklearn_stub.__path__ = []
        feature_stub = types.ModuleType("sklearn.feature_extraction")
        feature_stub.__path__ = []
        text_stub = types.ModuleType("sklearn.feature_extraction.text")
        metrics_stub = types.ModuleType("sklearn.metrics")
        metrics_stub.__path__ = []
        pairwise_stub = types.ModuleType("sklearn.metrics.pairwise")

        class FakeMatrix:
            def __getitem__(self, _index):
                return self

        class FakeVectorizer:
            def __init__(self, **_kwargs):
                pass

            def fit_transform(self, corpus):
                self.corpus = corpus
                FakeVectorizer.query = corpus[0]
                FakeVectorizer.papers = corpus[1:]
                return FakeMatrix()

        def fake_cosine_similarity(_query, _papers):
            scores = [
                float("foxo" in paper.lower() and "foxo" in FakeVectorizer.query.lower())
                for paper in FakeVectorizer.papers
            ]
            return [scores]

        text_stub.TfidfVectorizer = FakeVectorizer
        pairwise_stub.cosine_similarity = fake_cosine_similarity
        stubs = {
            "anthropic": anthropic_stub,
            "Bio": bio_stub,
            "Bio.Entrez": bio_stub.Entrez,
            "requests": requests_stub,
            "numpy": numpy_stub,
            "sklearn": sklearn_stub,
            "sklearn.feature_extraction": feature_stub,
            "sklearn.feature_extraction.text": text_stub,
            "sklearn.metrics": metrics_stub,
            "sklearn.metrics.pairwise": pairwise_stub,
        }
        with patch.dict(sys.modules, stubs):
            assistant_module = importlib.import_module("drosophila_assistant")

        assistant = object.__new__(assistant_module.DrosophilaAssistant)
        papers = [
            {
                "title": "Drosophila general biology",
                "abstract": "A general study of fruit fly biology.",
                "source": "PubMed",
            },
            {
                "title": "Drosophila foxo function",
                "abstract": "A study of foxo in the fruit fly.",
                "source": "PubMed",
            },
        ]

        with patch.dict(sys.modules, stubs):
            ranked = assistant.filter_relevant_papers(
                papers, "zzzzunmatchedtopic", ["foxo"]
            )

        self.assertEqual(ranked[0]["title"], "Drosophila foxo function")
        self.assertGreater(ranked[0]["_relevance_score"], ranked[1]["_relevance_score"])
        self.assertIn("foxo foxo foxo", FakeVectorizer.query)

    def test_gene_extraction_matches_whole_tokens_only(self):
        anthropic_stub = types.ModuleType("anthropic")
        anthropic_stub.Anthropic = object
        bio_stub = types.ModuleType("Bio")
        bio_stub.__path__ = []
        bio_stub.Entrez = types.ModuleType("Bio.Entrez")
        stubs = {
            "anthropic": anthropic_stub,
            "Bio": bio_stub,
            "Bio.Entrez": bio_stub.Entrez,
            "requests": types.ModuleType("requests"),
        }
        with patch.dict(sys.modules, stubs):
            assistant_module = importlib.import_module("drosophila_assistant")
        assistant = object.__new__(assistant_module.DrosophilaAssistant)
        extract = assistant.extract_gene_names

        # Substrings of ordinary words ("Drosophila" -> dr, "intestinal" -> al) are not genes.
        self.assertEqual(
            extract("Hippo pathway dysregulation in intestinal stem cell tumorigenesis in Drosophila midgut"),
            ["hpo"],
        )
        self.assertEqual(extract("Role of yki and wts in Drosophila wing growth"), ["yki", "wts"])
        self.assertEqual(extract("How does the foxo gene affect lifespan?"), ["foxo"])
        self.assertEqual(extract("Notch signaling in neuroblast fate"), ["n"])
        self.assertEqual(extract("Does the how gene regulate muscle attachment?"), ["how"])
        self.assertEqual(extract("what is sleep"), [])

        # Chat questions become PubMed terms without filler words that PubMed would AND together.
        query = assistant.build_chat_search_query
        self.assertEqual(query("tell me about FOX genes in Drosophila"), "forkhead")
        self.assertEqual(query("tell me about fox gene in drosophilia"), "forkhead")
        self.assertEqual(query("tell me abt fox genes"), "forkhead")
        # Recognized family terms lead, so dropping trailing terms removes unknown words first.
        self.assertEqual(query("insulin regulation of FOX genes"), "forkhead insulin regulation")
        self.assertEqual(query("What is the role of Hox genes in segment identity?"), "homeotic segment identity")
        self.assertEqual(query("How does insulin affect lifespan?"), "insulin affect lifespan")


class CitationTests(unittest.TestCase):
    def test_chat_citations_only_resolve_to_retrieved_records(self):
        answer = (
            "The supplied study reports the finding [PMID:12345678]. "
            "Another claim cites a nonexistent paper [PMID:99999999].\n\n"
            "### References\n\n[1] Invented authors. Fake paper. [Source](https://fake.example)"
        )
        papers = [{
            "pmid": "12345678",
            "title": "Retrieved title",
            "authors": "Author A et al.",
            "year": "2024",
            "url": "https://pubmed.ncbi.nlm.nih.gov/12345678/",
        }]

        result = planner_module.validate_answer_citations(answer, papers)

        self.assertIn("finding [1]", result)
        self.assertNotIn("99999999", result)
        self.assertNotIn("fake.example", result)
        self.assertIn("Invalid citation marker(s) were removed", result)
        self.assertIn("[1] Author A et al. (2024). Retrieved title.", result)
        self.assertIn("https://pubmed.ncbi.nlm.nih.gov/12345678/", result)

    def test_chat_rejects_reused_numeric_citation_from_prior_turn(self):
        result = planner_module.validate_answer_citations(
            "A claim from this turn [1].",
            [{
                "pmid": "12345678",
                "title": "Retrieved title",
                "authors": "Author A et al.",
                "year": "2024",
            }],
        )

        self.assertNotIn("[1]", result)
        self.assertIn("Invalid citation marker(s) were removed", result)

    def test_chat_without_retrieved_papers_discloses_unverified_background(self):
        result = planner_module.validate_answer_citations(
            "Here is a brief general overview.", []
        )

        self.assertIn("No relevant papers were retrieved", result)
        self.assertIn("is unverified", result)

    def test_bibliography_is_numbered_from_retrieved_metadata(self):
        references = planner_module._canonical_references([
            {
                "pmid": "12345678",
                "title": "Retrieved title",
                "authors": "Author A et al.",
                "year": "2024",
            },
            {
                "title": "FlyBase record",
                "url": "https://flybase.org/reports/FBgn0000001",
            }
        ])

        self.assertEqual([ref["number"] for ref in references], [1, 2])
        self.assertIn("Retrieved title", references[0]["citation"])
        self.assertEqual(references[0]["url"], "https://pubmed.ncbi.nlm.nih.gov/12345678/")
        self.assertEqual(references[0]["pmid"], "12345678")
        self.assertEqual(references[1]["url"], "https://flybase.org/reports/FBgn0000001")

    def test_generated_bibliography_ignores_model_supplied_reference_metadata(self):
        client = Mock()
        client.messages.create.return_value = types.SimpleNamespace(
            content=[types.SimpleNamespace(
                text='{"background":"Prior work [1]. Unverified claim [99].","references":[{"number":1,"citation":"Invented source"},{"number":99,"citation":"Unretrieved source"}]}'
            )]
        )
        planner = planner_module.ResearchPlanner(client)
        proposal = planner.generate_proposal(
            topic="test topic",
            genes=["foxo"],
            papers=[{
                "pmid": "12345678",
                "title": "Retrieved title",
                "authors": "Author A et al.",
                "year": "2024",
            }]
        )

        self.assertEqual(proposal["references"], [{
            "number": 1,
            "citation": "Author A et al. (2024). Retrieved title. PMID: 12345678.",
            "url": "https://pubmed.ncbi.nlm.nih.gov/12345678/",
            "pmid": "12345678",
        }])
        self.assertEqual(proposal["evidence_status"], "retrieved_literature")
        self.assertIn(
            "Do not invent additional sources",
            client.messages.create.call_args_list[0].kwargs["system"],
        )
        self.assertEqual(proposal["background"], "Prior work [1]. Unverified claim .")

    def test_proposal_without_literature_is_marked_unverified(self):
        client = Mock()
        client.messages.create.return_value = types.SimpleNamespace(
            content=[types.SimpleNamespace(
                text='{"background":"A proposed model.","references":[]}'
            )]
        )
        planner = planner_module.ResearchPlanner(client)
        proposal = planner.generate_proposal(
            topic="test topic",
            genes=[],
            papers=[],
        )

        self.assertEqual(proposal["evidence_status"], "no_retrieved_literature")
        self.assertEqual(proposal["references"], [])
        self.assertIn(
            "No papers were retrieved",
            planner.format_proposal_for_chat(proposal),
        )


class EvidenceVerifierTests(unittest.TestCase):
    def test_markdown_answer_claims_are_checked_with_sdk_parameters_only(self):
        responses = [
            {"claims": [{
                "claim": "The study reports that foxo regulates sleep in Drosophila [1].",
                "evidence_type": "finding",
                "source_ids": ["12345678"],
            }]},
            {"reviews": [{
                "claim_id": 1,
                "status": "supported",
                "source_id": "12345678",
                "evidence_quote": "We find that foxo regulates sleep duration in Drosophila.",
                "reason": "Directly reported.",
            }]},
        ]

        # Mirrors the SDK signature: Mock() accepted an unsupported temperature argument
        # that the real client rejects, which disabled every evidence check.
        def create(*, model, max_tokens, messages, system=None, tools=None, tool_choice=None):
            self.assertEqual(tool_choice["name"], tools[0]["name"])
            return types.SimpleNamespace(
                stop_reason="tool_use",
                content=[types.SimpleNamespace(type="tool_use", input=responses.pop(0))],
            )

        client = types.SimpleNamespace(messages=types.SimpleNamespace(create=create))
        papers = [{
            "pmid": "12345678",
            "title": "Sleep regulation in flies",
            "abstract": "We find that foxo regulates sleep duration in Drosophila.",
        }]

        review = verifier_module.verify_claims(
            "The study reports that ***foxo*** regulates sleep in *Drosophila* [1].",
            papers,
            client,
        )

        self.assertEqual(review["coverage"]["assessed_claims"], 1)
        self.assertEqual(review["claims"][0]["status"], "supported")

    def test_verification_runs_in_batches_and_a_failed_batch_is_reported_incomplete(self):
        claims = [f"Fact number {n} about flies." for n in range(1, 7)]

        def create(*, model, max_tokens, messages, system=None, tools=None, tool_choice=None):
            prompt = messages[0]["content"]
            if prompt.startswith("Extract every"):
                payload = {"claims": [{"claim": c, "evidence_type": "general_background", "source_ids": []} for c in claims]}
            else:
                batch = json.loads(prompt.split("FACTUAL CLAIMS TO ASSESS:\n", 1)[1])
                if any(item["claim_id"] == 6 for item in batch):
                    raise TimeoutError("batch timed out")
                payload = {"reviews": [{"claim_id": item["claim_id"], "status": "not_found", "source_id": "",
                                        "evidence_quote": "", "reason": "x"} for item in batch]}
            return types.SimpleNamespace(stop_reason="tool_use", content=[types.SimpleNamespace(type="tool_use", input=payload)])

        client = types.SimpleNamespace(messages=types.SimpleNamespace(create=create))
        review = verifier_module.verify_claims(" ".join(claims), [{"pmid": "1", "abstract": "x"}], client)

        self.assertEqual(review["status"], "incomplete")
        self.assertEqual(review["coverage"], {"factual_claims": 6, "assessed_claims": 4, "percent": 67})

    def test_paraphrased_claims_are_assessed_and_survive_revision_recheck(self):
        answer = "| Gene | Role |\n|---|---|\n| *jumu* | Heart |\n\nThe *fkh* gene controls gut development."
        table_claim = "jumu plays a role in the heart"
        prose_claim = "The fkh gene controls gut development."
        calls = []

        def create(*, model, max_tokens, messages, system=None, tools=None, tool_choice=None):
            prompt = messages[0]["content"]
            calls.append(prompt)
            if prompt.startswith("Extract every"):
                # The recheck's extraction omits the table claim, as the model did in practice.
                claims = [{"claim": table_claim, "evidence_type": "general_background", "source_ids": []}]
                if len(calls) == 1:
                    claims.append({"claim": prose_claim, "evidence_type": "general_background", "source_ids": []})
                else:
                    claims = [{"claim": "Whether the fkh gene controls gut development is not established.",
                               "evidence_type": "general_background", "source_ids": []}]
                payload = {"claims": claims}
            elif prompt.startswith("Assess EVERY"):
                ids = [item["claim_id"] for item in json.loads(prompt.split("FACTUAL CLAIMS TO ASSESS:\n", 1)[1])]
                payload = {"reviews": [{"claim_id": i, "status": "not_found", "source_id": "", "evidence_quote": "", "reason": "x"} for i in ids]}
            else:
                ids = [item["claim_id"] for item in json.loads(prompt.split("CLAIMS:\n", 1)[1])]
                payload = {"revisions": [{"claim_id": i, "revised_claim": "Whether the fkh gene controls gut development is not established."} for i in ids]}
            return types.SimpleNamespace(stop_reason="tool_use", content=[types.SimpleNamespace(type="tool_use", input=payload)])

        client = types.SimpleNamespace(messages=types.SimpleNamespace(create=create))
        papers = [{"pmid": "1", "title": "t", "abstract": "Unrelated."}]
        revised, review = verifier_module.verify_and_revise_claims(answer, papers, client)

        self.assertIn("is not established", revised)
        self.assertIn("could not be located", review["revision"]["message"])
        statuses = {claim["claim"]: claim["status"] for claim in review["claims"]}
        self.assertEqual(statuses[table_claim], "not_found")
        self.assertEqual(review["coverage"]["assessed_claims"], review["coverage"]["factual_claims"])

    def test_supported_claim_requires_matching_verbatim_abstract_evidence(self):
        client = Mock()
        client.messages.create.side_effect = [
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "claims": [{
                    "claim": "The study reports a role in sleep regulation [1].",
                    "evidence_type": "finding",
                    "source_ids": ["12345678"],
                }]
            }))]),
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "reviews": [{
                    "claim_id": 1,
                    "status": "supported",
                    "source_id": "12345678",
                    "evidence_quote": "We find that foxo regulates sleep duration in Drosophila.",
                    "reason": "The abstract directly reports this result.",
                }]
            }))]),
        ]
        papers = [{
            "pmid": "12345678",
            "title": "Sleep regulation in flies",
            "abstract": "We find that foxo regulates sleep duration in Drosophila.",
            "url": "https://pubmed.ncbi.nlm.nih.gov/12345678/",
        }]

        _, review = verifier_module.verify_and_revise_claims(
            "The study reports a role in sleep regulation [1].",
            papers,
            client,
        )

        self.assertEqual(review["status"], "checked")
        claim = review["claims"][0]
        self.assertEqual(claim["status"], "supported")
        self.assertEqual(claim["evidence"]["source_id"], "12345678")
        self.assertEqual(
            claim["evidence"]["quote"],
            "We find that foxo regulates sleep duration in Drosophila.",
        )
        self.assertEqual(review["coverage"], {
            "factual_claims": 1,
            "assessed_claims": 1,
            "percent": 100,
        })
        self.assertEqual(client.messages.create.call_count, 2)
        self.assertEqual(review["revision"]["status"], "not_needed")
        extraction_prompt = client.messages.create.call_args_list[0].kwargs["messages"][0]["content"]
        self.assertIn('"citation_number": 1', extraction_prompt)
        self.assertIn('"source_id": "12345678"', extraction_prompt)
        formatted = verifier_module.format_evidence_review(review)
        self.assertIn("✅ **SUPPORTED** —", formatted)

    def test_quote_within_a_marked_up_abstract_sentence_is_accepted(self):
        sentence = ("In addition, knockdown of <i>TBPH</i> in <i>Drosophila</i> muscle precursor "
                    "cells reduced lifespan.")
        papers = [{"pmid": "39283487", "abstract": f"Background sentence here. {sentence} Closing sentence."}]

        evidence = verifier_module._find_verbatim_evidence(
            "Knockdown of TBPH in Drosophila muscle precursor cells reduced lifespan.", "39283487", papers
        )

        self.assertIsNotNone(evidence)
        self.assertEqual(evidence["quote"], sentence)
        # Short fragments and text absent from the abstract are still rejected.
        self.assertIsNone(verifier_module._find_verbatim_evidence("reduced lifespan.", "39283487", papers))
        self.assertIsNone(verifier_module._find_verbatim_evidence(
            "Knockdown of TBPH in Drosophila neurons extended lifespan.", "39283487", papers
        ))

    def test_claim_is_uncertain_when_reported_quote_is_not_in_abstract(self):
        client = Mock()
        client.messages.create.side_effect = [
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "claims": [{
                    "claim": "The study found lifespan extension [1].",
                    "evidence_type": "finding",
                    "source_ids": ["12345678"],
                }]
            }))]),
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "reviews": [{
                    "claim_id": 1,
                    "status": "supported",
                    "source_id": "12345678",
                    "evidence_quote": "foxo extends lifespan by fifty percent.",
                    "reason": "The abstract describes a lifespan effect.",
                }]
            }))]),
        ]

        review = verifier_module.verify_claims(
            "The study found lifespan extension [1].",
            [{
                "pmid": "12345678",
                "title": "Sleep regulation in flies",
                "abstract": "We investigate sleep behavior in Drosophila.",
            }],
            client,
        )

        self.assertEqual(review["status"], "needs_review")
        self.assertEqual(review["claims"][0]["status"], "uncertain")
        self.assertIsNone(review["claims"][0]["evidence"])

    def test_omitted_assessment_is_reported_in_claim_coverage(self):
        client = Mock()
        client.messages.create.side_effect = [
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "claims": [
                    {"claim": "Claim one [1].", "evidence_type": "finding", "source_ids": ["12345678"]},
                    {"claim": "Claim two [1].", "evidence_type": "finding", "source_ids": ["12345678"]},
                ]
            }))]),
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "reviews": [{
                    "claim_id": 1,
                    "status": "not_found",
                    "source_id": "",
                    "evidence_quote": "",
                    "reason": "No support was found.",
                }]
            }))]),
        ]

        review = verifier_module.verify_claims(
            "Claim one [1]. Claim two [1].",
            [{"pmid": "12345678", "abstract": "A sufficiently long sentence in this abstract."}],
            client,
        )

        self.assertEqual(review["status"], "incomplete")
        self.assertEqual(review["coverage"], {
            "factual_claims": 2,
            "assessed_claims": 1,
            "percent": 50,
        })
        self.assertEqual(review["claims"][1]["status"], "uncertain")
        self.assertIn("Coverage:** 1/2 factual claims assessed (50%).", verifier_module.format_evidence_review(review))

    def test_structured_claims_are_merged_and_unknown_sources_are_discarded(self):
        client = Mock()
        client.messages.create.side_effect = [
            types.SimpleNamespace(content=[types.SimpleNamespace(text='{"claims":[]}')]),
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "reviews": [{
                    "claim_id": 1,
                    "status": "supported",
                    "source_id": "12345678",
                    "evidence_quote": "We find that foxo regulates sleep duration in Drosophila.",
                    "reason": "The abstract directly reports this result.",
                }]
            }))]),
        ]
        claim = "The study reports a role in sleep regulation [1]."

        review = verifier_module.verify_claims(
            claim,
            [{
                "pmid": "12345678",
                "abstract": "We find that foxo regulates sleep duration in Drosophila.",
            }],
            client,
            structured_claims=[{
                "claim": claim,
                "evidence_type": "finding",
                "source_ids": ["12345678", "99999999"],
            }],
        )

        self.assertEqual(review["status"], "checked")
        self.assertEqual(review["claims"][0]["source_ids"], ["12345678"])
        self.assertEqual(review["coverage"]["percent"], 100)

    def test_unverified_claim_is_qualified_and_checked_again(self):
        client = Mock()
        original_claim = "The study found lifespan extension [1]."
        revised_claim = "The supplied abstracts do not establish that the study found lifespan extension [1]."
        client.messages.create.side_effect = [
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "claims": [{
                    "claim": original_claim,
                    "evidence_type": "finding",
                    "source_ids": ["12345678"],
                }]
            }))]),
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "reviews": [{
                    "claim_id": 1,
                    "status": "not_found",
                    "source_id": "",
                    "evidence_quote": "",
                    "reason": "The abstract does not report a lifespan result.",
                }]
            }))]),
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "revisions": [{"claim_id": 1, "revised_claim": revised_claim}]
            }))]),
            # Only the rewritten claim is re-verified; there is no second extraction.
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "reviews": [{
                    "claim_id": 1,
                    "status": "not_found",
                    "source_id": "",
                    "evidence_quote": "",
                    "reason": "The abstract does not establish this point.",
                }]
            }))]),
        ]
        original = {
            "answer": original_claim,
            "evidence_claims": [{
                "claim": original_claim,
                "evidence_type": "finding",
                "source_ids": ["12345678"],
            }],
            "references": [{"citation": "Author et al. (2024). Verified source."}],
        }

        revised, review = verifier_module.verify_and_revise_claims(
            original["answer"],
            [{"pmid": "12345678", "abstract": "The study investigates sleep behavior in Drosophila."}],
            client,
            structured_claims=original["evidence_claims"],
            content=original,
        )

        self.assertIn(revised_claim, revised["answer"])
        self.assertEqual(revised["evidence_claims"][0]["claim"], revised_claim)
        self.assertEqual(revised["references"], original["references"])
        self.assertEqual(review["revision"]["status"], "applied")
        self.assertEqual(review["claims"][0]["status"], "not_found")
        self.assertEqual(client.messages.create.call_count, 4)
        self.assertIn("checked again", verifier_module.format_evidence_review(review))

    def test_failed_revision_keeps_original_claim_and_reports_failure(self):
        client = Mock()
        original_claim = "The study found lifespan extension [1]."
        client.messages.create.side_effect = [
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "claims": [{
                    "claim": original_claim,
                    "evidence_type": "finding",
                    "source_ids": ["12345678"],
                }]
            }))]),
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "reviews": [{
                    "claim_id": 1,
                    "status": "not_found",
                    "source_id": "",
                    "evidence_quote": "",
                    "reason": "No matching evidence.",
                }]
            }))]),
            types.SimpleNamespace(content=[types.SimpleNamespace(text='{"revisions":[]}')]),
        ]

        revised, review = verifier_module.verify_and_revise_claims(
            original_claim,
            [{"pmid": "12345678", "abstract": "The study investigates sleep behavior in Drosophila."}],
            client,
        )

        self.assertEqual(revised, original_claim)
        self.assertEqual(review["revision"]["status"], "failed")
        self.assertEqual(review["claims"][0]["status"], "not_found")
        self.assertEqual(client.messages.create.call_count, 3)

    def test_disabled_auto_revision_leaves_reviewed_claim_unchanged(self):
        client = Mock()
        claim = "The study found lifespan extension [1]."
        client.messages.create.side_effect = [
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "claims": [{
                    "claim": claim,
                    "evidence_type": "finding",
                    "source_ids": ["12345678"],
                }]
            }))]),
            types.SimpleNamespace(content=[types.SimpleNamespace(text=json.dumps({
                "reviews": [{
                    "claim_id": 1,
                    "status": "not_found",
                    "source_id": "",
                    "evidence_quote": "",
                    "reason": "No matching abstract evidence.",
                }]
            }))]),
        ]

        revised, review = verifier_module.verify_and_revise_claims(
            claim,
            [{"pmid": "12345678", "abstract": "The study investigates sleep behavior in Drosophila."}],
            client,
            auto_revise=False,
        )

        self.assertEqual(revised, claim)
        self.assertEqual(review["revision"]["status"], "disabled")
        self.assertEqual(review["claims"][0]["status"], "not_found")
        self.assertEqual(client.messages.create.call_count, 2)

    def test_auto_revision_limit_prevents_large_revision_batches(self):
        client = Mock()
        review = {
            "status": "needs_review",
            "claims": [
                {
                    "claim_id": index,
                    "claim": f"Unverified claim {index}.",
                    "evidence_type": "finding",
                    "source_ids": ["12345678"],
                    "status": "not_found",
                }
                for index in range(1, verifier_module.MAX_AUTO_REVISED_CLAIMS + 2)
            ],
        }

        replacements, revision = verifier_module.create_claim_revisions(
            review,
            [{"pmid": "12345678", "abstract": "A retrieved abstract."}],
            client,
        )

        self.assertEqual(replacements, {})
        self.assertEqual(revision["status"], "limit_reached")
        self.assertEqual(revision["revision_limit"], 10)
        self.assertEqual(client.messages.create.call_count, 0)

    def test_checker_failure_is_returned_as_visible_unavailable_state(self):
        client = Mock()
        client.messages.create.side_effect = RuntimeError("service unavailable")

        review = verifier_module.verify_claims(
            "Claim [1].",
            [{"pmid": "12345678", "title": "Title", "abstract": "Long enough abstract text to check."}],
            client,
        )

        self.assertEqual(review["status"], "unavailable")
        self.assertIn("could not be completed", review["message"])
        self.assertIn("Treat the answer's factual claims as unverified.", verifier_module.format_evidence_review(review))


class ModelChoiceTests(unittest.TestCase):
    def test_unknown_or_missing_model_falls_back_to_default(self):
        agents_module = importlib.import_module("agents")
        self.assertEqual(agents_module.resolve_model("claude-opus-5-5"), "claude-opus-5-5")
        self.assertEqual(agents_module.resolve_model("not-a-model"), agents_module.CLAUDE_MODEL)
        self.assertEqual(agents_module.resolve_model(None), agents_module.CLAUDE_MODEL)

    def test_response_text_skips_a_leading_thinking_block(self):
        agents_module = importlib.import_module("agents")
        response = types.SimpleNamespace(content=[
            types.SimpleNamespace(type="thinking", thinking=""),
            types.SimpleNamespace(type="text", text='{"claims": []}'),
        ])
        self.assertEqual(agents_module.response_text(response), '{"claims": []}')

    def test_json_call_forces_the_tool_only_on_models_that_allow_it(self):
        seen = []

        def create(*, model, max_tokens, messages, system=None, tools=None, tool_choice=None):
            seen.append(tool_choice["type"])
            return types.SimpleNamespace(
                stop_reason="tool_use",
                content=[types.SimpleNamespace(type="tool_use", input={"claims": []})],
            )

        client = types.SimpleNamespace(messages=types.SimpleNamespace(create=create))
        for model in ("claude-sonnet-4-6", "claude-opus-5-5"):
            verifier_module._call_json_model(client, "prompt", 100, "test", "claims", model)

        self.assertEqual(seen, ["tool", "auto"])
