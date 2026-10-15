import importlib
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


class CitationTests(unittest.TestCase):
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
        self.assertEqual(proposal["background"], "Prior work [1]. Unverified claim .")
