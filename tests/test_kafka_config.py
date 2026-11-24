import ast
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import uuid

from ruamel.yaml import YAML, YAMLError


STARTUP = Path(__file__).resolve().parents[1] / "startup" / "00-startup.py"


def kafka_startup_code():
    module = ast.parse(STARTUP.read_text())
    start = next(
        index for index, node in enumerate(module.body)
        if isinstance(node, ast.FunctionDef) and node.name == "load_kafka_config"
    )
    end = next(
        index for index in range(start, len(module.body))
        if isinstance(module.body[index], ast.Assign)
        and any(
            isinstance(target, ast.Name) and target.id == "kafka_publisher"
            for target in module.body[index].targets
        )
    )
    return compile(ast.Module(body=module.body[start:end + 1], type_ignores=[]), str(STARTUP), "exec")


class KafkaConfigTests(unittest.TestCase):
    def run_startup(self, config_path, environment=None):
        publisher = Mock()
        namespace = {
            "YAML": YAML,
            "YAMLError": YAMLError,
            "Path": lambda path: config_path,
            "certifi": SimpleNamespace(where=lambda: "/default/ca.pem"),
            "Publisher": publisher,
            "uuid": uuid,
            "os": os,
        }
        with patch.dict(os.environ, environment or {}, clear=True):
            exec(kafka_startup_code(), namespace)
        return namespace, publisher

    def test_yaml_settings_are_used_by_publisher(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "kafka.yml"
            path.write_text(
                "bootstrap_servers: ['broker1:9092', 'broker2:9092']\n"
                "runengine_producer_config:\n"
                "  sasl.password: yaml-password\n"
                "  sasl.username: yaml-user\n"
                "  acks: all\n"
            )
            namespace, publisher = self.run_startup(path, {
                "BLUESKY_KAFKA_BOOTSTRAP_SERVERS": "ignored:9092",
                "BLUESKY_KAFKA_PASSWORD": "ignored-password",
            })

        kwargs = publisher.call_args.kwargs
        self.assertEqual(kwargs["bootstrap_servers"], "broker1:9092,broker2:9092")
        self.assertEqual(kwargs["producer_config"]["sasl.password"], "yaml-password")
        self.assertEqual(kwargs["producer_config"]["sasl.username"], "yaml-user")
        self.assertEqual(kwargs["producer_config"]["acks"], "all")
        self.assertEqual(kwargs["producer_config"]["message.timeout.ms"], 3000)
        self.assertEqual(kwargs["producer_config"]["ssl.ca.location"], "/default/ca.pem")
        self.assertEqual(kwargs["topic"], "hxn.bluesky.datum.documents")
        self.assertTrue(kwargs["flush_on_stop_doc"])
        self.assertIs(namespace["kafka_publisher"], publisher.return_value)

    def test_invalid_configuration_fails_without_exposing_values(self):
        configs = [
            ("", "must be a mapping"),
            ("bootstrap_servers: []\nrunengine_producer_config: {sasl.password: dummy}\n", "bootstrap_servers"),
            ("bootstrap_servers: ['broker:9092,other:9092']\nrunengine_producer_config: {sasl.password: dummy}\n", "bootstrap_servers"),
            ("bootstrap_servers: ['broker:9092']\nrunengine_producer_config: {}\n", "sasl.password"),
            ("bootstrap_servers: ['broker:9092']\nrunengine_producer_config: []\n", "runengine_producer_config"),
            ("bootstrap_servers: ['broker:9092']\nrunengine_producer_config: {sasl.password: secret-test-value\n", "Invalid YAML"),
        ]
        with TemporaryDirectory() as directory:
            path = Path(directory) / "kafka.yml"
            for config, message in configs:
                with self.subTest(message=message):
                    path.write_text(config)
                    with self.assertRaisesRegex(RuntimeError, message) as error:
                        self.run_startup(path)
                    self.assertNotIn("secret-test-value", str(error.exception))
            path.unlink()
            with self.assertRaisesRegex(RuntimeError, "Cannot read Kafka configuration"):
                self.run_startup(path, {
                    "BLUESKY_KAFKA_BOOTSTRAP_SERVERS": "ignored:9092",
                    "BLUESKY_KAFKA_PASSWORD": "ignored-password",
                })

    def test_azure_testing_still_requires_configuration(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "kafka.yml"
            with self.assertRaisesRegex(RuntimeError, "Cannot read Kafka configuration"):
                self.run_startup(path, {"AZURE_TESTING": "1"})
            path.write_text(
                "bootstrap_servers: ['broker:9092']\n"
                "runengine_producer_config: {sasl.password: dummy}\n"
            )
            namespace, publisher = self.run_startup(path, {"AZURE_TESTING": "1"})
        publisher.assert_not_called()
        self.assertIsNone(namespace["kafka_publisher"])


if __name__ == "__main__":
    unittest.main()
