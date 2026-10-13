"""Regression coverage for geometry, audit integrity, API and browser assets."""
import json
import math
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError

from agents.api import app
from agents.base import AuditTrail, SecurityException
from agents.models import SystemTaskPayload
from agents.supervisor import SystemSupervisor
from evofold_geometry.cli import _parse_residues_from_json
from evofold_geometry.engine import (
    Residue, analyze_structure, calculate_hydrogen_bonds, generate_contact_map
)
from evofold_geometry.server import create_app


def test_audit_hmac_detects_modified_record():
    ledger = AuditTrail(secret_key="test-key-for-unit-testing")
    ledger.log("worker", "test", "EVENT", {"value": 1})
    assert ledger.verify_integrity()
    copy = ledger.get_trail()
    copy[0]["actor"] = "forged"
    assert ledger.verify_integrity()  # A consumer cannot edit the ledger.
    ledger.logs[0]["event_type"] = "FORGED"
    assert not ledger.verify_integrity()


def test_audit_hmac_detects_replaced_signature():
    ledger = AuditTrail(secret_key="test-key-for-unit-testing")
    ledger.log("worker", "test", "EVENT", {"value": 1})
    ledger.logs[0]["current_hash"] = "0" * 64
    assert not ledger.verify_integrity()


def test_supervisor_rejects_identifiers_in_attributes():
    sup = SystemSupervisor(model_provider="mock")
    payload = SystemTaskPayload(
        task_id="T1", target_identifier="SPECIMEN-01", primary_metric=1,
        attributes={"note": "Patient MRN-994827"},
    )
    with pytest.raises(SecurityException):
        sup.process_task(payload)


def test_numerical_schema_rejects_nonfinite():
    for invalid in (math.inf, -math.inf, math.nan):
        with pytest.raises(ValidationError):
            SystemTaskPayload(task_id="T1", target_identifier="KEY", primary_metric=invalid)


def test_partial_torsions_preserved_and_completed():
    residues = [
        Residue("ALA", 1, phi=-57, psi=-47,
                n_coord=(0, 0, 0), ca_coord=(1, 0, 0), c_coord=(1, 1, 0)),
        Residue("ALA", 2, psi=-42,
                n_coord=(2, 1, 1), ca_coord=(2, 2, 1), c_coord=(3, 2, 1)),
        Residue("ALA", 3,
                n_coord=(4, 3, 0), ca_coord=(5, 3, 0), c_coord=(5, 4, 0)),
    ]
    result = analyze_structure(residues)
    assert result.num_residues == 3
    assert residues[0].phi == -57
    assert residues[0].psi == -47
    assert residues[1].psi == -42
    assert residues[1].phi is not None
    assert residues[2].phi is not None


@pytest.mark.parametrize("cutoff", [-1, 0, math.inf, math.nan])
def test_invalid_cutoffs_rejected(cutoff):
    with pytest.raises(ValueError):
        calculate_hydrogen_bonds([], cutoff)
    with pytest.raises(ValueError):
        generate_contact_map([], cutoff)


@pytest.mark.parametrize("data", [
    [{"name":"ALA","index":1,"phi":float("nan"),"psi":-47}],
    [{"name":"ALA","index":1,"ca":[1,2]}],
    [{"name":"ALA","index":1,"n":["oops",0,1]}],
    [{"name":"ALA","index":"first","phi":-57,"psi":-47}],
])
def test_cli_rejects_invalid_residue_data(tmp_path, data):
    input_file = tmp_path / "invalid.json"
    input_file.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        _parse_residues_from_json(str(input_file))


def test_fastapi_launch_and_static_assets(monkeypatch):
    assert create_app() is app
    client = TestClient(app)
    assert client.get("/health").status_code == 200
    assert "Protein Geometry Analyzer" in client.get("/").text
    assert "def analyze_structure(" in client.get("/evofold_geometry/engine.py").text
    assert client.get("/metrics").headers["content-type"].startswith("text/plain")
    monkeypatch.delenv("AUDIT_LOG_ACCESS_TOKEN", raising=False)
    assert client.get("/api/audit/logs").status_code == 403
    monkeypatch.setenv("AUDIT_LOG_ACCESS_TOKEN", "unit-test-private-token")
    assert client.get("/api/audit/logs").status_code == 403
    assert client.get("/api/audit/logs",headers={"Authorization":"Bearer unit-test-private-token"}).status_code == 200


def test_fastapi_audit_and_guard():
    client = TestClient(app)
    good = {"task_id":"T1","target_identifier":"KEY-1","primary_metric":10}
    result = client.post("/api/audit", json=good)
    assert result.status_code == 200
    assert result.json()["task_id"] == "T1"
    bad = {**good,"attributes":{"note":"MRN-123456"}}
    assert client.post("/api/audit", json=bad).status_code == 422


def test_mock_provider_does_not_claim_validation():
    from agents.llm_factory import LLMFactory
    output = LLMFactory.create("mock").invoke("Describe specifications")
    assert "No model inference" in output
    with pytest.raises(ValueError):
        LLMFactory.create("openai")


def test_static_builder():
    from scripts.build_pages import build, DIST
    build()
    for name in ("index.html","app.js","operations.html",
                 "evofold_geometry/cli.py","evofold_geometry/engine.py"):
        assert (DIST/name).is_file()


def test_short_secondary_structure_runs_even_for_short_sequences():
    from evofold_geometry.engine import _smooth_secondary_structure
    assert _smooth_secondary_structure(['H']) == ['C']
    assert _smooth_secondary_structure(['E', 'E']) == ['C', 'C']
    assert _smooth_secondary_structure(['H', 'H', 'H']) == ['H', 'H', 'H']
