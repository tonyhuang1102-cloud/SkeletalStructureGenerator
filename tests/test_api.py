"""Exercise the app through its actual HTTP boundary."""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture
def client():
    with TestClient(app) as session:
        yield session


def test_frontend_and_health(client):
    assert client.get("/").status_code == 200
    assert client.get("/static/styles.css").status_code == 200
    assert client.get("/static/app.js").status_code == 200
    assert client.get("/api/health").json()["status"] == "ok"


def test_real_analysis(client):
    response = client.post("/api/analyze", json={"input": "CC(=O)O"})
    assert response.status_code == 200
    data = response.json()
    assert data["molecular_formula"] == "C2H4O2"
    assert data["atom_count"] == 4
    assert data["bond_count"] == 3
    assert "<svg" in data["svg"]
    assert any(group["name"].lower() == "carboxylic acid" for group in data["functional_groups"])
    assert {atom["hybridization"] for atom in data["atoms"]} == {"sp2", "sp3"}


def test_condensed_request(client):
    response = client.post("/api/analyze", json={"input": "CH3CH2OH", "format": "condensed"})
    assert response.status_code == 200
    assert response.json()["canonical_smiles"] == "CCO"


@pytest.mark.parametrize("payload", [
    {"input": ""},
    {"input": "C" * 501},
    {"input": "CCO", "format": "unknown"},
    {"input": "C2H6O", "format": "condensed"},
    {"input": "not a molecule"},
    {"input": "C", "unexpected": True},
    {"input": "C", "highlight_groups": ["alcohol"] * 31},
])
def test_bad_requests_are_rejected(client, payload):
    assert client.post("/api/analyze", json=payload).status_code == 422


def test_highlighting_and_labels(client):
    molecule = client.post("/api/analyze", json={"input": "CCO"}).json()
    group_id = molecule["functional_groups"][0]["id"]
    response = client.post("/api/analyze", json={
        "input": "CCO", "show_hybridization": True, "highlight_groups": [group_id]
    })
    assert response.status_code == 200
    assert response.json()["svg"] != molecule["svg"]


def test_user_text_cannot_be_served_as_svg(client):
    response = client.post("/api/analyze", json={"input": "<script>alert(1)</script>"})
    assert response.status_code == 422
