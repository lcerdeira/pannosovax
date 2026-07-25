"""Testes mínimos da API do PanNosoVax Studio.

Cobrem o contrato dos endpoints sem depender de o pipeline rodar de verdade:
estrutura das etapas, estado inicial e validação de entrada.
"""
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from app.backend import app, STAGES  # noqa: E402

client = TestClient(app)


def test_stages_config_bem_formada():
    ids = [s[0] for s in STAGES]
    assert ids, "deve haver etapas configuradas"
    assert len(ids) == len(set(ids)), "ids de etapa não podem repetir"
    for sid, name, help_ in STAGES:
        assert sid and name and help_, "toda etapa tem id, nome e ajuda"


def test_get_stages_retorna_lista():
    r = client.get("/api/stages")
    assert r.status_code == 200
    body = r.json()
    assert "stages" in body
    ids = {s["id"] for s in body["stages"]}
    assert ids == {s[0] for s in STAGES}
    for s in body["stages"]:
        assert set(s) >= {"id", "name", "help", "pending"}
        assert isinstance(s["pending"], bool)


def test_status_inicial_ocioso():
    r = client.get("/api/status")
    assert r.status_code == 200
    assert r.json()["state"] in {"idle", "running", "done", "error", "cancelled"}


def test_log_retorna_linhas():
    r = client.get("/api/log")
    assert r.status_code == 200
    assert isinstance(r.json()["lines"], list)


def test_run_rejeita_etapa_invalida():
    r = client.post("/api/run", json={"stages": ["etapa_que_nao_existe"]})
    assert r.status_code == 200
    assert r.json()["ok"] is False


def test_cancel_sem_job():
    r = client.post("/api/cancel")
    assert r.json()["ok"] is False
