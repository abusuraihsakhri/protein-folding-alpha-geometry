"""REST API and optional local web frontend."""
import hmac
import os
import time
from pathlib import Path

from fastapi import FastAPI, Header, HTTPException
from fastapi.responses import PlainTextResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from .base import AuditLogger, SecurityException
from .metrics import GLOBAL_METRICS
from .models import SystemTaskPayload
from .supervisor import SystemSupervisor

supervisor = SystemSupervisor(model_provider="mock")

app = FastAPI(
    title="Protein Folding Alpha Geometry API",
    description="Heuristic protein-geometry analysis and demonstration task-audit API",
    version="3.0.1",
)


class ChatRequest(BaseModel):
    query: str


@app.get("/health")
def health():
    return {"status": "HEALTHY", "service": "protein-folding-alpha-geometry", "version": "3.0.1"}


@app.get("/metrics", response_class=PlainTextResponse)
def metrics():
    return GLOBAL_METRICS.export_prometheus_text()


@app.post("/api/audit")
def api_audit(payload: SystemTaskPayload):
    started = time.monotonic()
    try:
        dossier = supervisor.process_task(payload)
    except SecurityException:
        GLOBAL_METRICS.record_phi_block()
        raise HTTPException(status_code=422, detail="Sensitive identifier detected in input")
    GLOBAL_METRICS.record_task(dossier.overall_urgency.value, time.monotonic() - started)
    return dossier.to_dict()


@app.post("/api/chat")
def api_chat(req: ChatRequest):
    try:
        return {"response": supervisor.query_supervisory_chat(req.query)}
    except SecurityException:
        GLOBAL_METRICS.record_phi_block()
        raise HTTPException(status_code=422, detail="Sensitive identifier detected in input")


@app.get("/api/audit/logs")
def api_audit_logs(authorization: str | None = Header(default=None)):
    # Audit metadata must not be exposed anonymously. Disabled unless a
    # separate bearer token has been configured server-side.
    token = os.getenv("AUDIT_LOG_ACCESS_TOKEN")
    expected = f"Bearer {token}" if token else None
    if not expected or not hmac.compare_digest(authorization or "", expected):
        raise HTTPException(status_code=403, detail="Audit log access denied")
    return {"audit_trail": AuditLogger.get_trail(), "verified": AuditLogger.verify_integrity()}


# Keep this catch-all mount below API routes.
app.mount("/evofold_geometry", StaticFiles(directory=Path(__file__).resolve().parent.parent / "evofold_geometry"), name="engine")
app.mount("/", StaticFiles(directory=Path(__file__).resolve().parent.parent / "web", html=True), name="web")
