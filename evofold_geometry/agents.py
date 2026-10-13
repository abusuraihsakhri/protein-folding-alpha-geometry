"""Compatibility wrappers for the task-audit demonstration.

These are NOT protein structure prediction or SE(3)-equivariant agents.
The maintained implementation is in the top-level agents package.
Install the optional API dependencies before importing this module.
"""
from agents.models import SystemTaskPayload
from agents.supervisor import SystemSupervisor
from agents.workers import (
    InvariantQCWorker,
    ProtocolConformanceWorker,
    SafetyEscalationWorker,
)


class SE3EquivariantFrameAgent:
    """Legacy name for the demonstration's primary-metric worker."""

    @staticmethod
    def audit(payload: SystemTaskPayload):
        return InvariantQCWorker.evaluate(payload)


class DihedralRamachandranAuditorAgent:
    """Legacy name for the demonstration's critical-flag worker."""

    @staticmethod
    def audit(payload: SystemTaskPayload):
        return SafetyEscalationWorker.evaluate(payload)


class ClashPenaltyCalculatorAgent:
    """Legacy name for the demonstration's descriptor worker."""

    @staticmethod
    def audit(payload: SystemTaskPayload):
        return ProtocolConformanceWorker.evaluate(payload)


class EvoFoldCoordinator:
    """Adapter preserving the old process / query method names."""

    def __init__(self):
        self._supervisor = SystemSupervisor(model_provider="mock")
        self.sub_1 = SE3EquivariantFrameAgent()
        self.sub_2 = DihedralRamachandranAuditorAgent()
        self.sub_3 = ClashPenaltyCalculatorAgent()
        self.execution_ledger = {}

    def process(self, payload: SystemTaskPayload):
        if isinstance(payload, dict):
            payload = SystemTaskPayload.model_validate(payload)
        dossier = self._supervisor.process_task(payload)
        result = dossier.to_dict()
        self.execution_ledger[payload.task_id] = result
        return result

    def query_supervisory_chat(self, query: str):
        return self._supervisor.query_supervisory_chat(query)
