"""
Data models for Protein Structure Analysis.
Re-exports from engine for convenience.
"""
from .engine import (
    Residue,
    HydrogenBond,
    RamachandranClassification,
    ContactMap,
    StructureAnalysisResult,
)

__all__ = [
    'Residue',
    'HydrogenBond',
    'RamachandranClassification',
    'ContactMap',
    'StructureAnalysisResult',
]
