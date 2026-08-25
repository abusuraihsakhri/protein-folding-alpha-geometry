"""
Protein Structure Analysis: Ramachandran plots, secondary structure assignment,
hydrogen bond calculation, and contact map generation.
"""
__version__ = "3.0.0"

from .engine import (
    Residue,
    HydrogenBond,
    RamachandranClassification,
    ContactMap,
    StructureAnalysisResult,
    classify_ramachandran,
    analyze_ramachandran,
    ramachandran_summary,
    assign_secondary_structure,
    calculate_hydrogen_bonds,
    generate_contact_map,
    compute_dihedral_angles,
    analyze_structure,
)
