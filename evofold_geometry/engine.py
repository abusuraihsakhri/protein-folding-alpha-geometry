"""
Protein Structure Analysis Engine
Real implementations of Ramachandran plot analysis, secondary structure assignment,
hydrogen bond calculation, and contact map generation.
"""
import math
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field


# --- Constants ---
DEG_TO_RAD = math.pi / 180.0
RAD_TO_DEG = 180.0 / math.pi

# Ramachandran plot region definitions (phi, psi) in degrees
ALPHA_HELIX_CENTER = (-57.0, -47.0)
BETA_SHEET_CENTER = (-135.0, 135.0)
LEFT_HELIX_CENTER = (57.0, 47.0)
GLYCINE_CENTER = (0.0, 0.0)

# Tolerances for region classification
FAVORED_TOLERANCE = 30.0
ALLOWED_TOLERANCE = 50.0
GLYCINE_FAVORED_TOLERANCE = 45.0
GLYCINE_ALLOWED_TOLERANCE = 70.0

# Hydrogen bond parameters
HBOND_DISTANCE_CUTOFF = 3.5  # Angstroms
HBOND_ENERGY_COEFFICIENT = -2.7  # kcal/mol typical H-bond energy

# Van der Waals radii (Angstroms)
VDW_RADII = {
    'N': 1.55, 'CA': 1.70, 'C': 1.70, 'O': 1.52,
    'CB': 1.70, 'CG': 1.70, 'CD': 1.70, 'NE': 1.55,
    'CZ': 1.70, 'NH': 1.55, 'OH': 1.52, 'SG': 1.80,
    'SD': 1.80, 'CE': 1.70, 'NZ': 1.55,
}
VDW_DEFAULT = 1.70


@dataclass
class Residue:
    """Represents a single amino acid residue with backbone coordinates."""
    name: str
    index: int
    phi: Optional[float] = None    # degrees
    psi: Optional[float] = None    # degrees
    omega: Optional[float] = None  # degrees
    # Backbone atom coordinates (N, CA, C, O)
    n_coord: Optional[Tuple[float, float, float]] = None
    ca_coord: Optional[Tuple[float, float, float]] = None
    c_coord: Optional[Tuple[float, float, float]] = None
    o_coord: Optional[Tuple[float, float, float]] = None
    secondary_structure: Optional[str] = None  # H, E, T, C


@dataclass
class HydrogenBond:
    """Represents an N-H···O=C hydrogen bond."""
    donor_residue: int
    acceptor_residue: int
    distance: float
    energy: float
    donor_atom: str = "N"
    acceptor_atom: str = "O"


@dataclass
class RamachandranClassification:
    """Classification of a residue's phi/psi angles."""
    residue_index: int
    residue_name: str
    phi: float
    psi: float
    region: str       # 'favored', 'allowed', 'outlier'
    structure_type: str  # 'alpha_helix', 'beta_sheet', 'left_helix', 'ppII', 'other'


@dataclass
class ContactMap:
    """CA-CA contact map for a protein."""
    size: int
    contacts: List[Tuple[int, int]]
    distance_matrix: List[List[float]]
    cutoff: float = 8.0


@dataclass
class StructureAnalysisResult:
    """Complete structure analysis result."""
    num_residues: int
    secondary_structure: List[str]
    ramachandran: List[RamachandranClassification]
    hbonds: List[HydrogenBond]
    contact_map: Optional[ContactMap] = None
    summary: Dict[str, Any] = field(default_factory=dict)


# --- Ramachandran Plot Analysis ---

def _angle_distance(phi: float, psi: float, center_phi: float, center_psi: float) -> float:
    """Compute circular distance between two phi/psi pairs in degrees.
    Accounts for the periodic nature of dihedral angles."""
    dphi = phi - center_phi
    dpsi = psi - center_psi
    # Wrap to [-180, 180]
    dphi = ((dphi + 180) % 360) - 180
    dpsi = ((dpsi + 180) % 360) - 180
    return math.sqrt(dphi ** 2 + dpsi ** 2)


def classify_ramachandran(residue_name: str, phi: float, psi: float) -> RamachandranClassification:
    """Classify a residue's phi/psi angles into Ramachandran regions.
    
    Regions:
    - Favored: within tight tolerance of known secondary structure centers
    - Allowed: within wider tolerance
    - Outlier: outside allowed regions
    
    Structure types:
    - alpha_helix: near (-57, -47)
    - beta_sheet: near (-135, 135)
    - left_helix: near (57, 47)
    - ppII: polyproline II helix near (-75, 145)
    - other: none of the above
    """
    is_glycine = residue_name.upper() in ('GLY', 'G')
    
    favored_tol = GLYCINE_FAVORED_TOLERANCE if is_glycine else FAVORED_TOLERANCE
    allowed_tol = GLYCINE_ALLOWED_TOLERANCE if is_glycine else ALLOWED_TOLERANCE
    
    # Calculate distances to known centers
    dist_alpha = _angle_distance(phi, psi, *ALPHA_HELIX_CENTER)
    dist_beta = _angle_distance(phi, psi, *BETA_SHEET_CENTER)
    dist_left = _angle_distance(phi, psi, *LEFT_HELIX_CENTER)
    dist_ppII = _angle_distance(phi, psi, -75.0, 145.0)
    
    # Determine structure type (closest center)
    distances = {
        'alpha_helix': dist_alpha,
        'beta_sheet': dist_beta,
        'left_helix': dist_left,
        'ppII': dist_ppII,
    }
    min_dist_name = min(distances, key=distances.get)
    min_dist = distances[min_dist_name]
    
    # Determine region classification
    if min_dist <= favored_tol:
        region = 'favored'
    elif min_dist <= allowed_tol:
        region = 'allowed'
    else:
        region = 'outlier'
        min_dist_name = 'other'
    
    return RamachandranClassification(
        residue_index=0,  # set by caller
        residue_name=residue_name,
        phi=phi,
        psi=psi,
        region=region,
        structure_type=min_dist_name,
    )


def analyze_ramachandran(residues: List[Residue]) -> List[RamachandranClassification]:
    """Analyze all residues and classify their Ramachandran angles."""
    results = []
    for res in residues:
        if res.phi is not None and res.psi is not None:
            classification = classify_ramachandran(res.name, res.phi, res.psi)
            classification.residue_index = res.index
            results.append(classification)
    return results


def ramachandran_summary(classifications: List[RamachandranClassification]) -> Dict[str, Any]:
    """Summarize Ramachandran plot statistics."""
    total = len(classifications)
    if total == 0:
        return {'total': 0, 'favored_pct': 0, 'allowed_pct': 0, 'outlier_pct': 0}
    
    favored = sum(1 for c in classifications if c.region == 'favored')
    allowed = sum(1 for c in classifications if c.region == 'allowed')
    outlier = sum(1 for c in classifications if c.region == 'outlier')
    
    structure_counts = {}
    for c in classifications:
        structure_counts[c.structure_type] = structure_counts.get(c.structure_type, 0) + 1
    
    return {
        'total': total,
        'favored': favored,
        'allowed': allowed,
        'outlier': outlier,
        'favored_pct': 100.0 * favored / total,
        'allowed_pct': 100.0 * allowed / total,
        'outlier_pct': 100.0 * outlier / total,
        'structure_distribution': structure_counts,
    }


# --- Secondary Structure Assignment ---

def assign_secondary_structure(residues: List[Residue]) -> List[str]:
    """Assign secondary structure based on Ramachandran angles.
    
    Uses phi/psi angle criteria:
    - H (alpha helix): phi in [-90, -30], psi in [-70, -10]
    - E (beta strand): phi in [-170, -80], psi in [80, 170]
    - T (turn): specific turn criteria
    - C (coil): everything else
    """
    assignments = []
    for res in residues:
        if res.phi is None or res.psi is None:
            assignments.append('C')
            continue
        
        phi, psi = res.phi, res.psi
        
        # Alpha helix region
        if -90 <= phi <= -30 and -70 <= psi <= -10:
            assignments.append('H')
        # Beta sheet region
        elif -170 <= phi <= -80 and 80 <= psi <= 170:
            assignments.append('E')
        # Turn region (Type I and Type II turns)
        elif -90 <= phi <= -30 and -10 <= psi <= 70:
            assignments.append('T')
        # Polyproline II
        elif -90 <= phi <= -50 and 100 <= psi <= 180:
            assignments.append('P')
        else:
            assignments.append('C')
    
    # Smooth: require at least 3 consecutive residues for H or E
    smoothed = _smooth_secondary_structure(assignments)
    return smoothed


def _smooth_secondary_structure(assignments: List[str]) -> List[str]:
    """Remove isolated secondary structure assignments (require >= 3 consecutive)."""
    if len(assignments) < 3:
        return assignments[:]
    
    result = assignments[:]
    
    for ss_type in ('H', 'E'):
        i = 0
        while i < len(result):
            if result[i] == ss_type:
                # Find run length
                j = i
                while j < len(result) and result[j] == ss_type:
                    j += 1
                run_length = j - i
                if run_length < 3:
                    # Replace short run with coil
                    for k in range(i, j):
                        result[k] = 'C'
                i = j
            else:
                i += 1
    
    return result


# --- Hydrogen Bond Calculation ---

def _distance_3d(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    """Euclidean distance between two 3D points."""
    return math.sqrt(
        (a[0] - b[0]) ** 2 +
        (a[1] - b[1]) ** 2 +
        (a[2] - b[2]) ** 2
    )


def calculate_hydrogen_bonds(residues: List[Residue],
                              distance_cutoff: float = HBOND_DISTANCE_CUTOFF) -> List[HydrogenBond]:
    """Calculate backbone N-H···O=C hydrogen bonds.
    
    A hydrogen bond exists when the N-O distance is < distance_cutoff (default 3.5 Å).
    Energy is estimated using a simple distance-dependent formula:
        E = k * (1/r) where k = -2.7 kcal/mol (typical H-bond)
    
    Excludes bonds between residues closer than 2 in sequence (i, i+1, i+2).
    """
    hbonds = []
    
    for i, donor_res in enumerate(residues):
        if donor_res.n_coord is None:
            continue
        for j, acceptor_res in enumerate(residues):
            if acceptor_res.o_coord is None:
                continue
            # Exclude sequential neighbors (|i-j| < 3)
            if abs(i - j) < 3:
                continue
            
            dist = _distance_3d(donor_res.n_coord, acceptor_res.o_coord)
            
            if dist < distance_cutoff and dist > 0.5:  # 0.5 to avoid artifacts
                energy = HBOND_ENERGY_COEFFICIENT / dist
                hbonds.append(HydrogenBond(
                    donor_residue=donor_res.index,
                    acceptor_residue=acceptor_res.index,
                    distance=round(dist, 3),
                    energy=round(energy, 3),
                ))
    
    return hbonds


# --- Contact Map Generation ---

def generate_contact_map(residues: List[Residue],
                          cutoff: float = 8.0) -> ContactMap:
    """Generate a CA-CA contact map.
    
    Two residues are in contact if their CA atoms are within cutoff distance (default 8.0 Å).
    Returns the full distance matrix and list of contacts.
    """
    n = len(residues)
    distance_matrix = [[0.0] * n for _ in range(n)]
    contacts = []
    
    for i in range(n):
        for j in range(i + 1, n):
            if residues[i].ca_coord is None or residues[j].ca_coord is None:
                continue
            dist = _distance_3d(residues[i].ca_coord, residues[j].ca_coord)
            distance_matrix[i][j] = round(dist, 3)
            distance_matrix[j][i] = round(dist, 3)
            if dist <= cutoff:
                contacts.append((residues[i].index, residues[j].index))
    
    return ContactMap(
        size=n,
        contacts=contacts,
        distance_matrix=distance_matrix,
        cutoff=cutoff,
    )


# --- Dihedral Angle Calculation from Coordinates ---

def _dihedral_angle(p1: Tuple[float, float, float],
                     p2: Tuple[float, float, float],
                     p3: Tuple[float, float, float],
                     p4: Tuple[float, float, float]) -> float:
    """Calculate dihedral angle (in degrees) from four 3D points.
    
    Uses the atan2 method for numerical stability.
    """
    b1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    b2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])
    b3 = (p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2])
    
    # Cross products
    n1 = _cross(b1, b2)
    n2 = _cross(b2, b3)
    
    # Normalize b2
    b2_mag = math.sqrt(b2[0]**2 + b2[1]**2 + b2[2]**2)
    if b2_mag < 1e-10:
        return 0.0
    b2_hat = (b2[0]/b2_mag, b2[1]/b2_mag, b2[2]/b2_mag)
    
    # m1 = n1 x b2_hat
    m1 = _cross(n1, b2_hat)
    
    x = _dot(n1, n2)
    y = _dot(m1, n2)
    
    angle = math.atan2(y, x) * RAD_TO_DEG
    return angle


def _cross(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """Cross product of two 3D vectors."""
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    """Dot product of two 3D vectors."""
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def compute_dihedral_angles(residues: List[Residue]) -> List[Residue]:
    """Compute phi, psi, omega dihedral angles from backbone coordinates.
    
    Phi: C(i-1) - N(i) - CA(i) - C(i)
    Psi: N(i) - CA(i) - C(i) - N(i+1)
    Omega: CA(i-1) - C(i-1) - N(i) - CA(i)
    """
    for i, res in enumerate(residues):
        if res.n_coord is None or res.ca_coord is None or res.c_coord is None:
            continue
        
        # Phi: requires C of previous residue
        if i > 0 and residues[i-1].c_coord is not None:
            res.phi = _dihedral_angle(
                residues[i-1].c_coord, res.n_coord, res.ca_coord, res.c_coord
            )
        
        # Psi: requires N of next residue
        if i < len(residues) - 1 and residues[i+1].n_coord is not None:
            res.psi = _dihedral_angle(
                res.n_coord, res.ca_coord, res.c_coord, residues[i+1].n_coord
            )
        
        # Omega: requires CA and C of previous residue
        if i > 0 and residues[i-1].ca_coord is not None and residues[i-1].c_coord is not None:
            res.omega = _dihedral_angle(
                residues[i-1].ca_coord, residues[i-1].c_coord, res.n_coord, res.ca_coord
            )
    
    return residues


# --- Full Analysis Pipeline ---

def analyze_structure(residues: List[Residue],
                       hbond_cutoff: float = HBOND_DISTANCE_CUTOFF,
                       contact_cutoff: float = 8.0) -> StructureAnalysisResult:
    """Run complete protein structure analysis.
    
    1. Compute dihedral angles from coordinates (if angles not already set)
    2. Assign secondary structure
    3. Classify Ramachandran angles
    4. Calculate hydrogen bonds
    5. Generate contact map
    """
    # Compute angles from coordinates if not set
    has_coords = any(r.n_coord is not None for r in residues)
    has_angles = any(r.phi is not None for r in residues)
    
    if has_coords and not has_angles:
        residues = compute_dihedral_angles(residues)
    
    # Assign secondary structure
    ss_assignments = assign_secondary_structure(residues)
    for i, ss in enumerate(ss_assignments):
        if i < len(residues):
            residues[i].secondary_structure = ss
    
    # Ramachandran analysis
    rama = analyze_ramachandran(residues)
    
    # Hydrogen bonds
    hbonds = calculate_hydrogen_bonds(residues, hbond_cutoff)
    
    # Contact map
    contact_map = None
    if has_coords:
        contact_map = generate_contact_map(residues, contact_cutoff)
    
    # Summary
    rama_summary = ramachandran_summary(rama)
    ss_counts = {}
    for ss in ss_assignments:
        ss_counts[ss] = ss_counts.get(ss, 0) + 1
    
    summary = {
        'ramachandran': rama_summary,
        'secondary_structure_counts': ss_counts,
        'num_hbonds': len(hbonds),
        'num_contacts': len(contact_map.contacts) if contact_map else 0,
    }
    
    return StructureAnalysisResult(
        num_residues=len(residues),
        secondary_structure=ss_assignments,
        ramachandran=rama,
        hbonds=hbonds,
        contact_map=contact_map,
        summary=summary,
    )
