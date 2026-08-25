"""
Real tests for Protein Structure Analysis.
Tests Ramachandran classification, secondary structure assignment,
hydrogen bond calculation, contact maps, and dihedral angle computation.
"""
import sys
import math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from evofold_geometry.engine import (
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
    _angle_distance,
    _distance_3d,
    _dihedral_angle,
    _cross,
    _dot,
    _smooth_secondary_structure,
    ALPHA_HELIX_CENTER,
    BETA_SHEET_CENTER,
    LEFT_HELIX_CENTER,
)


# --- Helper to create test residues ---

def _make_residue(name, index, phi=None, psi=None, n=None, ca=None, c=None, o=None):
    return Residue(
        name=name, index=index, phi=phi, psi=psi,
        n_coord=n, ca_coord=ca, c_coord=c, o_coord=o,
    )


def _make_alpha_helix_residues(count=10):
    """Create residues with typical alpha helix phi/psi angles."""
    residues = []
    for i in range(count):
        # Slight variation around ideal helix angles
        phi = -57.0 + (i % 3 - 1) * 5.0
        psi = -47.0 + (i % 3 - 1) * 5.0
        residues.append(_make_residue('ALA', i, phi=phi, psi=psi))
    return residues


def _make_beta_sheet_residues(count=10):
    """Create residues with typical beta sheet phi/psi angles."""
    residues = []
    for i in range(count):
        phi = -135.0 + (i % 3 - 1) * 5.0
        psi = 135.0 + (i % 3 - 1) * 5.0
        residues.append(_make_residue('VAL', i, phi=phi, psi=psi))
    return residues


def _make_mixed_residues():
    """Create a mix of helix and sheet residues."""
    residues = []
    # Helix residues (0-4)
    for i in range(5):
        residues.append(_make_residue('ALA', i, phi=-57.0, psi=-47.0))
    # Sheet residues (5-9)
    for i in range(5, 10):
        residues.append(_make_residue('VAL', i, phi=-135.0, psi=135.0))
    return residues


def _make_residues_with_coords():
    """Create residues with 3D coordinates for H-bond and contact map testing."""
    residues = []
    n_res = 8
    for i in range(n_res):
        # Simple linear chain along x-axis with slight z-offset for backbone
        x = i * 3.8  # ~3.8 A per residue (CA-CA distance)
        n_coord = (x, 0.0, 1.0)
        ca_coord = (x + 0.5, 0.0, 0.0)
        c_coord = (x + 1.5, 0.0, 0.5)
        o_coord = (x + 1.5, 1.2, 0.5)
        residues.append(_make_residue('ALA', i, n=n_coord, ca=ca_coord, c=c_coord, o=o_coord))
    return residues


# --- Test Ramachandran Classification ---

def test_classify_alpha_helix():
    """Residue at ideal alpha helix center should be classified as favored alpha_helix."""
    result = classify_ramachandran('ALA', -57.0, -47.0)
    assert result.region == 'favored'
    assert result.structure_type == 'alpha_helix'
    assert result.phi == -57.0
    assert result.psi == -47.0


def test_classify_beta_sheet():
    """Residue at ideal beta sheet center should be classified as favored beta_sheet."""
    result = classify_ramachandran('ALA', -135.0, 135.0)
    assert result.region == 'favored'
    assert result.structure_type == 'beta_sheet'


def test_classify_left_helix():
    """Residue at left-handed helix center should be classified as favored left_helix."""
    result = classify_ramachandran('ALA', 57.0, 47.0)
    assert result.region == 'favored'
    assert result.structure_type == 'left_helix'


def test_classify_outlier():
    """Residue far from any known center should be classified as outlier."""
    result = classify_ramachandran('ALA', 100.0, -100.0)
    assert result.region == 'outlier'


def test_classify_glycine_wider_regions():
    """Glycine should have wider allowed regions than other residues."""
    # (10, 10) is ~59.8 deg from left helix center (57, 47)
    # For glycine (allowed tolerance 70): within allowed
    # For ALA (allowed tolerance 50): outside allowed -> outlier
    gly_result = classify_ramachandran('GLY', 10.0, 10.0)
    ala_result = classify_ramachandran('ALA', 10.0, 10.0)
    # Glycine should be more permissive
    assert gly_result.region in ('favored', 'allowed')
    # ALA at same angle should be outlier
    assert ala_result.region == 'outlier'


def test_classify_glycine_alias():
    """Both 'GLY' and 'G' should be recognized as glycine."""
    r1 = classify_ramachandran('GLY', 20.0, 20.0)
    r2 = classify_ramachandran('G', 20.0, 20.0)
    assert r1.region == r2.region


def test_analyze_ramachandran_multiple():
    """analyze_ramachandran should classify all residues with angles."""
    residues = _make_alpha_helix_residues(5)
    results = analyze_ramachandran(residues)
    assert len(results) == 5
    for r in results:
        assert r.region == 'favored'
        assert r.structure_type == 'alpha_helix'


def test_analyze_ramachandran_skips_no_angles():
    """Residues without angles should be skipped."""
    residues = [
        _make_residue('ALA', 0, phi=-57.0, psi=-47.0),
        _make_residue('ALA', 1),  # no angles
        _make_residue('ALA', 2, phi=-135.0, psi=135.0),
    ]
    results = analyze_ramachandran(residues)
    assert len(results) == 2


def test_ramachandran_summary_statistics():
    """Summary should report correct percentages."""
    residues = _make_alpha_helix_residues(8) + _make_beta_sheet_residues(2)
    classifications = analyze_ramachandran(residues)
    summary = ramachandran_summary(classifications)
    assert summary['total'] == 10
    assert summary['favored'] == 10  # all should be favored
    assert abs(summary['favored_pct'] - 100.0) < 0.1


def test_ramachandran_summary_empty():
    """Empty input should return zero summary."""
    summary = ramachandran_summary([])
    assert summary['total'] == 0
    assert summary['favored_pct'] == 0


# --- Test Secondary Structure Assignment ---

def test_assign_ss_alpha_helix():
    """Residues in alpha helix region should be assigned 'H'."""
    residues = _make_alpha_helix_residues(5)
    assignments = assign_secondary_structure(residues)
    # With 5 consecutive residues in helix region, all should be H
    assert all(ss == 'H' for ss in assignments)


def test_assign_ss_beta_sheet():
    """Residues in beta sheet region should be assigned 'E'."""
    residues = _make_beta_sheet_residues(5)
    assignments = assign_secondary_structure(residues)
    assert all(ss == 'E' for ss in assignments)


def test_assign_ss_coil():
    """Residues with unusual angles should be assigned 'C'."""
    residues = [
        _make_residue('ALA', i, phi=100.0, psi=-100.0)
        for i in range(5)
    ]
    assignments = assign_secondary_structure(residues)
    assert all(ss == 'C' for ss in assignments)


def test_assign_ss_short_runs_removed():
    """Isolated H or E assignments (< 3 consecutive) should be smoothed to C."""
    # Mix: H H E E E H H -> the 2 Hs at start and end should become C
    residues = []
    # 2 helix (too short)
    for i in range(2):
        residues.append(_make_residue('ALA', i, phi=-57.0, psi=-47.0))
    # 3 sheet (long enough)
    for i in range(2, 5):
        residues.append(_make_residue('VAL', i, phi=-135.0, psi=135.0))
    # 2 helix (too short)
    for i in range(5, 7):
        residues.append(_make_residue('ALA', i, phi=-57.0, psi=-47.0))
    
    assignments = assign_secondary_structure(residues)
    assert assignments[0] == 'C'  # short helix run smoothed
    assert assignments[1] == 'C'
    assert assignments[2] == 'E'
    assert assignments[3] == 'E'
    assert assignments[4] == 'E'
    assert assignments[5] == 'C'  # short helix run smoothed
    assert assignments[6] == 'C'


def test_assign_ss_no_angles():
    """Residues without angles should be assigned 'C'."""
    residues = [_make_residue('ALA', i) for i in range(3)]
    assignments = assign_secondary_structure(residues)
    assert all(ss == 'C' for ss in assignments)


# --- Test Hydrogen Bond Calculation ---

def test_hbonds_basic():
    """Test basic hydrogen bond detection with close N-O pairs."""
    # Place residues in a helical arrangement where N(i) is close to O(i-3)
    # Typical alpha helix: N of residue i is ~2.8 A from O of residue i-3
    residues = []
    for i in range(8):
        # Helical arrangement: each residue offset by ~1.5 A in z
        z = i * 1.5
        n_coord = (0.0, 0.0, z)
        ca_coord = (1.0, 0.5, z)
        c_coord = (2.0, 0.0, z)
        o_coord = (2.0, 1.2, z)
        residues.append(_make_residue('ALA', i, n=n_coord, ca=ca_coord, c=c_coord, o=o_coord))
    
    hbonds = calculate_hydrogen_bonds(residues, distance_cutoff=3.5)
    # With this arrangement, N(i) at z=i*1.5 and O(i-3) at z=(i-3)*1.5
    # Distance = sqrt(4 + 1.44 + (3*1.5)^2) = sqrt(4+1.44+20.25) = sqrt(25.69) ≈ 5.07
    # That's too far. Let's use a tighter arrangement.
    residues2 = []
    for i in range(8):
        z = i * 0.8
        n_coord = (0.0, 0.0, z)
        ca_coord = (1.0, 0.0, z)
        c_coord = (1.5, 0.0, z)
        o_coord = (1.5, 1.0, z)
        residues2.append(_make_residue('ALA', i, n=n_coord, ca=ca_coord, c=c_coord, o=o_coord))
    
    hbonds = calculate_hydrogen_bonds(residues2, distance_cutoff=3.5)
    # N(i) at (0,0,i*0.8), O(j) at (1.5,1.0,j*0.8)
    # For i=3, j=0: dist = sqrt(2.25+1+5.76) = sqrt(9.01) ≈ 3.0 < 3.5 ✓
    assert len(hbonds) > 0
    for hb in hbonds:
        assert abs(hb.donor_residue - hb.acceptor_residue) >= 3
        assert hb.distance < 3.5
        assert hb.energy < 0  # negative = favorable


def test_hbonds_no_close_pairs():
    """Residues far apart should have no H-bonds."""
    residues = []
    for i in range(4):
        x = i * 100.0  # very far apart
        residues.append(_make_residue('ALA', i,
            n=(x, 0, 0), ca=(x+1, 0, 0), c=(x+2, 0, 0), o=(x+2, 1, 0)))
    
    hbonds = calculate_hydrogen_bonds(residues)
    assert len(hbonds) == 0


def test_hbonds_excludes_adjacent():
    """H-bonds should not form between residues closer than 3 in sequence."""
    # Place N and O of adjacent residues very close
    residues = [
        _make_residue('ALA', 0, n=(0, 0, 0), ca=(1, 0, 0), c=(2, 0, 0), o=(2, 1, 0)),
        _make_residue('ALA', 1, n=(0.5, 0, 0), ca=(1.5, 0, 0), c=(2.5, 0, 0), o=(2.5, 1, 0)),
        _make_residue('ALA', 2, n=(1.0, 0, 0), ca=(2.0, 0, 0), c=(3.0, 0, 0), o=(3.0, 1, 0)),
    ]
    hbonds = calculate_hydrogen_bonds(residues, distance_cutoff=5.0)
    # All pairs are within 2 of each other, so no bonds
    assert len(hbonds) == 0


def test_hbonds_custom_cutoff():
    """Custom cutoff should affect number of bonds found."""
    residues = _make_residues_with_coords()
    hbonds_tight = calculate_hydrogen_bonds(residues, distance_cutoff=2.0)
    hbonds_loose = calculate_hydrogen_bonds(residues, distance_cutoff=10.0)
    assert len(hbonds_loose) >= len(hbonds_tight)


def test_hbonds_energy_sign():
    """H-bond energy should be negative (favorable)."""
    residues = _make_residues_with_coords()
    hbonds = calculate_hydrogen_bonds(residues, distance_cutoff=10.0)
    for hb in hbonds:
        assert hb.energy < 0


# --- Test Contact Map ---

def test_contact_map_basic():
    """Test basic contact map generation."""
    residues = _make_residues_with_coords()
    cmap = generate_contact_map(residues, cutoff=8.0)
    assert cmap.size == len(residues)
    assert cmap.cutoff == 8.0
    # Adjacent residues should be in contact
    assert len(cmap.contacts) > 0


def test_contact_map_distance_matrix_symmetry():
    """Distance matrix should be symmetric."""
    residues = _make_residues_with_coords()
    cmap = generate_contact_map(residues)
    for i in range(cmap.size):
        for j in range(cmap.size):
            assert abs(cmap.distance_matrix[i][j] - cmap.distance_matrix[j][i]) < 0.001


def test_contact_map_diagonal_zero():
    """Diagonal of distance matrix should be zero."""
    residues = _make_residues_with_coords()
    cmap = generate_contact_map(residues)
    for i in range(cmap.size):
        assert cmap.distance_matrix[i][i] == 0.0


def test_contact_map_tight_cutoff():
    """Tighter cutoff should give fewer contacts."""
    residues = _make_residues_with_coords()
    cmap_wide = generate_contact_map(residues, cutoff=15.0)
    cmap_tight = generate_contact_map(residues, cutoff=3.0)
    assert len(cmap_tight.contacts) <= len(cmap_wide.contacts)


def test_contact_map_no_coords():
    """Residues without coordinates should produce empty contact map."""
    residues = [_make_residue('ALA', i) for i in range(5)]
    cmap = generate_contact_map(residues)
    assert len(cmap.contacts) == 0


# --- Test Dihedral Angle Calculation ---

def test_dihedral_angle_trans():
    """Trans dihedral (180 degrees) should be ~180 or ~-180."""
    # Proper trans configuration: 4 points in a zigzag in the same plane
    p1 = (0, 1, 0)
    p2 = (0, 0, 0)
    p3 = (1, 0, 0)
    p4 = (1, -1, 0)
    angle = _dihedral_angle(p1, p2, p3, p4)
    assert abs(abs(angle) - 180.0) < 1.0


def test_dihedral_angle_cis():
    """Cis dihedral (0 degrees) should be ~0."""
    p1 = (0, 0, 0)
    p2 = (1, 0, 0)
    p3 = (2, 0, 0)
    p4 = (2, 0, 1)
    angle = _dihedral_angle(p1, p2, p3, p4)
    assert abs(angle) < 1.0 or abs(abs(angle) - 360) < 1.0


def test_compute_dihedral_angles_from_coords():
    """compute_dihedral_angles should set phi/psi from coordinates."""
    residues = _make_residues_with_coords()
    residues = compute_dihedral_angles(residues)
    # At least some angles should be computed
    has_phi = sum(1 for r in residues if r.phi is not None)
    has_psi = sum(1 for r in residues if r.psi is not None)
    assert has_phi > 0
    assert has_psi > 0


# --- Test Utility Functions ---

def test_angle_distance_same_point():
    """Distance from a point to itself should be 0."""
    assert _angle_distance(-57, -47, -57, -47) < 0.01


def test_angle_distance_wrapping():
    """Angle distance should handle wrapping around 180/-180."""
    d = _angle_distance(179, 0, -179, 0)
    assert d < 5.0  # should be ~2 degrees apart, not ~358


def test_distance_3d_basic():
    """Test 3D distance calculation."""
    d = _distance_3d((0, 0, 0), (1, 0, 0))
    assert abs(d - 1.0) < 0.001


def test_distance_3d_symmetric():
    """Distance should be symmetric."""
    a = (1.5, 2.3, 4.7)
    b = (7.1, 3.2, 1.9)
    assert abs(_distance_3d(a, b) - _distance_3d(b, a)) < 0.001


def test_cross_product():
    """Test cross product of unit vectors."""
    x = (1, 0, 0)
    y = (0, 1, 0)
    z = _cross(x, y)
    assert abs(z[0]) < 0.001
    assert abs(z[1]) < 0.001
    assert abs(z[2] - 1.0) < 0.001


def test_dot_product():
    """Test dot product."""
    assert abs(_dot((1, 0, 0), (1, 0, 0)) - 1.0) < 0.001
    assert abs(_dot((1, 0, 0), (0, 1, 0))) < 0.001
    assert abs(_dot((1, 2, 3), (4, 5, 6)) - 32.0) < 0.001


# --- Test Full Analysis Pipeline ---

def test_analyze_structure_with_angles():
    """Full analysis with pre-set angles."""
    residues = _make_mixed_residues()
    result = analyze_structure(residues)
    assert result.num_residues == 10
    assert len(result.secondary_structure) == 10
    assert len(result.ramachandran) == 10
    assert 'ramachandran' in result.summary
    assert 'secondary_structure_counts' in result.summary


def test_analyze_structure_with_coords():
    """Full analysis with coordinates should compute angles and contact map."""
    residues = _make_residues_with_coords()
    result = analyze_structure(residues)
    assert result.num_residues == 8
    assert result.contact_map is not None
    assert result.contact_map.size == 8


def test_analyze_structure_empty():
    """Analysis of empty input should not crash."""
    result = analyze_structure([])
    assert result.num_residues == 0
    assert len(result.secondary_structure) == 0


def test_smooth_secondary_structure():
    """Short runs should be removed."""
    assignments = ['H', 'H', 'E', 'E', 'E', 'H', 'H']
    smoothed = _smooth_secondary_structure(assignments)
    assert smoothed == ['C', 'C', 'E', 'E', 'E', 'C', 'C']


def test_smooth_preserves_long_runs():
    """Runs of 3+ should be preserved."""
    assignments = ['H', 'H', 'H', 'H', 'E', 'E', 'E']
    smoothed = _smooth_secondary_structure(assignments)
    assert smoothed == ['H', 'H', 'H', 'H', 'E', 'E', 'E']


# --- Test CLI ---

def test_cli_rama(tmp_path):
    """Test CLI ramachandran command."""
    import json
    data = [
        {'name': 'ALA', 'index': 0, 'phi': -57.0, 'psi': -47.0},
        {'name': 'VAL', 'index': 1, 'phi': -135.0, 'psi': 135.0},
    ]
    f = tmp_path / 'test.json'
    f.write_text(json.dumps(data))
    
    from evofold_geometry.cli import main
    assert main(['rama', '-i', str(f)]) == 0


def test_cli_rama_json(tmp_path):
    """Test CLI ramachandran with JSON output."""
    import json
    data = [
        {'name': 'ALA', 'index': 0, 'phi': -57.0, 'psi': -47.0},
    ]
    f = tmp_path / 'test.json'
    f.write_text(json.dumps(data))
    
    from evofold_geometry.cli import main
    assert main(['rama', '-i', str(f), '--json']) == 0


def test_cli_ss(tmp_path):
    """Test CLI secondary structure command."""
    import json
    data = [
        {'name': 'ALA', 'index': i, 'phi': -57.0, 'psi': -47.0}
        for i in range(5)
    ]
    f = tmp_path / 'test.json'
    f.write_text(json.dumps(data))
    
    from evofold_geometry.cli import main
    assert main(['ss', '-i', str(f)]) == 0


def test_cli_hbonds(tmp_path):
    """Test CLI hydrogen bond command."""
    import json
    data = []
    for i in range(6):
        x = i * 3.8
        data.append({
            'name': 'ALA', 'index': i,
            'n': [x, 0, 0], 'ca': [x+1, 0, 0],
            'c': [x+2, 0, 0], 'o': [x+2, 1.2, 0],
        })
    f = tmp_path / 'test.json'
    f.write_text(json.dumps(data))
    
    from evofold_geometry.cli import main
    assert main(['hbonds', '-i', str(f)]) == 0


def test_cli_contacts(tmp_path):
    """Test CLI contact map command."""
    import json
    data = []
    for i in range(6):
        x = i * 3.8
        data.append({
            'name': 'ALA', 'index': i,
            'n': [x, 0, 0], 'ca': [x+0.5, 0, 0],
            'c': [x+1.5, 0, 0], 'o': [x+1.5, 1.2, 0],
        })
    f = tmp_path / 'test.json'
    f.write_text(json.dumps(data))
    
    from evofold_geometry.cli import main
    assert main(['contacts', '-i', str(f)]) == 0


def test_cli_analyze(tmp_path):
    """Test CLI full analysis command."""
    import json
    data = [
        {'name': 'ALA', 'index': i, 'phi': -57.0, 'psi': -47.0}
        for i in range(5)
    ]
    f = tmp_path / 'test.json'
    f.write_text(json.dumps(data))
    
    from evofold_geometry.cli import main
    assert main(['analyze', '-i', str(f)]) == 0
