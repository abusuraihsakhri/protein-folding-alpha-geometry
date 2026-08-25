# Protein Structure Analysis

Real implementations of protein structure analysis algorithms: Ramachandran plot classification, secondary structure assignment, hydrogen bond calculation, and contact map generation.

## Features

- **Ramachandran Plot Classification**: Classify phi/psi angles into favored, allowed, and outlier regions for alpha helices, beta sheets, left-handed helices, and polyproline II helices. Glycine gets wider allowed regions.
- **Secondary Structure Assignment**: Assign H (helix), E (strand), T (turn), P (PPII), C (coil) based on dihedral angle criteria with smoothing to remove short runs.
- **Hydrogen Bond Detection**: Calculate backbone N-H···O=C hydrogen bonds using N-O distance cutoff (default 3.5 Å) with energy estimation.
- **Contact Map Generation**: Generate CA-CA contact maps with configurable distance cutoff (default 8.0 Å).
- **Dihedral Angle Computation**: Calculate phi, psi, omega angles from backbone 3D coordinates using the atan2 method.

## Quick Start

```bash
# Ramachandran plot analysis
python cli.py rama -i residues.json

# Secondary structure assignment
python cli.py ss -i residues.json

# Hydrogen bond calculation
python cli.py hbonds -i residues.json --cutoff 3.5

# Contact map generation
python cli.py contacts -i residues.json --cutoff 8.0

# Full analysis
python cli.py analyze -i residues.json

# JSON output (add --json to any command)
python cli.py rama -i residues.json --json
```

## Input Format

JSON file with residue data:

```json
[
  {"name": "ALA", "index": 0, "phi": -57.0, "psi": -47.0},
  {"name": "VAL", "index": 1, "phi": -135.0, "psi": 135.0}
]
```

Or with 3D coordinates (angles computed automatically):

```json
[
  {"name": "ALA", "index": 0, "n": [0,0,0], "ca": [1,0,0], "c": [2,0,0], "o": [2,1.2,0]}
]
```

## Python API

```python
from evofold_geometry import (
    Residue, classify_ramachandran, assign_secondary_structure,
    calculate_hydrogen_bonds, generate_contact_map, analyze_structure,
)

# Create residues with angles
residues = [
    Residue(name='ALA', index=0, phi=-57.0, psi=-47.0),
    Residue(name='VAL', index=1, phi=-135.0, psi=135.0),
]

# Classify Ramachandran angles
for r in residues:
    result = classify_ramachandran(r.name, r.phi, r.psi)
    print(f"{r.name}: {result.region} ({result.structure_type})")

# Full analysis
result = analyze_structure(residues)
print(f"Favored: {result.summary['ramachandran']['favored_pct']:.1f}%")
```

## Ramachandran Regions

| Region | Center (phi, psi) | Tolerance |
|--------|-------------------|-----------|
| Alpha helix | (-57°, -47°) | ±30° (favored), ±50° (allowed) |
| Beta sheet | (-135°, 135°) | ±30° (favored), ±50° (allowed) |
| Left-handed helix | (57°, 47°) | ±30° (favored), ±50° (allowed) |
| Polyproline II | (-75°, 145°) | ±30° (favored), ±50° (allowed) |
| Glycine | any center | ±45° (favored), ±70° (allowed) |

## License

MIT License.
