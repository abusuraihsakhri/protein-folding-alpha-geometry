# Protein Folding Alpha Geometry

A Python package and browser interface for **protein backbone geometry analysis**. Given residue torsion angles or backbone coordinates, it classifies approximate Ramachandran regions, assigns heuristic secondary structure labels, finds backbone N–O proximity pairs, and generates Cα contact maps.

**Scope:** This code does not implement AlphaFold, Invariant Point Attention (IPA), protein sequence folding, or structure prediction. It is a deterministic geometry analyzer.

## Features

- Ramachandran classification using circular angular distance to reference centers
- Angle-based secondary-structure labels (H, E, T, P, C), with short helix/strand runs smoothed
- N–O proximity pairs and Cα contact maps with configurable distance thresholds
- Computation of missing phi/psi/omega torsions from coordinates
- CLI commands, local Python-in-browser visualization/export, and an optional FastAPI task-audit demonstration

## Browser analyzer

The static interface under the web directory executes the **same Python engine** in Pyodide. Paste residue JSON or upload a JSON file, run analysis, inspect the Ramachandran scatter plot and numerical output, then download the result as JSON.

Build the static site using:

~~~bash
python scripts/build_pages.py
~~~

The Pages workflow publishes the generated dist folder from master when GitHub Pages is enabled. Pyodide is loaded from an external CDN; an internet connection is needed to initialize the Python runtime. User residue input remains in the browser (there is no application-operated upload endpoint). The browser limits JSON to 1 MB and analysis to 1,000 residues. The CLI supports up to 10,000.

## Installation and CLI

Python 3.10+ is required. The geometry CLI has no third-party runtime dependencies.

~~~bash
git clone https://github.com/abusuraihsakhri/protein-folding-alpha-geometry.git
cd protein-folding-alpha-geometry
python -m pip install -e .

python -m evofold_geometry.cli rama -i residues.json --json
python -m evofold_geometry.cli ss -i residues.json --json
python -m evofold_geometry.cli hbonds -i residues.json --cutoff 3.5 --json
python -m evofold_geometry.cli contacts -i residues.json --cutoff 8 --json
python -m evofold_geometry.cli analyze -i residues.json --json
~~~

Minimal example:

~~~json
[
  {"name": "ALA", "index": 1, "phi": -57, "psi": -47},
  {"name": "ALA", "index": 2, "phi": -62, "psi": -42},
  {"name": "ALA", "index": 3, "phi": -54, "psi": -50}
]
~~~

Coordinate input supplies n, ca, c, and o as three-element lists, in Å:

~~~json
[
  {"name": "ALA", "index": 1, "n": [0,0,0], "ca": [1,0,0], "c": [1,1,0], "o": [1,2,0]}
]
~~~

All supplied angles and coordinates must be finite. For the full analysis, supplied angles take precedence, and only missing torsions are calculated from adjacent coordinate records where possible.

## Optional FastAPI server

~~~bash
python -m pip install -e ".[api]"
uvicorn agents.api:app --host 127.0.0.1 --port 8000
~~~

The geometry analyzer is available at http://127.0.0.1:8000/ and the separate task-audit demonstration at http://127.0.0.1:8000/operations.html. API documentation is at http://127.0.0.1:8000/docs.

Endpoints: GET /health, GET /metrics (Prometheus text), POST /api/audit, POST /api/chat, and GET /api/audit/logs. The logs endpoint is disabled unless the server has an AUDIT_LOG_ACCESS_TOKEN; supply that token in an Authorization: Bearer header.

Docker:

~~~bash
export AUDIT_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_hex(32))')"
docker compose up --build
~~~

The audit demonstration is distinct from the geometry calculation. Its threshold-based workers are not a protein prediction system. The only implemented model provider is the deterministic mock. No external model API integration is available.

## Testing

~~~bash
python -m pip install -e ".[api,dev]"
python -m pytest -q
python -m compileall -q agents evofold_geometry scripts
python scripts/build_pages.py
node --check web/app.js
~~~

CI tests Python 3.10, 3.11 and 3.12 and checks the static-site build.

## Scientific and privacy limitations

Ramachandran classifications use simplified angular regions, **not** residue-specific validated probability distributions. Secondary structure is not computed using DSSP. N–O proximity and inverse-distance energy scores do not establish genuine hydrogen bonds. Cα contacts include nearby sequence neighbors. Unknown atom-pair distances in the matrix are represented as zero, not as measured zero-length separation.

The task-audit service has a regex-based sensitive-identifier guard, not comprehensive de-identification or HIPAA certification. Do not submit personal or clinical identifiers. Audit records are held in memory only; HMAC validation detects mutation of retained records but does not anchor the log externally or guarantee persistence across restarts.

The browser is designed for recent WebAssembly-capable desktop and mobile browsers, subject to available device memory and network access for Pyodide.

## License

[MIT](LICENSE).
