# Protein Folding Alpha Geometry

> **Domain:** Computational Biology & AI Drug Discovery
> **Reference Guidelines & Standards:** `wwPDB, IUPAC & CLSI Computational Guidelines`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

**Protein Folding Alpha Geometry** is an advanced analytical and computational platform implementing Invariant Point Attention (IPA) SE(3)-equivariant protein 3D backbone generator. It provides:

- **Protein Structure Analysis**: Ramachandran plot classification, secondary structure assignment, hydrogen bond calculation, and contact map generation
- **Distributed Component Coordination**: Multi-worker evaluation system with consensus-based urgency classification
- **Enterprise Security**: Zero-PHI outbound guard and tamper-evident HMAC-SHA256 audit trail
- **FastAPI REST API**: OpenAPI endpoints for programmatic access

---

## ⚙️ Key Capabilities & Algorithmic Modules

- **Deterministic Calculation Engine**: Strict compliance with standard reference formulations and thresholds
- **Risk & Urgency Classification**: Multi-tier categorization with automated clinical/operational action recommendations
- **Validation & Guardrails**: Rigorous input bounds checking and anomaly detection
- **Ramachandran Plot Analysis**: Classify phi/psi angles into favored/allowed/outlier regions
- **Secondary Structure Assignment**: H/E/T/P/C classification with smoothing
- **Hydrogen Bond Detection**: N-H···O=C backbone hydrogen bond identification
- **Contact Map Generation**: CA-CA distance-based contact mapping

---

## 💻 Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/protein-folding-alpha-geometry.git
cd protein-folding-alpha-geometry

# Install dependencies
pip install -e .

# For development (includes test dependencies)
pip install -e ".[dev]"
```

---

## 🖥️ CLI Usage

### Ramachandran Plot Analysis
```bash
python -m evofold_geometry.cli rama -i residues.json
python -m evofold_geometry.cli rama -i residues.json --json
```

### Secondary Structure Assignment
```bash
python -m evofold_geometry.cli ss -i residues.json
python -m evofold_geometry.cli ss -i residues.json --json
```

### Hydrogen Bond Calculation
```bash
python -m evofold_geometry.cli hbonds -i residues.json --cutoff 3.5
python -m evofold_geometry.cli hbonds -i residues.json --cutoff 3.5 --json
```

### Contact Map Generation
```bash
python -m evofold_geometry.cli contacts -i residues.json --cutoff 8.0
python -m evofold_geometry.cli contacts -i residues.json --cutoff 8.0 --json
```

### Full Structure Analysis
```bash
python -m evofold_geometry.cli analyze -i residues.json
python -m evofold_geometry.cli analyze -i residues.json --hbond-cutoff 3.5 --contact-cutoff 8.0 --json
```

### Input Data Schema

Residues JSON format:
```json
[
  {"name": "ALA", "index": 1, "phi": -57.0, "psi": -47.0},
  {"name": "VAL", "index": 2, "phi": -135.0, "psi": 135.0}
]
```

Or with 3D coordinates:
```json
[
  {"name": "ALA", "index": 1, "n": [x,y,z], "ca": [x,y,z], "c": [x,y,z], "o": [x,y,z]},
  ...
]
```

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`)

### Environment Configuration

Copy `.env.example` to `.env` and configure:

```bash
cp .env.example .env
```

Key variables:
- `AUDIT_SECRET_KEY`: HMAC-SHA256 signing key (generate with `python -c "import secrets; print(secrets.token_hex(32))"`)
- `MODEL_PROVIDER`: LLM provider (`mock`, `ollama`, `claude`, `openai`)

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t protein-folding-alpha-geometry .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY=your-secret-key protein-folding-alpha-geometry
```

Or using Docker Compose:

```bash
# Set your audit key
export AUDIT_SECRET_KEY=your-secret-key

docker-compose up -d
```

---

## 📁 Project Structure

```
protein-folding-alpha-geometry/
├── agents/                  # Enterprise distributed component agents
│   ├── api.py              # FastAPI REST server
│   ├── base.py             # Security, PHI guard, HMAC audit
│   ├── learning.py         # Bayesian calibration engine
│   ├── llm_factory.py      # LLM provider factory
│   ├── metrics.py          # Prometheus metrics exporter
│   ├── models.py           # Pydantic data models
│   ├── streamer.py         # WebSocket telemetry broadcaster
│   ├── supervisor.py       # Master coordinator
│   └── workers.py          # Specialized domain workers
├── evofold_geometry/        # Core protein structure analysis
│   ├── agents.py           # SE(3)-equivariant frame agents
│   ├── cli.py              # Command-line interface
│   ├── engine.py           # Structure analysis engine
│   ├── models.py           # Data model re-exports
│   └── server.py           # FastAPI application factory
├── tests/                   # Test suite
├── web/                     # Operations console (HTML)
├── enrichment.py            # Enrichment feature suite
├── simulator.py             # High-throughput simulation
├── Dockerfile               # Container build
├── docker-compose.yml       # Container orchestration
└── pyproject.toml           # Project configuration
```

---

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.
