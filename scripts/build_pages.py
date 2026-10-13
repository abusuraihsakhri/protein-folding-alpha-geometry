"""Assemble the static Pyodide site without duplicating Python source."""
from pathlib import Path
from shutil import copy2, copytree, rmtree

ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "dist"


def build() -> None:
    if DIST.exists():
        rmtree(DIST)
    copytree(ROOT / "web", DIST)
    package = DIST / "evofold_geometry"
    package.mkdir(exist_ok=True)
    for filename in ("__init__.py", "engine.py", "cli.py"):
        copy2(ROOT / "evofold_geometry" / filename, package / filename)
    for path in ("index.html", "app.js", "operations.html",
                 "evofold_geometry/engine.py", "evofold_geometry/cli.py",
                 "evofold_geometry/__init__.py"):
        if not (DIST / path).is_file():
            raise RuntimeError("Missing static-site asset: " + path)


if __name__ == "__main__":
    build()
    print("Static site built at", DIST)
