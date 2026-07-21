"""Shared pytest setup: make src/ and the repo root importable, silence RDKit."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
for p in (ROOT, ROOT / "src"):
    sp = str(p)
    if sp not in sys.path:
        sys.path.insert(0, sp)

from rdkit import RDLogger  # noqa: E402
RDLogger.DisableLog("rdApp.*")

FIXTURES = Path(__file__).resolve().parent / "fixtures"
