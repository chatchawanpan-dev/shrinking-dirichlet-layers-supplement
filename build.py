#!/usr/bin/env python3
"""Recompute all numerical illustrations and regenerate their LaTeX artifacts."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
(root / 'manuscript/figures').mkdir(parents=True, exist_ok=True)
(root / 'supplementary').mkdir(exist_ok=True)
for name in ['verify.py', 'continuum_galerkin.py', 'make_table.py', 'make_continuum_artifacts.py']:
    subprocess.run([sys.executable, str(root / 'verification' / name)], cwd=root, check=True)
print('All numerical checks and artifact generation completed.')
