#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
python src/phaseB_closure_analysis.py
python src/make_plots.py
pytest -q
