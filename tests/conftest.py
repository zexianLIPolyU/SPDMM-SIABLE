"""
conftest.py
===========
Ensures the local ``spdmm`` package is importable when running pytest
from any directory, without requiring ``pip install -e .``.

pytest auto-discovers and executes this file before any test module is
imported, so it is the canonical place to put sys.path manipulation
that applies to every test.
"""

import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
