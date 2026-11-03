"""Guarded detection of the optional third-party stack.

``numpy`` and ``scikit-learn`` are permitted in this layer but **nothing here
requires them**. Every number produced anywhere under ``syn4rch.models`` is
computed with the standard library only, so a metric computed on a machine
without numpy is bit-for-bit the same metric as one computed with it. These
imports exist so that:

* ``import syn4rch.models`` can never fail because an optional wheel is
  missing (every import below is inside ``try``), and
* the benchmark report can state which optional libraries were present when
  it was written, which is provenance, not a dependency.

The alternative -- computing accuracy with numpy and a fallback in pure
Python -- would mean two implementations that have to be proven equal. One
implementation that is always used is cheaper and more honest.

Python 3.10 compatible.
"""

from typing import Any, Dict, Optional

__all__ = [
    "HAVE_NUMPY",
    "HAVE_SKLEARN",
    "NUMPY_VERSION",
    "SKLEARN_VERSION",
    "numpy_module",
    "sklearn_module",
    "optional_dependency_report",
]

_np: Optional[Any] = None
_NUMPY_VERSION: Optional[str] = None
try:  # pragma: no cover - exercised by whichever branch the machine is in
    import numpy as _np_mod  # type: ignore

    _np = _np_mod
    _NUMPY_VERSION = getattr(_np_mod, "__version__", "unknown")
except Exception:  # pragma: no cover - ImportError, or a broken build
    _np = None
    _NUMPY_VERSION = None

_sk: Optional[Any] = None
_SKLEARN_VERSION: Optional[str] = None
try:  # pragma: no cover
    import sklearn as _sk_mod  # type: ignore

    _sk = _sk_mod
    _SKLEARN_VERSION = getattr(_sk_mod, "__version__", "unknown")
except Exception:  # pragma: no cover
    _sk = None
    _SKLEARN_VERSION = None

#: True when ``numpy`` imported. No metric depends on this flag.
HAVE_NUMPY = _np is not None
#: True when ``scikit-learn`` imported. No metric depends on this flag.
HAVE_SKLEARN = _sk is not None

NUMPY_VERSION = _NUMPY_VERSION
SKLEARN_VERSION = _SKLEARN_VERSION


def numpy_module() -> Optional[Any]:
    """The imported ``numpy`` module, or ``None``. Never raises."""
    return _np


def sklearn_module() -> Optional[Any]:
    """The imported ``sklearn`` module, or ``None``. Never raises."""
    return _sk


def optional_dependency_report() -> Dict[str, Optional[str]]:
    """What was importable, for the provenance block of a report."""
    return {
        "numpy": NUMPY_VERSION,
        "scikit-learn": SKLEARN_VERSION,
        "metrics_implementation": "stdlib (no optional library affects any number)",
    }
