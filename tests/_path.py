"""Put ``src/`` on ``sys.path`` so the suite runs from a bare checkout.

``python3 -m unittest discover -s tests`` must work with no install step and
no network. Import this module (``import _path  # noqa: F401``) before
importing ``syn4rch``.
"""

import os
import sys

_SRC = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"
)
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
