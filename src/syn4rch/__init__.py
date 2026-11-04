"""SYN4RCH: curated hypothalamic neuronal perturbation -> physiology dataset.

A hypothesis-generation platform, not a validated simulator of anything.
See README.md for framing and docs/SCHEMA.md for the binding data contract.

Standard library only; targets Python 3.10.
"""

__version__ = "0.3.0"

#: Version of the Tier A record schema this tree implements (docs/SCHEMA.md).
#:
#: Re-exported, never restated. This constant used to hold its own literal, and
#: at the v0.2.0 bump the two copies diverged: the schema package said 0.2.0
#: while this one still said 0.1.0, and the QC CLI refused to load the dataset
#: it had just written. One source of truth, so that cannot recur.
from .schema.record import SCHEMA_VERSION  # noqa: E402,F401
