"""The QC finding type and its three severities.

Severity has a precise operational meaning here:

``HARD``
    A schema violation (SCHEMA.md 2 or structural). Blocks admission to
    ``data/curated/``; the CLI exits non-zero.
``WARNING``
    Something a curator should look at: a duplicate, a unit mismatch, a
    missing experimental condition. Does not block anything automatically.
``INFO``
    A *property of the data*, correctly recorded. Conflicts between studies
    and possible non-independence live here. They are findings about the
    literature, not defects in the dataset, and nothing downstream may treat
    them as something to fix.

Python 3.10 compatible, standard library only.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

__all__ = ["QCSeverity", "QCFinding"]


class QCSeverity(str, Enum):
    HARD = "HARD"
    WARNING = "WARNING"
    INFO = "INFO"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return str(self.value)


@dataclass(frozen=True)
class QCFinding:
    """One QC observation, possibly spanning several records."""

    code: str
    severity: QCSeverity
    message: str
    record_ids: List[str] = field(default_factory=list)
    field_name: Optional[str] = None
    details: Dict[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        ids = ", ".join(self.record_ids) if self.record_ids else "-"
        return "[%s] %s: %s  records={%s}" % (
            self.severity.value,
            self.code,
            self.message,
            ids,
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "code": self.code,
            "severity": self.severity.value,
            "message": self.message,
            "record_ids": list(self.record_ids),
            "field": self.field_name,
            "details": dict(self.details),
        }
