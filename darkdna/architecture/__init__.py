"""Mode B: sequence-indifferent amount, length, copy, spacing, and occupancy architecture."""

from .comparison import compare_sequence_vs_quantity
from .pipeline import run_sequence_indifferent_architecture, write_architecture_outputs

__all__ = [
    "compare_sequence_vs_quantity",
    "run_sequence_indifferent_architecture",
    "write_architecture_outputs",
]
