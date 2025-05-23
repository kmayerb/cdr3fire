"""
CDR3 Fire: Efficient regex pattern matching for CDR3β sequences.

This package provides tools for:
- Generating synthetic CDR3β sequences
- Building k-mer indices for fast pattern matching
- Matching regex patterns against large sequence datasets
"""

from .core import (
    build_3mer_index,
    query_3mers,
    regex_literal_3mers,
    match_regexes,
    match_regexes_against_cdr3s,
)

from .generators import (
    generate_cdr3b_realistic,
    generate_random_strings,
    random_deletion,
    random_insertion,
)

__version__ = "0.1.0"
__author__ = "Koshlan Mayer-Blackwell"
__email__ = "your.email@example.com"

__all__ = [
    "build_3mer_index",
    "query_3mers", 
    "regex_literal_3mers",
    "match_regexes",
    "match_regexes_against_cdr3s",
    "generate_cdr3b_realistic",
    "generate_random_strings",
    "random_deletion",
    "random_insertion",
]
