"""Tests for sequence generators."""

import pytest
from cdr3fire.generators import (
    random_deletion, random_insertion, generate_cdr3b_realistic,
    generate_random_strings, AMINO_ACIDS
)


def test_random_deletion():
    seq = "ABCDEF"
    result = random_deletion(seq, max_del=2)
    assert len(result) >= len(seq) - 2
    assert len(result) <= len(seq)
    assert seq.startswith(result)


def test_random_insertion():
    insertion = random_insertion(max_len=5)
    assert len(insertion) <= 5
    assert all(aa in AMINO_ACIDS for aa in insertion)


def test_generate_cdr3b_realistic():
    cdr3b, genes = generate_cdr3b_realistic()
    assert isinstance(cdr3b, str)
    assert len(cdr3b) > 0
    assert "V" in genes
    assert "D" in genes
    assert "J" in genes


def test_generate_random_strings():
    strings = generate_random_strings(n=10, length=5)
    assert len(strings) == 10
    assert all(len(s) == 5 for s in strings)
    assert all(all(aa in AMINO_ACIDS for aa in s) for s in strings)