"""Tests for core functionality."""

import pytest
from cdr3fire.core import (
    build_3mer_index, query_3mers, regex_literal_3mers, 
    match_regexes, match_regexes_against_cdr3s
)


def test_build_3mer_index():
    strings = ["ABCDE", "CDEFG", "XYZAB"]
    index = build_3mer_index(strings)
    assert "ABC" in index
    assert 0 in index["ABC"]
    assert "CDE" in index
    assert {0, 1} == index["CDE"]


def test_query_3mers():
    strings = ["ABCDE", "CDEFG", "ABCFG"]
    index = build_3mer_index(strings)
    result = query_3mers(index, ["ABC"])
    assert result == {0, 2}
    result = query_3mers(index, ["ABC", "CDE"])
    assert result == {0}


def test_regex_literal_3mers():
    regex = "ABCDEF"
    kmers = regex_literal_3mers(regex)
    expected = {"ABC", "BCD", "CDE", "DEF"}
    assert kmers == expected


def test_match_regexes():
    regexes = ["ABC.*", ".*DEF"]
    strings = ["ABCXYZ", "XYZDEF", "ABCDEF", "XYZ"]
    mat = match_regexes(regexes, strings)
    
    # Check matches
    assert mat[0, 0] == 1  # "ABC.*" matches "ABCXYZ"
    assert mat[0, 2] == 1  # "ABC.*" matches "ABCDEF"
    assert mat[1, 1] == 1  # ".*DEF" matches "XYZDEF"
    assert mat[1, 2] == 1  # ".*DEF" matches "ABCDEF"


def test_match_regexes_against_cdr3s():
    regexes = ["CASS.*", ".*YEQYF"]
    cdr3s = ["CASSXYZF", "XYZTYEQYF", "CASSYEQYF", "RANDOM"]
    results = match_regexes_against_cdr3s(regexes, cdr3s)
    
    assert len(results) == 2
    assert results[0][0] == "CASS.*"
    assert set(results[0][1]) == {0, 2}  # indices of matches
    assert results[1][0] == ".*YEQYF"
    assert set(results[1][1]) == {1, 2}