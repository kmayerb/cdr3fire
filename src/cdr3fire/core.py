"""Core pattern matching functionality."""

import re
import sre_parse
from collections import defaultdict
from typing import Set, List, Tuple
import numpy as np
from scipy.sparse import lil_matrix


def build_3mer_index(strings: List[str]) -> defaultdict:
    """
    Build a k-mer index for 3-mers from a list of strings.

    Parameters
    ----------
    strings : list of str
        List of strings to index.

    Returns
    -------
    index : collections.defaultdict
        Dictionary mapping 3-mers to sets of string indices containing them.
    """
    index = defaultdict(set)
    for i, s in enumerate(strings):
        for j in range(len(s) - 2):
            kmer = s[j:j+3]
            index[kmer].add(i)
    return index


def query_3mers(index: defaultdict, kmers: List[str]) -> Set[int]:
    """
    Query the k-mer index for strings containing all specified k-mers.

    Parameters
    ----------
    index : collections.defaultdict
        k-mer index from build_3mer_index.
    kmers : list of str
        List of 3-mers to search for.

    Returns
    -------
    indices : set of int
        Set of string indices containing all specified k-mers.
    """
    sets = [index.get(kmer, set()) for kmer in kmers]
    return set.intersection(*sets) if sets else set()


def regex_literal_3mers(regex: str) -> Set[str]:
    """
    Extract literal 3-mers from a regex pattern for pre-filtering.

    Parameters
    ----------
    regex : str
        Regular expression pattern.

    Returns
    -------
    kmers : set of str
        Set of 3-mer strings that must be present for the regex to match.
    """
    parsed = sre_parse.parse(regex)
    kmers = set()
    buffer = []

    for token in parsed:
        if token[0] == sre_parse.LITERAL:
            buffer.append(chr(token[1]))
        else:
            if len(buffer) >= 3:
                for i in range(len(buffer) - 2):
                    kmers.add(''.join(buffer[i:i+3]))
            buffer = []
    
    if len(buffer) >= 3:
        for i in range(len(buffer) - 2):
            kmers.add(''.join(buffer[i:i+3]))
    
    return kmers


def match_regexes(regexes: List[str], strings: List[str]):
    """
    Match multiple regex patterns against multiple strings efficiently.

    Parameters
    ----------
    regexes : list of str
        List of regex patterns.
    strings : list of str
        List of strings to match against.

    Returns
    -------
    mat : scipy.sparse.csr_matrix
        Sparse CSR matrix where mat[i, j] = 1 if regexes[i] matches strings[j].
    """
    index = build_3mer_index(strings)
    mat = lil_matrix((len(regexes), len(strings)), dtype=np.uint8)

    for i, regex in enumerate(regexes):
        try:
            kmers = regex_literal_3mers(regex)
            candidate_indices = query_3mers(index, kmers) if kmers else set(range(len(strings)))
            pattern = re.compile(regex)
            for j in candidate_indices:
                if pattern.fullmatch(strings[j]):
                    mat[i, j] = 1
        except Exception as e:
            print(f"Regex error ({regex}): {e}")
            continue

    return mat.tocsr()


def match_regexes_against_cdr3s(regexes: List[str], cdr3_list: List[str]) -> List[Tuple[str, List[int]]]:
    """
    Match regex patterns against CDR3 sequences and return matches with indices.

    Parameters
    ----------
    regexes : list of str
        List of regex patterns.
    cdr3_list : list of str
        List of CDR3 sequences.

    Returns
    -------
    results : list of tuple
        List of tuples (regex, list_of_matching_indices).
    """
    index = build_3mer_index(cdr3_list)
    results = []

    for regex in regexes:
        kmers = regex_literal_3mers(regex)
        candidate_indices = query_3mers(index, kmers) if kmers else set(range(len(cdr3_list)))

        pattern = re.compile(regex)
        matched = [i for i in candidate_indices if pattern.fullmatch(cdr3_list[i])]
        results.append((regex, matched))

    return results
