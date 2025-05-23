"""CDR3β sequence generation utilities."""

import random
from typing import Tuple, Dict, List

# Conserved CDR3β segment motifs (simplified)
TRBV_MOTIFS = {
    "TRBV5-1": "CASS",
    "TRBV7-2": "CASR",
    "TRBV20-1": "CASS"
}

TRBD_SEGMENTS = {
    "TRBD1": "GGGGG",
    "TRBD2": "NAGGG"
}

TRBJ_SUFFIXES = {
    "TRBJ1-1": "YEQYF",
    "TRBJ2-3": "FGTQYF",
    "TRBJ2-7": "EQYF"
}

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def random_deletion(seq: str, max_del: int) -> str:
    """
    Trim up to max_del characters from the end of a sequence.
    
    Args:
        seq: Input sequence
        max_del: Maximum number of characters to delete
        
    Returns:
        Sequence with 0 to max_del characters removed from the end
    """
    n = random.randint(0, min(max_del, len(seq)))
    return seq[:-n] if n > 0 else seq


def random_insertion(max_len: int = 5) -> str:
    """
    Insert a random string of amino acids (N/P region).
    
    Args:
        max_len: Maximum length of insertion
        
    Returns:
        Random amino acid string of length 0 to max_len
    """
    n = random.randint(0, max_len)
    return ''.join(random.choices(AMINO_ACIDS, k=n))


def generate_cdr3b_realistic() -> Tuple[str, Dict[str, str]]:
    """
    Generate a realistic CDR3β sequence using V(D)J recombination simulation.
    
    Returns:
        Tuple of (cdr3b_sequence, gene_usage_dict)
    """
    # Randomly select V, D, J genes
    v_name, v_seq = random.choice(list(TRBV_MOTIFS.items()))
    d_name, d_seq = random.choice(list(TRBD_SEGMENTS.items()))
    j_name, j_seq = random.choice(list(TRBJ_SUFFIXES.items()))

    # Simulate exonuclease trimming
    v_trimmed = random_deletion(v_seq, max_del=2)
    d_trimmed = random_deletion(random_deletion(d_seq, max_del=2), max_del=2)
    j_trimmed = j_seq[random.randint(0, 2):]  # trim from front

    # N/P insertions
    n1 = random_insertion(3)
    n2 = random_insertion(3)

    # Assemble final CDR3β
    cdr3b = v_trimmed + n1 + d_trimmed + n2 + j_trimmed
    return cdr3b, {"V": v_name, "D": d_name, "J": j_name}


def generate_random_strings(n: int = 200000, length: int = 10) -> List[str]:
    """
    Generate random amino acid strings.
    
    Args:
        n: Number of strings to generate
        length: Length of each string
        
    Returns:
        List of random amino acid strings
    """
    alphabet = list("ACDEFGHIKLMNPQRSTVWY")
    return [''.join(random.choices(alphabet, k=length)) for _ in range(n)]
