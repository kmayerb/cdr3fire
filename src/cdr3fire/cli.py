"""Command line interface for CDR3 matcher."""

import argparse
import pandas as pd
from scipy.sparse import save_npz
from .core import match_regexes


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Match regex patterns against strings with k-mer acceleration"
    )
    parser.add_argument('--patterns', required=True, 
                       help='CSV file containing regex patterns')
    parser.add_argument('--strings', required=True,
                       help='CSV file containing strings')
    parser.add_argument('--pattern_col', required=True,
                       help='Column name for regex patterns')
    parser.add_argument('--string_col', required=True,
                       help='Column name for strings')
    parser.add_argument('--output', required=True,
                       help='Output .npz filename for CSR matrix')
    
    args = parser.parse_args()

    # Load data
    print("Loading data...")
    df_patterns = pd.read_csv(args.patterns)
    df_strings = pd.read_csv(args.strings)

    regexes = df_patterns[args.pattern_col].astype(str).tolist()
    strings = df_strings[args.string_col].astype(str).tolist()

    # Match patterns
    print(f"Matching {len(regexes)} patterns against {len(strings)} strings...")
    mat = match_regexes(regexes, strings)

    # Save results
    save_npz(args.output, mat)
    print(f"Done. Saved sparse matrix to: {args.output}")


if __name__ == "__main__":
    main()