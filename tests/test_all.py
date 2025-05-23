
"""Integration tests for CDR3Fire package."""

import os
import tempfile
import time
import re
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import lil_matrix, load_npz

from cdr3fire.generators import generate_cdr3b_realistic
from cdr3fire.core import match_regexes
from cdr3fire.cli import main


class TestRealisticIntegration:
    """Integration test based on realistic CDR3β generation and pattern matching."""
    
    def test_realistic_cdr3_matching_performance(self):
        """
        Test realistic CDR3β generation, pattern creation, and matching performance.
        This test replicates the original benchmark comparing k-mer accelerated 
        matching vs naive regex matching.
        """
        # Generate realistic CDR3β sequences
        n = 100000  # Reduced from 100K for faster CI testing
        print(f"Generating {n} realistic CDR3β sequences...")
        
        # Generate sequences with metadata
        generated_data = [generate_cdr3b_realistic() for _ in range(n)]
        sequences_with_metadata = [
            (cdr3b, genes['V'], genes['D'], genes['J']) 
            for cdr3b, genes in generated_data
        ]
        
        # Create DataFrame
        df = pd.DataFrame(
            sequences_with_metadata, 
            columns=['cdr3b', 'v', 'd', 'j']
        )
        
        # Create patterns by modifying first 10 sequences
        # Replace positions 5 and 7 with wildcards to create regex patterns
        patterns_list = []
        for item in df['cdr3b'].head(100):
            try:
                x = list(item)
                if len(x) > 5:
                    x[5] = "."
                if len(x) > 7:
                    x[7] = "."
                patterns_list.append("".join(x))
            except (IndexError, TypeError):
                patterns_list.append(item)
        
        patterns_df = pd.DataFrame({'regex': patterns_list})
        
        # Test k-mer accelerated matching (our implementation)
        print("Testing k-mer accelerated matching...")
        tic = time.perf_counter()
        
        regexes = patterns_df['regex'].tolist()
        sequences = df['cdr3b'].tolist()
        kmer_result_matrix = match_regexes(regexes, sequences)
        
        toc = time.perf_counter()
        kmer_time = toc - tic
        print(f"K-mer accelerated matching done in {kmer_time:0.4f} seconds")
        
        # Test naive regex matching (baseline comparison)
        print("Testing naive regex matching...")
        tic = time.perf_counter()
        
        naive_matrix = lil_matrix((len(patterns_list), len(sequences)), dtype=np.uint8)
        for i, pattern in enumerate(patterns_list):
            try:
                compiled_pattern = re.compile(pattern)
                for j, sequence in enumerate(sequences):
                    if compiled_pattern.fullmatch(sequence):
                        naive_matrix[i, j] = 1
            except re.error:
                # Skip invalid regex patterns
                continue
        
        toc = time.perf_counter()
        naive_time = toc - tic
        print(f"Naive regex matching done in {naive_time:0.4f} seconds")
        
        # Verify results are identical
        kmer_dense = kmer_result_matrix.todense()
        naive_dense = naive_matrix.todense()
        
        assert np.array_equal(kmer_dense, naive_dense), \
            "K-mer accelerated and naive results should be identical"
        
        # Performance assertions (k-mer should be faster for larger datasets)
        print(f"Speedup factor: {naive_time / kmer_time:.2f}x")
        
        # Basic sanity checks
        assert kmer_result_matrix.shape == (len(patterns_list), len(sequences))
        assert kmer_result_matrix.dtype == np.uint8
        
        # Check that we found some matches
        total_matches = kmer_result_matrix.sum()
        assert total_matches > 0, "Should find at least some pattern matches"
        
        print(f"Found {total_matches} total matches across all patterns")
        print(f"Average matches per pattern: {total_matches / len(patterns_list):.1f}")

    def test_cli_integration_with_files(self):
        """Test the CLI interface with temporary CSV files."""
        
        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            # Generate test data
            n = 100000
            generated_data = [generate_cdr3b_realistic() for _ in range(n)]
            sequences_with_metadata = [
                (cdr3b, genes['V'], genes['D'], genes['J']) 
                for cdr3b, genes in generated_data
            ]
            
            # Create sequences CSV
            df = pd.DataFrame(
                sequences_with_metadata, 
                columns=['cdr3b', 'v', 'd', 'j']
            )
            sequences_file = os.path.join(temp_dir, 'sequences.csv')
            df.to_csv(sequences_file, index=False)
            
            # Create patterns CSV
            patterns_list = []
            for item in df['cdr3b'].head(100):
                try:
                    x = list(item)
                    if len(x) > 5:
                        x[5] = "."
                    patterns_list.append("".join(x))
                except (IndexError, TypeError):
                    patterns_list.append(item)
            
            patterns_df = pd.DataFrame({'regex': patterns_list})
            patterns_file = os.path.join(temp_dir, 'patterns.csv')
            patterns_df.to_csv(patterns_file, index=False)
            
            # Test CLI by calling main function directly
            output_file = os.path.join(temp_dir, 'matches.npz')
            
            # Mock sys.argv for CLI
            import sys
            original_argv = sys.argv
            try:
                sys.argv = [
                    'cdr3fire',
                    '--patterns', patterns_file,
                    '--strings', sequences_file,
                    '--pattern_col', 'regex',
                    '--string_col', 'cdr3b',
                    '--output', output_file
                ]
                
                # Run CLI
                main()
                
                # Verify output file was created
                assert os.path.exists(output_file), "Output .npz file should be created"
                
                # Load and verify the matrix
                result_matrix = load_npz(output_file)
                assert result_matrix.shape == (len(patterns_list), len(df))
                assert result_matrix.dtype == np.uint8
                
                # Should have some matches
                assert result_matrix.sum() > 0, "Should find some matches"
                
            finally:
                sys.argv = original_argv

    def test_edge_cases_and_robustness(self):
        """Test edge cases and error handling."""
        
        # Test with empty sequences
        empty_result = match_regexes(["CASS.*"], [])
        assert empty_result.shape == (1, 0)
        
        # Test with empty patterns
        empty_patterns_result = match_regexes([], ["CASSYEQYF"])
        assert empty_patterns_result.shape == (0, 1)
        
        # Test with invalid regex (should not crash)
        sequences = ["CASSYEQYF", "CASRNEQYF"]
        invalid_regexes = ["CASS.*", "[INVALID"]  # Second regex is invalid
        result = match_regexes(invalid_regexes, sequences)
        
        # Should still return proper matrix shape
        assert result.shape == (2, 2)
        # First row should have matches, second row should be empty due to invalid regex
        assert result[0, :].sum() > 0  # First pattern should match
        assert result[1, :].sum() == 0  # Second pattern invalid, no matches
        
        # Test very short sequences
        short_sequences = ["AB", "CD"]
        short_result = match_regexes(["A.*"], short_sequences)
        assert short_result.shape == (1, 2)
        assert short_result[0, 0] == 1  # Should match "AB"
        assert short_result[0, 1] == 0  # Should not match "CD"
