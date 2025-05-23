# CDR3 Fire

Efficient regex pattern matching for CDR3β sequences using k-mer indexing acceleration.

## Features

- **Fast Pattern Matching**: Uses 3-mer indexing to accelerate regex matching against large sequence datasets
- **CDR3β Generation**: Realistic synthetic CDR3β sequence generation using V(D)J recombination simulation
- **Sparse Matrix Output**: Memory-efficient sparse matrix representation of match results
- **Command Line Interface**: Easy-to-use CLI for batch processing
- **Python API**: Full programmatic access to all functionality

## Installation

```bash
pip install cdr3fire
```

For development:
```bash
git clone https://github.com/yourusername/cdr3fire.git
cd cdr3fire
pip install -e ".[dev]"
```

## Quick Start

### Command Line Usage

```bash
# Match patterns from CSV files
cdr3fire --patterns patterns.csv --strings sequences.csv \
         --pattern_col regex --string_col cdr3 \
         --output matches.npz
```

### Python API

```python
from cdr3fire import match_regexes, generate_cdr3b_realistic

# Generate some synthetic CDR3β sequences
sequences = []
for _ in range(1000):
    cdr3, genes = generate_cdr3b_realistic()
    sequences.append(cdr3)

# Define some patterns to search for
patterns = [
    "CASS.*YEQYF",  # TRBV5-1 + TRBJ1-1
    "CASR.*EQYF",   # TRBV7-2 + TRBJ2-7
]

# Perform matching
match_matrix = match_regexes(patterns, sequences)

# match_matrix[i,j] = 1 if patterns[i] matches sequences[j]
print(f"Pattern 0 matches {match_matrix[0].sum()} sequences")
```

### Advanced Usage

```python
from cdr3fire import build_3mer_index, query_3mers, regex_literal_3mers

# Build k-mer index for large sequence dataset
sequences = ["CASSYEQYF", "CASRNEQYF", "CASSPYEQYF"]
index = build_3mer_index(sequences)

# Extract literal k-mers from regex for pre-filtering
kmers = regex_literal_3mers("CASS.*YEQYF")  # Returns {"CAS", "ASS"}

# Find candidate sequences
candidates = query_3mers(index, kmers)
print(f"Found {len(candidates)} candidate sequences")
```

## Input Format

### CSV Files
- **Patterns file**: CSV with regex patterns in specified column
- **Strings file**: CSV with sequences in specified column

### Example patterns.csv:
```csv
regex
CASS.*YEQYF
CASR.*EQYF
```

### Example sequences.csv:
```csv
cdr3,frequency
CASSYEQYF,0.001
CASRNEQYF,0.0005
```

## Output

The tool outputs a sparse CSR matrix in `.npz` format where:
- Rows correspond to patterns
- Columns correspond to sequences  
- `matrix[i,j] = 1` if pattern `i` matches sequence `j`

Load the results:
```python
from scipy.sparse import load_npz
match_matrix = load_npz('matches.npz')
```

## Performance

The k-mer indexing provides significant speedup for:
- Large numbers of sequences (>10K)
- Patterns with literal substrings ≥3 characters
- Multiple pattern matching

Typical performance: ~1M sequences × 100 patterns in <30 seconds.

## Development

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Format code
black src/ tests/

# Type checking
mypy src/
```

## License

MIT License - see LICENSE file.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

=== LICENSE ===
MIT License

Copyright (c) 2024 Koshlan Mayer-BLackwell

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=== MANIFEST.in ===
include README.md
include LICENSE
recursive-include src/cdr3fire *.py
recursive-include tests *.py