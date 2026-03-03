# molecular-fragmentation-rdkit
RDKit-based Rotatable Bond Molecular Fragmentation

## Overview
This project implements a recursive molecular fragmentation workflow based on rotatable bonds using RDKit.

Given a CSV file containing compound IDs and SMILES strings, the script iteratively fragments molecules along rotatable single bonds until no further rotatable bonds remain.

The output is a CSV file containing all generated terminal fragments in canonical SMILES format.

---

## Features

- RDKit-based rotatable bond detection (SMARTS definition)
- Recursive fragmentation until terminal fragments
- Canonical SMILES normalization
- Atom mapping removal
- Command-line interface (CLI) support
- CSV input/output workflow

---

## Requirements

- Python 3.9+
- RDKit
- pandas

Install dependencies:

```bash
pip install pandas
conda install -c rdkit rdkit
```

---

## Input Format

Input CSV must contain the following columns:

| id | mol |
|----|------|
| Compound_1 | CCO |
| Compound_2 | CCC(=O)O |

- `id`: Compound identifier
- `mol`: SMILES string

---

## Usage

```bash
python3 fragmentation.py input.csv output.csv
```

Example:

```bash
python3 fragmentation.py test_compounds.csv output.csv
```

---

## Output

Output CSV format:

| ID | SMILES |
|----|--------|
| Compound_1 | CC |
| Compound_1 | O |
| Compound_2 | CCC |
| Compound_2 | O |

All fragments are returned as canonical SMILES without atom mapping.

---

## Algorithm Overview

1. Parse SMILES into RDKit molecule objects
2. Identify rotatable bonds using RDKit SMARTS definition
3. Fragment molecule along rotatable bonds
4. Normalize fragments to canonical SMILES
5. Repeat until no rotatable bonds remain

---

## Future Improvements

- Parallel processing for large compound libraries
- Visualization of fragmentation trees
- Unit tests with pytest
- PyPI package support

---

## Author

Your Name  
GitHub: https://github.com/yourusername
