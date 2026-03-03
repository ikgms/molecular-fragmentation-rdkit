#!/usr/bin/env python3

import argparse
import pandas as pd
from rdkit import Chem
import re
import sys


def remove_atom_map(smiles: str) -> str:
    """Remove atom mapping from SMILES."""
    return re.sub(r"\[\d+\*\]", "*", smiles)


def neutralize_oxygen(mol: Chem.Mol) -> Chem.Mol:
    """Neutralize negatively charged oxygen atoms."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
    Chem.SanitizeMol(mol)
    return mol

def fragment_molecule(mol: Chem.Mol):
    """Fragment molecule on rotatable bonds."""
    rotatable_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and not bond.IsInRing():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetDegree() > 1 and end_atom.GetDegree() > 1:
                rotatable_bonds.append(bond.GetIdx())
    if rotatable_bonds:
        fragments = Chem.FragmentOnSomeBonds(mol, rotatable_bonds)
        return fragments if isinstance(fragments, tuple) else [fragments]
    else:
        return [mol]


def is_fragment_divisible(mol):
    """Check if fragment still contains rotatable bonds"""
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and not bond.IsInRing():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetDegree() > 1 and end_atom.GetDegree() > 1:
                return True  
    return False


def iterate_fragmentation(initial_smiles, max_iterations=50):
    """Iteratively fragment until no rotatable bonds remain."""
    smi_list = initial_smiles
    all_fragments = set()

    for _ in range(max_iterations):
        new_fragments = []

        for smi in smi_list:
            mol = Chem.MolFromSmiles(smi)
            if not mol:
                continue

            mol = neutralize_oxygen(mol)
            fragments = fragment_molecule(mol)

            for frag in fragments:
                frag_smi = Chem.MolToSmiles(
                    frag, canonical=True, isomericSmiles=False
                )
                frag_smi = remove_atom_map(frag_smi)
                new_fragments.append(frag_smi)

        new_fragments_set = set(new_fragments)
        all_fragments.update(new_fragments_set)

        smi_list = [
            smi
            for smi in new_fragments_set
            if is_fragment_divisible(Chem.MolFromSmiles(smi))
        ]

        if not smi_list:
            break

    return list(all_fragments)


def main():
    parser = argparse.ArgumentParser(
        description="RDKit-based molecular fragmentation using rotatable bonds"
    )
    parser.add_argument("input_csv", help="Input CSV file with 'id' and 'mol' columns")
    parser.add_argument("output_csv", help="Output CSV file")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_csv)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    results = []

    for _, row in df.iterrows():
        smi = row["mol"]
        compound_id = row["id"]

        fragments = iterate_fragmentation([smi])

        for frag in fragments:
            results.append({"ID": compound_id, "SMILES": frag})

    result_df = pd.DataFrame(results)
    result_df.to_csv(args.output_csv, index=False)

    print(f"Fragmentation complete. Results saved to {args.output_csv}")


if __name__ == "__main__":
    main()
