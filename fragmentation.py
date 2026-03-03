#!/usr/bin/env python3

import argparse
import pandas as pd
from rdkit import Chem
import re
import sys

def remove_atom_map(smiles):
    """Remove atom map annotations from a SMILES string."""
    return re.sub(r"\[\d+\*\]", "*", smiles)

def neutralize_oxygen(mol):
    """Convert negatively charged oxygen atoms (O-) to neutral O or OH."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:  # O and negative 
            atom.SetFormalCharge(0)  # neutral
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)  # add Oxigen
    return mol

def fragment_molecule(mol):
    """Fragment a molecule based on rotatable (non-ring single) bonds."""
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
    """Check whether a fragment can be further divided."""
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and not bond.IsInRing():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetDegree() > 1 and end_atom.GetDegree() > 1:
                return True  # splitable bond existent
    return False  # no more splitable bonds

def process_smiles_list(smi_list):
    """Process a list of SMILES strings and fragment them."""
    fragment_mols = []
    c = []
    for smi in smi_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:  # mol is not None ...checking 
            mol = neutralize_oxygen(mol)  # fixing O negative 
            fragments = fragment_molecule(mol)
            fragment_mols.extend(fragments)
    
    smi = [Chem.MolToSmiles(frag, canonical=True, isomericSmiles=False) for frag in fragment_mols]
    split = [smi.split('.') for smi in smi]
    split_smiles = [item for sublist in split for item in sublist]
    
    for smiles in split_smiles:
        b = remove_atom_map(smiles)
        c.append(b)
    
    return c

def iterate_fragmentation(initial_smi_list, max_iterations=50):
    """Iteratively fragment molecules until no further divisible fragments remain."""
    smi_list = initial_smi_list
    all_fragments = set()

    for iteration in range(max_iterations):
        print(f"Iteration {iteration + 1}")
        new_fragments = process_smiles_list(smi_list)
        new_fragments_set = set(new_fragments)
        
        # add new fragments
        all_fragments.update(new_fragments_set)

        # nest processing...splitable fragments
        smi_list = [
            smi for smi in (new_fragments_set)
            if is_fragment_divisible(Chem.MolFromSmiles(smi))
        ]

        # no more splitable, finising
        if not smi_list:
            print("No more divisible fragments to process.")
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
