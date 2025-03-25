import sys
if sys.version_info[0] < 3:
    sys.exit("Error: This script requires Python 3. Please run it with python3.")
    
import argparse
import subprocess
import os
from Bio import PDB
import warnings
from Bio import BiopythonWarning

# Suppress all Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Renumber TCR chains in a TCR or TCR-peptide-MHC complex PDB file using ANARCI (AHo numbering).\n"
                    "This script requires Biopython (version 1.85) and the ANARCI program installed.\n"
                    "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/anarci/\n"
                    "https://biopython.org/\n\n"
                    "The output PDB file will have the following standard chain IDs:\n"
                    "  - MHC-A   -> A\n"
                    "  - MHC-B   -> B\n"
                    "  - Peptide -> C\n"
                    "  - TCR-A   -> D\n"
                    "  - TCR-B   -> E\n",
        epilog="Example usage:\n"
               "  ./renumber_tcr.py -i 3e2h.pdb --tcr_a B --tcr_b C --peptide Q --mhc_a A \n\n"
               "For more tools and TCR-related databases, visit: https://piercelab.ibbr.umd.edu/tools.html \n",
        formatter_class=argparse.RawTextHelpFormatter
    )


    parser.add_argument("-i", "--input", required=True, help="Input PDB file containing the TCR-peptide-MHC complex.")
    parser.add_argument("-o", "--output", default="renumbered_tcr.pdb",
                        help="Output PDB file with renumbered TCR chains (default: renumbered_tcr.pdb).")
    parser.add_argument("--tcr_a", help="Chain ID for TCR-alpha (default: D).")
    parser.add_argument("--tcr_b", help="Chain ID for TCR-beta (default: E).")
    parser.add_argument("--peptide", help="Chain ID for Peptide (default: C).")
    parser.add_argument("--mhc_a", help="Chain ID for MHC-A (default: A).")
    parser.add_argument("--mhc_b", help="Chain ID for MHC-B (default: B, used only for Class II MHC).")

    return parser.parse_args()

# Standard chain IDs
STANDARD_CHAINS = {
    "mhc_a": "A",
    "mhc_b": "B",
    "peptide": "C",
    "tcr_a": "D",
    "tcr_b": "E"
}

# Amino acid three-letter to one-letter conversion dictionary
three_to_one_dict = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

def three_to_one(res_name):
    """Convert three-letter amino acid code to one-letter code."""
    return three_to_one_dict.get(res_name, "X")  # Use 'X' for unknown residues

def parse_anarci_output(output_file):
    """Parse ANARCI output and return AHo numbering, handling insertion codes."""
    aho_mapping = []
    with open(output_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue  # Skip invalid lines
            
            chain_id = parts[0]
            aho_number = int(parts[1])  # Base residue number (integer)
            
            # Handle both cases (with/without insertion code)
            insertion_code = ""
            amino_acid = parts[3] if len(parts) >= 4 else parts[2]
            insertion_code = parts[2] if len(parts) >= 4 else " "
            
            if amino_acid == "-":
                continue  # Ignore gaps

            # Validate insertion_code (only A-Z allowed, otherwise use space)
            if insertion_code and insertion_code not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                insertion_code = " "  # Reset invalid insertion codes to space
                
            # Store tuple with number, insertion code, and amino acid
            aho_mapping.append((aho_number, insertion_code, amino_acid))
    
    return aho_mapping

def renumber_tcr_chain(structure, chain_id, aho_mappings, pdb_sequence):
    """Renumber a single TCR chain using AHo numbering from ANARCI."""
    """Renumber PDB residues and keep only those matching the ANARCI sequence using AHo numbering."""
    if chain_id not in aho_mappings or not aho_mappings[chain_id]:
        print(f"No AHo numbering found for chain {chain_id}. Skipping renumbering.")
        return

    # Extract ANARCI sequence and numbering
    aho_mapping = aho_mappings[chain_id]  # Get the list for this chain
    valid_aho_numbers = [(num, icode) for num, icode, _ in aho_mapping]  # Extract (number, insertion_code)
    anarci_sequence = "".join([aa for _, _, aa in aho_mapping])

    # Find the starting position of the ANARCI sequence within the PDB sequence
    pdb_start_idx = pdb_sequence.find(anarci_sequence)
    if pdb_start_idx == -1:
        print(f"Error: ANARCI sequence not found in PDB sequence for chain {chain_id}!")
        print(f"PDB sequence:   {pdb_sequence}")
        print(f"ANARCI sequence: {anarci_sequence}")
        print(f"Skipping renumbering for chain {chain_id}.")
        return

    pdb_end_idx = pdb_start_idx + len(anarci_sequence)
    print(f"ANARCI sequence found in PDB sequence at positions {pdb_start_idx + 1} to {pdb_end_idx} (1-based indexing).")
    # Iterate over structure to find the target chain
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                # Get all residues in the chain
                pdb_residues = list(chain.get_residues())

                if len(pdb_residues) < pdb_end_idx:
                    print(f"Error: PDB chain {chain_id} has fewer residues ({len(pdb_residues)}) than required ({pdb_end_idx}). Skipping.")
                    return

                # Remove residues before the matching start
                for res in pdb_residues[:pdb_start_idx]:
                    #print(f"Removing N-terminal residue {res.id[1]} from chain {chain_id}")
                    chain.detach_child(res.id)

                # Remove residues after the matching end
                for res in pdb_residues[pdb_end_idx:]:
                    #print(f"Removing C-terminal residue {res.id[1]} from chain {chain_id}")
                    chain.detach_child(res.id)

                # Refresh the residue list after removals
                pdb_residues = list(chain.get_residues())

                # Renumber the remaining residues with AHo numbers and insertion codes
                for residue, (new_num, new_icode) in zip(pdb_residues, valid_aho_numbers):
                    old_id = f"{residue.id[1]}{residue.id[2]}"
                    residue.id = (residue.id[0], new_num, new_icode)  # Update with number and insertion code
                    #print(f"Renumbered residue {old_id} to AHo {new_num}{new_icode} in chain {chain_id}")

                print(f"Chain {chain_id} renumbered with AHo numbering, extra residues removed.")

def renumber_tcr_chains(structure, in_chains, output_pdb):
    """Renumber TCR chains using AHo numbering"""
    # Run ANARCI for TCR chains
    aho_mappings = {}
    pdb_sequences = {}
    for tcr_chain in ["tcr_a", "tcr_b"]:
        chain_id = in_chains.get(tcr_chain)
        if chain_id:
            # Extract sequence
            seq = "".join([three_to_one(res.resname) for model in structure for chain in model if chain.id == chain_id for res in chain if res.id[0] == ' '])
            if not seq:
                print(f"Warning: No TCR sequence found for {tcr_chain} ({chain_id}). Skipping.")
                continue

            pdb_sequences[chain_id] = seq
            # Save sequence to a FASTA file
            fasta_file = f"{tcr_chain}.fasta"
            output_file = f"anarci_output_{tcr_chain}.txt"
            with open(fasta_file, "w") as f:
                f.write(f">{chain_id}\n{seq}\n")

            # Run ANARCI
            subprocess.run(["ANARCI", "-i", fasta_file, "-o", output_file, "-r", "tr", "-s", "a"], check=True)

            # Parse ANARCI output
            aho_mappings[chain_id] = parse_anarci_output(output_file)

    # Renumber TCR chains with sequence validation
    for model in structure:
        for chain in model:
            curr_chain_id = chain.id
            if curr_chain_id in in_chains.values():
                logical_name = next((key for key, value in in_chains.items() if value == curr_chain_id), None)
                if curr_chain_id in [in_chains.get('tcr_a'), in_chains.get('tcr_b')]:
                    print(f"Renumbering {logical_name} chain {curr_chain_id}")
                    renumber_tcr_chain(structure, curr_chain_id, aho_mappings, pdb_sequences.get(curr_chain_id, ""))

    #rename chain id's
    for model in structure:
        for chain in model:
            old_chain_id = chain.id
            # Initialize new_chain_id to the current one, assuming no change
            new_chain_id = old_chain_id
            # Rename chains according to the standard chain ID mapping
            for key, value in STANDARD_CHAINS.items():
                if in_chains.get(key) == old_chain_id:  # Check if current chain id matches the one in Chains
                    new_chain_id = value  # Assign the new standard chain ID
                    break  # Found the correct mapping, no need to check further
            
            # Rename only if the new chain ID is different from the old one
            if old_chain_id != new_chain_id:
                print(f"Renaming chain {old_chain_id} → {new_chain_id}")
                chain.id = new_chain_id  # Update the chain ID to the new one

    # Remove chains that are not in the provided chain list                
    keep_chains = set()
    for chain_key, chain_value in in_chains.items():
        if chain_key in STANDARD_CHAINS:
            keep_chains.add(STANDARD_CHAINS[chain_key])  # Add the chain ID to keep
    for model in structure:
        for chain in list(model):  # Use list(model) to avoid modifying during iteration
            if chain.id not in keep_chains:
                print(f"Removing chain {chain.id}")
                model.detach_child(chain.id)  # Remove unwanted chain

    # Save the renumbered structure
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    print(f"Renumbered PDB saved as {output_pdb}")

def check_chains_in_pdb(input_pdb, provided_chains):
    # Parse the input PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("input_structure", input_pdb)
    # Get a list of all chain IDs in the PDB structure
    chain_ids_in_pdb = {chain.id for model in structure for chain in model}
    # Check for each chain in provided_chains
    for chain_key, chain_value in provided_chains.items():
        if chain_value not in chain_ids_in_pdb:
            raise ValueError(f"Chain {chain_value} (from {chain_key}) not found in input PDB!")

def check_unique_chain_ids(provided_chains):
    if len(set(provided_chains.values())) != len(provided_chains):
        raise ValueError("Chain IDs must be unique.")

def clean_pdb_structure(input_pdb, provided_chains):
    """Remove unwanted chains and HETATM residues from the PDB structure."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("complex", input_pdb)

    for model in structure:
        chains_to_remove = []
        for chain in model:
            if chain.id not in provided_chains.values():
                print(f"Removing chain: {chain.id}")
                chains_to_remove.append(chain.id)
            else:
                # Remove HETATM residues
                residues_to_remove = [res.id for res in chain.get_residues() if res.id[0] != " "]
                for res_id in residues_to_remove:
                    #print(f"Removing HETATM: {chain[res_id].resname} in chain {chain.id}")
                    chain.detach_child(res_id)

        # Remove unwanted chains
        for chain_id in chains_to_remove:
            model.detach_child(chain_id)

    return structure
    
def main():
    args = parse_arguments()
    
    # Collect only explicitly provided chains
    provided_chains = {key: value for key, value in vars(args).items() if key in STANDARD_CHAINS and value is not None}

    # If the user provided no chain options, use all defaults
    if not provided_chains:
        print("No chain options provided. Using default chain IDs:", STANDARD_CHAINS)
        provided_chains = STANDARD_CHAINS.copy()
    else:
        print("Using user-provided chains:", provided_chains)

    # Ensure chain IDs are unique
    check_unique_chain_ids(provided_chains)

    # Ensure chain IDs are present in PDB structure
    check_chains_in_pdb(args.input, provided_chains)

    # **Clean the PDB structure**
    cleaned_pdb = clean_pdb_structure(args.input, provided_chains)

    # Pass the cleaned PDB to renumber_tcr_chains
    renumber_tcr_chains(cleaned_pdb, provided_chains, args.output)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
        
