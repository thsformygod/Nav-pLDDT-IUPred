# This is basic function to read residues from PDB

from Bio.PDB.MMCIFParser import MMCIFParser

def parse_cif(cif_filename):
    parser = MMCIFParser()
    structure = parser.get_structure("structure", cif_filename)

    # Initialize dictionaries to store information
    chain_sequences = {}
    chain_b_factors = {}

    for model in structure:
        for chain in model:
            chain_id = chain.id
            sequence = ""
            b_factors = {}

            for residue in chain:
                # Extract sequence information
                if residue.get_id()[0] == ' ':
                    sequence += residue.get_resname()

                # Extract B-factor information
                for atom in residue:
                    atom_id = atom.get_id()
                    if atom_id == 'CA':  # Check if missing
                        b_factors[residue.id[1]] = atom.get_bfactor()

            chain_sequences[chain_id] = sequence
            chain_b_factors[chain_id] = b_factors

    return chain_sequences, chain_b_factors

input_example = "D:\\Examples\\02_read_extract\\example\\inputs\\6s7s.cif"


chain_sequences,chain_b_factors = parse_cif(input_example)
print(chain_sequences)
print(chain_b_factors)