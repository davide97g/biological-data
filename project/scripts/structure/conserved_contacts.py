
from Bio.PDB import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO
import numpy as np

path = "../../data/structure/"
sequence_separation = 6


def create_contact_map(pdb_id, chain, seq_sep=6):
    print(f"Creating contact map for {pdb_id} on chain {chain}")
    structure = PDBParser(QUIET=True).get_structure(
        pdb_id, path+"pdb_cut/pdb{}_cut.ent".format(pdb_id))
    # select chain of first model (0)
    selected_residues = structure[0][chain]
    # Calculate the contact map using the NeighborSearch module (fast)
    # https://biopython.org/docs/1.75/api/Bio.PDB.NeighborSearch.html
    ns = NeighborSearch(
        [atom for residue in selected_residues for atom in residue.get_atoms()])
    # level="R" returns pairs of residues in contact considering all atoms
    contact_map = []
    for residue1, residue2 in ns.search_all(3.5, level="R"):
        if residue1.id[0] == " " and residue2.id[0] == " ":  # Exclude hetero/water residues
            # Sequence separation check
            if abs(residue1.id[1] - residue2.id[1]) >= seq_sep:
                print(residue1.id, residue2.id)


def get_distance_matrix(pdb_id, chain, seq_sep=6):
    structure = PDBParser(QUIET=True).get_structure(
        pdb_id, path+"pdb_cut/pdb{}_cut.ent".format(pdb_id))
    # select chain of first model (0)
    selected_residues = structure[0][chain]
    # Calculate the distance matrix
    distances = []
    for residue1 in selected_residues:
        if residue1.id[0] == " ":  # Exclude hetero/water residues
            row = []
            for residue2 in selected_residues:
                if residue2.id[0] == " ":  # Exclude hetero/water residues
                    if abs(residue1.id[1] - residue2.id[1]) >= seq_sep:
                        row.append(residue1["CA"] - residue2["CA"])
                    else:
                        row.append(None)
            distances.append(row)

    res = np.array(distances, dtype=float)
    np.nan_to_num(res, copy=False, nan=0.0)
    return res


def parse_multiple_structural_alignment():
    print("Parsing multiple structural alignment")
    for seq_record in SeqIO.parse('../../data/structure/multiple_structural_alignment.fasta', 'fasta'):
        print('ID: ', seq_record.id)
        print(seq_record.seq)
        print('Len: ', len(seq_record))

        id_chain = seq_record.id.replace(".pdb", "")
        pdb_id = id_chain[:-1]
        chain = id_chain[-1]
        dm = get_distance_matrix(pdb_id, chain)

        print(f"Distance Matrix {pdb_id} chain {chain} :", dm.shape)


parse_multiple_structural_alignment()
