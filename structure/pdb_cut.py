#!/usr/.binx/env python
'''
Basic parsing and iteration of Structure objects implemented by the BIO module of BioPython
Save selected residues into a new PDB file
Generate distance matrix

Bio.PDB module FAQ
https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
'''


from Bio.PDB import PDBList, is_aa, PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import IUPACData
from Bio.PDB.PDBIO import Select
from Bio.SeqIO.PdbIO import PdbSeqresIterator

# user
user = "darko"

# domains =

domains_davide = ['4npo', '3bm7', '1y0h', '1sqe',
                  '1tz0', '1iuj', '3hx9', '1vqs', '1lq9', '5k9f']

domains_darko = ['4co3', '4ush', '5d4n', '1uku',
                 '2nuh', '4e98', '1nza', '1p1l', '6gdx', '1naq']

domains = []

if user == 'darko':
    domains = domains_darko
else:
    domains = domains_davide

for domain in domains:
    # Input
    pdb_id = domain

    # Fetch a PDB file to the current dir
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, pdir='../data/structure/'+user,
                           file_format='pdb')  # Will save to pdbXXXX.ent

    # Get the SEQRES
    # with open("../data/structure/pdb{}.ent".format(pdb_id)) as f:
    #     seq_records = (PdbSeqresIterator(f))
    #     for seq_record in seq_records:
    #         print(seq_record)

    # Load the structure
    structure = PDBParser(QUIET=True).get_structure(
        pdb_id, "../data/structure/{}/pdb{}.ent".format(user, pdb_id))

    # Iterate structure
    # for model in structure:
    #     for chain in model:
    #         for residue in chain:
    #             if is_aa(residue):  # Filter hetero groups (returns only amino acids)
    #                 # residue.id tuple contains hetero_flag, position, insertion_code
    #                 print("model {} chain {} residue_id {} resname {} resname_3to1 {}".format(model.id, chain.id, residue.id, residue.get_resname(),
    #                                                        IUPACData.protein_letters_3to1.get(residue.get_resname().capitalize())))
    #                 for atom in residue:
    #                     print("atom {} {} {}".format(atom.id, atom.get_bfactor(), atom.get_coord()))

    # Extract a list of residues between start and end positions excluding hetero and water atoms
    # It assumes there are not insertion codes and residues numbers are not necessarily continuous
    residues = []
    start_flag = False
    domain_start = (" ", 17, " ")
    domain_end = (" ", 118, " ")
    total_residues = 0
    for residue in structure[0]['A'].get_residues():  # Model 0, chain A
        if residue.id[0] == " ":  # Exclude hetero and water atoms
            total_residues += 1
            # print(residue.id)
            # Get starting position
            if residue.id == domain_start:
                start_flag = True

            if start_flag:
                residues.append(residue)
                # print(residue.id)

            # Get ending position
            if residue.id == domain_end:
                break
    print(domain)
    print(total_residues)
    print(len(residues))
    print("---")
    # Save a PDB chain

    class Select(Select):
        def __init__(self, chain_ids=None, residues=None):
            self.chain_ids = chain_ids
            self.residues = residues

        def accept_chain(self, chain):
            return self.chain_ids is None or chain.id in self.chain_ids

        def accept_residue(self, residue):
            return self.residues is None or residue in self.residues

    # Save a PDB file containing only a list of selected residues
    io = PDBIO()
    io.set_structure(structure[0])
    io.save("../data/structure/{}/pdb{}_cut.ent".format(user, pdb_id),
            select=Select(residues=residues))
