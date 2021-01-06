#!/usr/.binx/env python

import gzip
from Bio import SeqIO
from Bio.PDB import PDBList, is_aa
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqIO.PdbIO import PdbSeqresIterator

# user
user = "davide"

# pdb ids

# ! '1vqs' not mapping to UNIPROT
pdb_ids_davide = ['4npo', '3bm7', '1y0h', '1sqe',
                  '1tz0', '1iuj', '3hx9', '1vqs', '1lq9', '5k9f']
pdb_ids_darko = ['4co3', '4ush', '5d4n', '1uku',
                 '2nuh', '4e98', '1nza', '1p1l', '6gdx', '1naq']
pdb_ids = []
if user == "davide":
    pdb_ids = pdb_ids_davide
else:
    pdb_ids = pdb_ids_darko

coverage_pdb = False
coverage_uniprot = True

for pdb_id in pdb_ids:
    # Fetch a PDB file to the current dir
    # pdbl = PDBList()
    # pdbl.retrieve_pdb_file(pdb_id, pdir='../data/structure/{}/'.format(user),
    #                        file_format='pdb')  # Will save to pdbXXXX.ent

    # Load the structure
    structure = PDBParser(QUIET=True).get_structure(
        pdb_id, "../data/structure/{}/pdb{}.ent".format(user, pdb_id))

    # # Iterate structure
    # for model in structure:
    #     for chain in model:
    #         for residue in chain:
    #             if is_aa(residue):
    #                 for atom in residue:
    #                     print("standard {} {}".format(residue.id, atom.id))
    #             else:
    #                 for atom in residue:
    #                     print("hetero {} {}".format(residue.id, atom.id))

    # Total number of hetero atoms
    # python3 exercise_result.py | grep -c "hetero"
    # 268

    # Total number of water atoms
    # python3 exercise_result.py | grep -c "hetero ('W"
    # 239

    # Total observed residues (count alpha carbons)
    # python3 exercise_result.py | grep "standard" | grep -c " CA"
    # 190

    #########################################
    # Calculate the coverage of observed residues over the SEQRES

    observed = None
    with open("../data/structure/{}/pdb{}.ent".format(user, pdb_id)) as f:
        seq_records = (PdbSeqresIterator(f))

        # Iterate over PDB chains
        for seq_record in seq_records:
            # print(seq_record)
            if seq_record.annotations['chain'] == 'A':

                # Get SEQRES
                # print(seq_record.seq)

                # Get the list of observed residues
                observed = [residue for residue in structure[0]
                            ['A'].get_residues() if is_aa(residue)]
                # print([residue.get_resname() for residue in observed])

                # Print coverage
                if coverage_pdb:
                    print(pdb_id + ' : ' +
                          str(len(observed) / len(seq_record.seq)))

                    # observed over SEQRES = 1.0

    #########################################
    # Observed over UniProt
    # 1 - Identify PDB chain / UniProt mapping using SIFTS data
    # wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz -O ../data/structure/pdb_chain_uniprot.tsv.gz

    mapping = {}  # { pdbid_chain : [uniprot_accession, ... ] }
    with gzip.open("../data/structure/pdb_chain_uniprot.tsv.gz") as f:
        next(f)  # skip first line
        next(f)  # skip second line
        for line in f:
            # print(line)
            pdb, chain, uniprot = line.decode().split("\t")[:3]
            if pdb == pdb_id and chain == "A":
                pdb_chain = "{}_{}".format(pdb, chain)
                mapping.setdefault(pdb_chain, []).append(uniprot)

    # for pdb_chain, uniprot in mapping.items():
    #     print("\n".join(uniprot))

    # 2 - Download corresponding UniProt sequences
    # Go to https://www.uniprot.org/uploadlists/
    # Download format FASTA --> ../data/structure/mapping_uniprot.fasta
    # 3 - Calculate observed / UniProt coverage
    seq_records = list(SeqIO.parse(
        "../data/structure/{}/mapping_uniprot.fasta".format(user), "fasta"))
    uniprot_sequences = {}  # { uniprot: sequence}
    for rec in seq_records:
        uniprot_sequences[rec.id.split("|")[1]] = rec.seq

    # Iterate over uniprot mapping to the PDB
    # The PDB chain can be a chimera, i.e. can map to multiple uniprot
    if coverage_uniprot:
        if mapping.get("{}_A".format(pdb_id)):
            for uniprot in mapping["{}_A".format(pdb_id)]:
                print(uniprot, "{}_A".format(pdb_id), len(observed) /
                      len(uniprot_sequences[uniprot]))  # 0.9947643979057592
        else:
            print("{} : Not found in UNIPROT".format(pdb_id))
