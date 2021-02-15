import pandas as pd
from os import path

model_df = pd.read_csv("../../data/hmm/v2/sequences.csv")
predictions = list(model_df['ID'])

print(f"Starting from {len(predictions)} sequences from our best model")

# Observed residues over UniProt: coverage
# 1 - Identify PDB chain / UniProt mapping using SIFTS data
# http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz

mapping = {}  # { pdbid_chain : [uniprot_accession, ... ] }
with open('../../data/mapping/pdb_chain_uniprot.csv') as f:
    next(f)  # skip first line
    next(f)  # skip second line
    keys = ["pdb", "chain", "uniprot", "res_begin", "res_end",
            "pdb_begin", "pdb_end", "sp_begin", "sp_end"]
    #keys = ["pdb"]
    for line in f:
        values = [value.strip() for value in line.split(",")]
        uniprot_id = values[2]
        #pdb_ids = values[1].split(";")
        # if not values[1] == 'B':
        mapping.setdefault(uniprot_id, [])
        mapping[uniprot_id].append(dict(zip(keys, values)))

unique_pdb_ids = set()
selected_mappings = {}

i = 0
count_mapped_preds = 0

for pred in predictions:
    if pred in mapping:  # check if we have structure data
        count_mapped_preds += 1
        selected_mappings[pred] = mapping[pred]
        for pdb_map in mapping[pred]:
            i += 1
            #print(i, pred, pdb_map)
            unique_pdb_ids.add(pdb_map['pdb'])
    else:
        pass

print(f"We have {len(unique_pdb_ids)} unique pdbs mapped with {count_mapped_preds} mapped SwissProt predictions")

# compute coverage
for up_id, mapping in selected_mappings.items():
    for pdb_map in mapping:
        # sometimes position is not available ('None')
        if 'None' not in [pdb_map['sp_begin'], pdb_map['sp_end'], pdb_map['pdb_begin'], pdb_map['pdb_end']]:
            sp_idxs = set(
                range(int(pdb_map['sp_begin']), int(pdb_map['sp_end'])))
            pdb_idxs = set(
                range(int(pdb_map['pdb_begin']), int(pdb_map['pdb_end'])))
            coverage = len(sp_idxs.intersection(pdb_idxs))/len(sp_idxs)
            pdb_map['coverage'] = coverage
            print(
                f"{pdb_map['pdb']}, chain {pdb_map['chain']} mapping coverage: {coverage}")
        else:
            print(
                f"{pdb_map['pdb']}, {pdb_map['chain']}, [{pdb_map['sp_begin']}, {pdb_map['sp_end']}], [{pdb_map['pdb_begin']}, {pdb_map['pdb_end']}]")

# select coverage > 0.8
coverage_mappings = {}
covered_pdb_ids = set()

for up_id, mapping in selected_mappings.items():
    for pdb_map in mapping:
        if 'coverage' in pdb_map:
            if pdb_map['coverage'] > 0.8:
                pdb_id = pdb_map['pdb']
                # we take always the "highest" chain available (A > B > C)
                if not pdb_id in coverage_mappings or coverage_mappings[pdb_id]['chain'] > pdb_map['chain']:
                    coverage_mappings[pdb_id] = pdb_map
                    covered_pdb_ids.add(pdb_map['pdb'])

print(
    f"There are {len(covered_pdb_ids)} unique pdbs with a coverage above 0.8")

pdb_cut_positions = []
for pdb_id in covered_pdb_ids:
    x = coverage_mappings.get(pdb_id)
    pdb_cut_positions.append(
        [x['pdb'], x['chain'], x['pdb_begin'], x['pdb_end']])
    # print(x['pdb']+":["+x['res_begin']+"," + x['res_end']+"],")
df = pd.DataFrame(data=pdb_cut_positions, columns=[
    'PDB ID', 'Chain', 'Begin', 'End'])
df.to_csv("../../data/structure/pdb_cut_positions.csv", index=False)
df.to_csv("../../data/datasets/family structures/family_structures.csv", index=False)
print("pdb_cut_positions.csv created")
