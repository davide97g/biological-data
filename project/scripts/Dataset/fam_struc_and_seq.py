''' Family Struct '''

prediction = {}

with open('data/hmm/v2/hmm_hits.csv') as f:
    next(f)
    keys = ['uniprot_id', 'Start_ali', 'End_ali', 'Start_env', 'End_env']
    for line in f:
        values = [value.strip() for value in line.split(',')]
        uniprot_id = values[0]
        prediction.setdefault(uniprot_id, [])
        prediction[uniprot_id].append(dict(zip(keys, values)))


path = 'data/pdb'
mapping = {}
with open(path + '/pdb_chain_uniprot.csv') as f:
    next(f)
    next(f)
    keys = ['PDB','CHAIN','SP_PRIMARY','RES_BEG','RES_END','PDB_BEG','PDB_END','SP_BEG','SP_END']
    for line in f:
        values = [value.strip() for value in line.split(",")]
        uniprot_id = values[2]
        # pdb_ids = values[1].split(";")
        # if not values[1] == 'B':
        mapping.setdefault(uniprot_id, [])
        mapping[uniprot_id].append(dict(zip(keys, values)))

unique_pdb_ids = set()
selected_mappings = {}

i = 0
count_mapped_preds = 0

for pred in prediction:
  if pred in mapping: # check if we have structure data
    count_mapped_preds += 1
    selected_mappings[pred] = mapping[pred]
    for pdb_map in mapping[pred]:
        i += 1
        #print(i, pred, pdb_map)
        unique_pdb_ids.add(pdb_map['PDB'])
  else:
    pass


for up_id, mapping in selected_mappings.items():
    for pdb_map in mapping:
        # sometimes position is not available ('None')
        if 'None' not in [pdb_map['SP_BEG'], pdb_map['SP_END'], pdb_map['PDB_BEG'], pdb_map['PDB_END']]:
            sp_idxs = set(range(int(pdb_map['SP_BEG']), int(pdb_map['SP_END'])))
            pdb_idxs = set(range(int(pdb_map['PDB_BEG']), int(pdb_map['PDB_END'])))
            coverage = len(sp_idxs.intersection(pdb_idxs))/len(sp_idxs)
            pdb_map['coverage'] = coverage
            print(f"{pdb_map['PDB']}, chain {pdb_map['CHAIN']} mapping coverage: {coverage}")
        else:
            print(f"{pdb_map['PDB']}, {pdb_map['CHAIN']}, [{pdb_map['SP_BEG']}, {pdb_map['SP_END']}], [{pdb_map['PDB_BEG']}, {pdb_map['PDB_END']}]")

coverage_mappings = {}
covered_pdb_ids = set()

for up_id, mapping in selected_mappings.items():
    for pdb_map in mapping:
        if 'coverage' in pdb_map:
            if pdb_map['coverage'] > 0.8:
                pdb_id = pdb_map['PDB']
                # we take always the "highest" chain available (A > B > C)
                if  not pdb_id in coverage_mappings  or  coverage_mappings[pdb_id]['CHAIN'] > pdb_map['CHAIN']:
                    coverage_mappings[pdb_id] = pdb_map
                    covered_pdb_ids.add(pdb_map['PDB'])

print(f"There are {len(covered_pdb_ids)} unique pdbs with a coverage above 0.8")

''' Family Sequences '''

'''
-> All UniRef90 sequences matching our best model (hmm)
'''
ids = {}
with open('data/hmm/v2/sequences.csv') as f:
    next(f)
    keys = ['ID', 'SEQUENCE', 'START', 'STOP']
    for line in f:
        values = [value.strip() for value in line.split(',')]
        id = values[0]
        ids.setdefault(id, [])
        ids[id].append(dict(zip(keys, values)))
with open ('data/Family Sequences/ids.txt','w') as f:
    for id in ids.keys():
        f.write(id+' ')