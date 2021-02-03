# setup.py
import json
import pandas as pd

# 1. Define your ground truth by finding all proteins in SwissProt annotated (and not annotated) with the assigned Pfam domain
# and collect the position of the Pfam domain for all sequences (use the InterPro API).
print("# 1. Defining ground truth")
# pfam id, start, end
ground_truth_map = {}
ground_truth = []

# ? the file domain_position.json was downloaded from this api results
# https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/pfam/PF01582?page_size=1000

with open('../data/domain_positions.json') as f:
    domain_positions = json.load(f)
    # print("count:", domain_positions['count'])
    for r in domain_positions['results']:
        entries = r['entries']
        accession = r['metadata']['accession']
        if ground_truth_map.get(accession) is None:
            ground_truth_map[accession] = []
        if entries:
            for e in entries:
                locations = e['entry_protein_locations']
                if locations:
                    for l in locations:
                        fragments = l['fragments']
                        if fragments:
                            for f in fragments:
                                start = f['start']
                                end = f['end']
                                ground_truth_map[accession].append(
                                    [start, end])
                                ground_truth.append([accession, start, end])
                        else:
                            print(accession, "no fragments")
                else:
                    print(accession, "no locations")
        else:
            print(accession, "no entries")

# ? check results
# print(len(ground_truth_map), ground_truth_map)

# ? check for multiple fragments inside one accession
# for accession in ground_truth_map:
#     if len(ground_truth_map.get(accession)) > 1:
#         print(accession, ground_truth_map.get(accession))

# ? create dataframe from array data and save to .csv file
gt = pd.DataFrame(data=ground_truth, columns=['Accession ID', 'start', 'end'])
gt.to_csv("../data/ground_truth.csv", index=False)
