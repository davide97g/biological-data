# setup.py
from urllib import request
import json
import pandas as pd
from os import path

# ? 1. Define your ground truth by finding all proteins in SwissProt annotated (and not annotated) with the assigned Pfam domain
# ? and collect the position of the Pfam domain for all sequences (use the InterPro API).

pfam_id = "PF01582"  # pfam id of our group
url = f"https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/pfam/{pfam_id}?page_size=1000"

req = request.Request(url)
response = request.urlopen(req)
encoded_response = response.read()
decoded_response = encoded_response.decode()
payload = json.loads(decoded_response)

print(f"Found {len(payload['results'])} proteins annotated with {pfam_id}")
with open('../../data/model/domain_positions.json', 'w') as f:
    json.dump(payload, f)
print("domain_positions.json created")

ground_truth_map = {}
ground_truth = []

with open('../../data/model/domain_positions.json') as f:
    domain_positions = json.load(f)
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
                                # do not save twice the same accession id
                                if len(ground_truth_map.get(accession)) == 0:
                                    ground_truth.append(
                                        [accession, True, start, end])
                                ground_truth_map[accession].append(
                                    [start, end])
                        else:
                            print(accession, "no fragments")
                else:
                    print(accession, "no locations")
        else:
            print(accession, "no entries")


# ? get the proteins that do not contain the Pfam domain assigned and are in SwissProt
# ? https://www.uniprot.org/uniprot/?query=NOT%20pf01582&fil=reviewed%3Ayes&sort=score
# ? Download > List > Uncompressed
# ? Add "Accession" as column name + save file as .csv

if path.isfile("../../data/model/uniprot_not_annotated.csv"):
    uniprot_not_annotated = pd.read_csv(
        "../../data/model/uniprot_not_annotated.csv")

    for i in range(len(uniprot_not_annotated)):
        x = uniprot_not_annotated['Accession'][i]
        if ground_truth_map.get(x) is not None:
            print(x, "is a problem")
        else:
            ground_truth.append([x, False, 0, 0])

    # create dataframe from array data and save to .csv file
    gt = pd.DataFrame(data=ground_truth, columns=[
        'Accession ID', 'Annotated', 'start', 'end'])
    gt.to_csv("../../data/model/ground_truth.csv", index=False)
    print("ground_truth.csv created")
else:
    print("Error: uniprot_not_annotated.csv not found.\n\n", "Please visit these link and follow these steps:\n",
          "Open this link https://www.uniprot.org/uniprot/?query=NOT%20pf01582&fil=reviewed%3Ayes&sort=score\n",
          "Download > format 'List' > 'Uncompressed' > download the file into 'project/data/model' folder\n",
          "Add 'Accession' as column name + rename file as uniprot_not_annotated.csv\n",
          "Rerun the script.\n")
