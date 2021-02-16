# Collect the taxonomic lineage for each protein of the family sequences: (ids.txt)
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd

df = pd.read_csv("../../data/datasets/family sequences/family_sequences.csv")
prot_list = list(df['GO ID'])

taxonomic_lin = {}
result_handle = open(
    '../../data/datasets/family sequences/family_sequences.xml', 'r')

# parse the XML file
for record in SeqIO.parse(result_handle, 'uniprot-xml'):
    keys = ['Taxonomy_id', 'Taxonomy']
    for rec in record.dbxrefs:
        if rec[:4] == 'NCBI':
            tax_id = rec
    values = [tax_id, record.annotations['taxonomy']]
    taxonomic_lin.setdefault(record.id, [])
    taxonomic_lin[record.id].append(dict(zip(keys, values)))

f = open('../../data/taxonomy/taxonomic_lin.txt', 'w+')
for k in taxonomic_lin.keys():
    f.write("{}\n".format(k, taxonomic_lin[k]))
    for key in taxonomic_lin[k][0].keys():
        f.write('{}: \t\t{}\n'. format(key, taxonomic_lin[k][0][key]))
f.close()

for rec in record.dbxrefs:
    if rec[:4] == 'NCBI':
        tax_id = rec
taxonomic_lin[k][0]['Taxonomy_id']
key
