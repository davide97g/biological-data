# Collect the taxonomic lineage for each protein of the family sequences: (ids.txt)
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd

df = pd.read_csv("../../data/datasets/family sequences/family_sequences.csv")
prot_list = list(df['GO ID'])

taxonomic_lin = {}
result_handle = open(
    '../../data/datasets/family sequences/family_sequences.xml', 'r')
for record in SeqIO.parse(result_handle, 'uniprot-xml'):
    '''    
    seq_record
    record.annotations['taxonomy']'''
    keys = ['Taxonomy_id', 'Taxonomy']
    for rec in record.dbxrefs:
        if rec[:4] == 'NCBI':
            tax_id = rec
    values = [tax_id, record.annotations['taxonomy']]
    taxonomic_lin[record.id] = dict(zip(keys, values))

f = open('../../data/taxonomy/taxonomic_lineage.txt', 'w+')
f.write('Record_ID | Taxonomy | Taxonomy Lineage\n')

for k in taxonomic_lin:
    f.write("{} | {} | {}\n".format(k, taxonomic_lin[k]['Taxonomy_id'].split(
        ':')[1], taxonomic_lin[k]['Taxonomy']))
f.close()

print("taxonomic_lineage.txt created")
