# Collect the taxonomic lineage for each protein of the family sequences: (ids.txt)
prot_list = []
with open('data/Family Sequences/ids.txt') as f:
    for line in f:
        prot_list = line.strip().split(' ')

prot_list

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
taxonomic_lin = {}
result_handle = open('data/Family Sequences/family_sequences.xml', 'r')
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
f = open('data/Family Sequences/Taxonomic/taxonomic_lin.txt', 'w+')
f.write('Record_ID | Taxonomy | Taxonomy Lineage\n')
for k in taxonomic_lin:
    f.write("{} | {} | {}\n".format(k, taxonomic_lin[k]['Taxonomy_id'].split(':')[1], taxonomic_lin[k]['Taxonomy']))
f.close()

for rec in record.dbxrefs:
    if rec[:4] == 'NCBI':
        tax_id = rec
taxonomic_lin[k][0]['Taxonomy_id']
taxonomic_lin[k]['Taxonomy']