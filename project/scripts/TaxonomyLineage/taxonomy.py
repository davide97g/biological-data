# Collect the taxonomic lineage for each protein of the family sequences: (ids.txt)
prot_list = []
with open('data/Family Sequences/ids.txt') as f:
    for line in f:
        prot_list = line.strip().split(' ')

prot_list

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
taxonomic_lin = {}
result_handle = open('data/Family Sequences/uniprot-yourlist%3AM20210214E5A08BB0B2D1C45B0C7BC3B55FD265561409D8J.xml', 'r')
for record in SeqIO.parse(result_handle, 'uniprot-xml'):
    '''    
    seq_record
    record.annotations['taxonomy']'''
    keys = ['Taxonomy_id', 'Taxonomy']
    for rec in record.dbxrefs:
        if rec[:4] == 'NCBI':
            tax_id = rec
    values = [tax_id, record.annotations['taxonomy']]
    taxonomic_lin.setdefault(record.id, [])
    taxonomic_lin[record.id].append(dict(zip(keys, values)))
f = open('data/Family Sequences/Taxonomic/taxonomic_lin.txt', 'w+')
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