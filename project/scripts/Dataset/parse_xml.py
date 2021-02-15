from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
result_handle = open('data/Family Sequences/uniprot-yourlist%3AM20210214E5A08BB0B2D1C45B0C7BC3B55FD265561409D8J.xml', 'r')
records = []
seq_records = []
for record in SeqIO.parse(result_handle, 'uniprot-xml'):
    records.append(record)
    seq_record = SeqRecord(seq=record.seq, id=record.id,
                            name=record.name, description=record.description,
                            annotations=record.annotations,
                            dbxrefs = record.dbxrefs)
    seq_records.append(seq_record)
go_annotations = {}
'''
For each sequence record in lis 'seq_record' 
add the record id and GO annotations
'''
go_annotations = {} # {sequence_id: go_ids}
for seq_record in seq_records:
    rec_id = seq_record.id

    go_id = []
    for ref in seq_record.dbxrefs:
        if ref[:3] == 'GO:' :            # if the initial letters of the string matching with 'GO:' then
            go_id.append(ref[3:])
    go_annotations[rec_id] = go_id

''' No. of GO:ID '''
go_id_num = {} # {go_id: no. of occurrences}
for id in go_annotations:
    for go_id in go_annotations[id]:
        go_id_num.setdefault(go_id, 0)
        go_id_num[go_id] += 1
