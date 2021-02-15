
# ? Collect GO annotations for each protein of the family_sequences dataset (entity/dbReference type="GO" in the UniProt XML).
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import json

result_handle = open(
    '../../data/datasets/family sequences/family_sequences.xml', 'r')
go_annotations = {}
go_id_num = {}
for record in SeqIO.parse(result_handle, 'uniprot-xml'):
    seq_record = SeqRecord(seq=record.seq, id=record.id,
                           name=record.name, description=record.description,
                           annotations=record.annotations,
                           dbxrefs=record.dbxrefs)
    rec_id = seq_record.id
    go_id = []
    for ref in seq_record.dbxrefs:
        if ref[:3] == 'GO:':            # if the initial letters of the string matching with 'GO:' then
            go_id.append(ref[3:])
    go_annotations[rec_id] = go_id

json.dumps(go_annotations)

''' No. of GO:ID '''
go_id_num = {}  # {go_id: no. of occurrences}
for id in go_annotations:
    for go_id in go_annotations[id]:
        go_id_num.setdefault(go_id, 0)
        go_id_num[go_id] += 1

# todo: convert to csv
f = open('../../data/function/go_annotations.txt', 'w+')
for k in go_annotations.keys():
    f.write("{} : {}\n".format(k, go_annotations[k]))
f.close()
# todo: convert to csv
f = open('../../data/function/go_annotations_count.txt', 'w+')
for k in go_id_num.keys():
    f.write("{} : {}\n".format(k, go_id_num[k]))
f.close()
# todo: save to csv the ids
