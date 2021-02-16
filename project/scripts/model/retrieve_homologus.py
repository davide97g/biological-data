from Bio import SeqIO

records = []
for seq_record in SeqIO.parse('../../data/model/blast_search.fasta', 'fasta'):
    records.append(seq_record.id)
    print('ID: ', seq_record.id)
    print(seq_record.seq)
    print(repr(seq_record.seq))
    print('Len: ', len(seq_record))

print(f"Blast search found {len(records)} hits")
