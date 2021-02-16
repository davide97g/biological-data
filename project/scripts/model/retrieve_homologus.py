# 2. Retrieve homologous proteins starting from your input sequence a BLAST search against UniProt or UniRef50 or UniRef90

from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import SearchIO
records = []
for seq_record in SeqIO.parse('../../data/model/blast_search.fasta', 'fasta'):
    records.append(seq_record.id)
    print('ID: ', seq_record.id)
    print(seq_record.seq)
    print(repr(seq_record.seq))
    print('Len: ', len(seq_record))

for record in SeqIO.parse("../../data/model/blast_search.fasta", "fasta"):
    print(record.id, len(record))
    result_handle = NCBIWWW.qblast(
        "blastp", "prot", seq_record.seq, megablast=False)
