from Bio.Seq import Seq
from Bio import SeqIO

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo

'''
TP53 in UniProt
https://www.uniprot.org/uniprot/P04637

Corresponding Fasta sequences from the ENSEMBL database (transcript 1)
https://www.ensembl.org/Homo_sapiens/Transcript/Summary?g=ENSG00000141510;r=17:7668402-7687538;t=ENST00000269305
"Export data" -> "Genomics None" -> "Next" -> "Text":

cDNA = complementary DNA, synthesized from RNA. 5'UTR + CDS + 3'UTR (UTRs regulate expression)
CDS = coding DNA sequence (contains only coding exons, does not contain introns and UTRs)
peptide = translated amino acid sequence
exon = coding CDS fragments (warning! Can have different phase, i.e. the reading frame can either start at 0, 1 or 2)
intron = non coding DNA fragments
utr3 = regulatory DNA region at the end
utr5 = regulatory DNA region at the beginning

'''

# Parse FASTA
seq_records = list(SeqIO.parse("../data/alignments/p53_ensemble.fa", "fasta"))

data = {}  # {molecule_type: [seq_record, ...]}
for rec in seq_records:
    # print(rec)

    # Change alphabet (only for old Python)
    # if 'peptide' in rec.description:
    #     rec.seq = Seq(str(rec.seq), generic_protein)
    # else:
    #     rec.seq = Seq(str(rec.seq), generic_dna)

    # Parse the sequence type (peptide, cds, cdna, utr5, utr3, <x>_exon, intron_<x>)
    seq_type = "_".join(rec.description.split(':')[0].split()[1:])
    data[seq_type] = rec

# for record_type in data:
#     print(record_type, data[record_type], sep="\n")
#     print()
#
#
# # Turn a nucleotide sequence into a protein sequence
print(data["cds"].seq, end="\n\n")
print(data["cds"].seq.transcribe(), end="\n\n")  # DNA to RNA
print(data["cds"].seq.transcribe().translate(),
      end="\n\n")  # DNA to RNA to AA (amino acid)
# DNA to AA (implicit transcribe)
print(data["cds"].seq.translate(), end="\n\n")
# DNA complement (other strand) to AA. Lot of stop codons (*)
print(data["cds"].seq.complement().transcribe().translate(), end="\n\n")


# *** Find the position of exons in the DNA coding sequence (CDS)


# *** Why some exons do not align with the cds?


'''
Align exons to the CDS

Pairwise2 returns the two aligned sequences, the score, the start and the end positions of the alignment

The match parameters are:
CODE  DESCRIPTION
x     No parameters. Identical characters have score of 1, otherwise 0.
m     A match score is the score of identical chars, otherwise mismatch
      score.
d     A dictionary returns the score of any pair of characters.
c     A callback function returns scores.

The gap penalty parameters are:
CODE  DESCRIPTION
x     No gap penalties.
s     Same open and extend gap penalties for both sequences.
d     The sequences have different open and extend gap penalties.
c     A callback function returns the gap penalties.

The align method return a list of alignments as tuples of:
seq1_aligned, seq2_aligned, score, alignment_begin, alignment_end
'''

alignments = []
for k in sorted(data.keys()):
    if "exon" in k:

        # Test different types of alignment
        # alignment = pairwise2.align.globalxx(data["cds"].seq, data[k].seq, one_alignment_only=True)[0]
        # alignment = pairwise2.align.localxx(data["cds"].seq, data[k].seq, one_alignment_only=True)[0]
        # alignment = pairwise2.align.localxs(data["cds"].seq, data[k].seq, -10, -0.5, one_alignment_only=True)[0]
        alignment = pairwise2.align.localds(
            data["cds"].seq, data[k].seq, MatrixInfo.ident, -200, -0.5, one_alignment_only=True)[0]

        # Calculate average score per alignment position
        seq1_aligned, seq2_aligned, score, alignment_begin, alignment_end = alignment
        score_norm = score / float(alignment_end - alignment_begin)

        alignments.append(
            (k, score, score_norm, alignment_begin, alignment_end, alignment))

# Sort by starting position of the alignment
for alignment in sorted(alignments, key=lambda ele: ele[3]):

    print("{} score:{} s_norm:{:.2f} start:{} end:{}".format(*alignment[:5]))
    print(pairwise2.format_alignment(*alignment[-1]))

print()


# *** Translate the exons, identify those which are coding for protein fragments and align against the peptide sequence.
#     Consider the DNA of the exons can have a non-zero phase
