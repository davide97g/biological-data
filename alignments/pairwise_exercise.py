#!/usr/bin/env python
from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna, generic_protein
from Bio import SeqIO

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo


# Parse FASTA
seq_records = list(SeqIO.parse("../data/alignments/p53_ensemble.fa", "fasta"))

data = {}
for rec in seq_records:
    # print rec

    # Change alphabet
    # if 'peptide' in rec.description:
    #     rec.seq = Seq(str(rec.seq), generic_protein)
    # else:
    #     rec.seq = Seq(str(rec.seq), generic_dna)

    # Parse the sequence type (peptide, cds, cdna, utr5, utr3, <x>_exon, intron_<x>)
    seq_type = "_".join(rec.description.split(':')[0].split()[1:])
    data[seq_type] = rec

# print("original data sequence")
# print(data)

# *** Find the position of exons in the DNA coding sequence (CDS)

# Find exons against the cds
for k in data.keys():
    if "exon" in k:
        exon_position_cds = data["cds"].seq.find(data[k].seq)
        print("{} {:>5} {:>5} {}...".format(
            k, exon_position_cds, len(data[k].seq), data[k].seq[0:10]))
print()


# *** Why some exons does not align with the cds?

# Find exons against the cdna
for k in data.keys():
    if "exon" in k:
        exon_position_cds = data["cds"].seq.find(data[k].seq)
        exon_position_cdna = data["cdna"].seq.find(data[k].seq)
        print("{} {:>5} {:>5} {}... (length {:>5})".format(
            k, exon_position_cds, exon_position_cdna, data[k].seq[0:10], len(data[k].seq)))
print()

# Find cds inside cdna (show first and last 10 residues of cds aligned with cdna)
# The first and last exons don't contain coding residues.
# The second is partially coding (it starts before the CDS in the cDNA)
cds_position_cdna = data["cdna"].seq.find(data["cds"].seq)
print("cdna_len:{} cds_len:{}".format(
    len(data["cdna"].seq), len(data["cds"].seq)))

print("{:>4} {}...\n   0 {}...".format(cds_position_cdna,
                                       data["cdna"].seq[cds_position_cdna:cds_position_cdna + 10],
                                       data["cds"].seq[:10]))

print("...{}... {}\n...{}    {}".format(data["cdna"].seq[cds_position_cdna + len(data["cds"].seq) - 10:cds_position_cdna + len(data["cds"].seq)],
                                        cds_position_cdna +
                                        len(data["cds"].seq),
                                        data["cds"].seq[-10:], len(data["cds"].seq)))
print()


# *** Translate the exons, identify those which are coding for protein fragments and align against the peptide sequence.
#     Consider the DNA of the exons can have a non-zero phase

alignments = []
for k in data.keys():
    if "exon" in k:
        for i in range(0, 3):
            print("seq[i:]", data[k].seq[i:])
            exon_peptide = data[k].seq[i:].translate(stop_symbol="")

            alignment = pairwise2.align.localds(
                data["peptide"].seq, exon_peptide, MatrixInfo.blosum62, -100, -0.5, one_alignment_only=True)[0]

            # Calculate average score per alignment position
            seq1_aligned, seq2_aligned, score, alignment_begin, alignment_end = alignment
            score_norm = score / float(alignment_end - alignment_begin)
            cov_exon = (alignment_end - alignment_begin) / len(exon_peptide)

            alignments.append((k, score, score_norm, alignment_begin, alignment_end, len(
                exon_peptide), cov_exon, alignment))


for alignment in alignments:
    # Filter out bad alignments, i.e. those with average bit score per position lower than 5. Lower scores
    # indicate there are gaps or unfavored mutations, here we are looking for an exact match.
    if alignment[2] > 5:
        print("{}\n  score:{}\n  score_normalized:{:.2f}\n  peptide_alignment_start:{}\n  peptide_alignment_end:{}\n  exon_len:{}\n  cov_exon:{:.2f}\n".format(
            *alignment[:-1]))
        print(pairwise2.format_alignment(*alignment[-1]))

print()
