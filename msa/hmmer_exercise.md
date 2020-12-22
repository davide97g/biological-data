
## Exercises 

Consider the MSA obtained with T-coffee (tcoffee-3.90.730.10.fasta) by aligning 10 domain sequences of the 3.90.730.10 CATH superfamily.


1. Remove empty columns and non-conserved columns at the beginning/end of the MSA using Jalview, rebuild the HMM and repeat the search.
    - How the new alignment performs compared a normal search (tcoffee-3.90.730.10.hmmer_domtblout, generated during practical)? Compare the number of returned hits and the alignment overlap.

2. Perform the same search using Psiblast and the same alignment as input. 
    - How many residues aligned with Psiblast overlap HMMER alignments (consider only significant hits)?

---

To parse domtblout (consider just: target, alifrom, alito, dievalue).
Available fields are:

    target, tacc, tlen, query, qacc, qlen, fevalue, fscore, fbias, _, _, dcevalue, dievalue, dscore, dbias, _, _ , alifrom, alito, evfrom, envto, accuracy, desc = line.strip().split()[:23]

To parse psiblast tabular output (consider just: sseqid, sstart, send, evalue)
Available fields are: 
    
    qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split()
