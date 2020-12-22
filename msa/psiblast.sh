#!/usr/binx/env bash

'''
The statistics of sequence similarity scores: https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html

Blast availability
  Web interface
    NCBI https://blast.ncbi.nlm.nih.gov/Blast.cgi
    EBI https://www.ebi.ac.uk/Tools/sss/ncbiblast/
  Executables
    NCBI ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
'''

# Install blast (linux build 64bit)
wget -P ../binx/ ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
tar -xf ../binx/ncbi-blast-2.11.0+-x64-linux.tar.gz -C ../binx/

# *** Prepare the sequence database ***
# Download SwissProt sequences in the data folder (-P data/)
wget -P ../data/msa/ ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip ../data/msa/uniprot_sprot.fasta.gz

# Create database (Blast indexes, phr+pin+psq files)
../binx/ncbi-blast-2.11.0+/bin/makeblastdb -dbtype prot -in ../data/msa/uniprot_sprot.fasta

# *** Run Blast ***
# Documentation: https://www.ncbi.nlm.nih.gov/books/NBK279690/pdf/Bookshelf_NBK279690.pdf

# Default, print alignemnts on screen
# ../binx/ncbi-blast-2.11.0+/bin/blastp -query ../data/msa/1cudA.fasta -db ../data/msa/uniprot_sprot.fasta

# Tabular output format (-outfmt 6) and E-value threshold 0.001
# ../binx/ncbi-blast-2.11.0+/bin/blastp -query ../data/msa/1cudA.fasta -db ../data/msa/uniprot_sprot.fasta -outfmt 6 -evalue 0.001

# Search the database with a MSA as input
# ../binx/ncbi-blast-2.11.0+/bin/psiblast -in_msa ../data/msa/tcoffee-3.90.730.10.fasta -db ../data/msa/uniprot_sprot.fasta


