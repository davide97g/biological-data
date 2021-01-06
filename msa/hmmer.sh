#!/usr/.binx/env bash

'''
HMMER availability
http://eddylab.org/software/hmmer/hmmer-3.3.1.tar.gz
'''

# Download and compile
# wget -P ../.binx/ http://eddylab.org/software/hmmer/hmmer-3.3.1.tar.gz
# tar -xf ../.binx/hmmer-3.3.1.tar.gz -C ../.binx/
# cd ../.binx/hmmer-3.3.1/ && ./configure && make


# Build a HMM from the MSA with hmmbuild
../.binx/hmmer-3.3.1/src/hmmbuild ../data/msa/davide/tcoffee-3.30.70.100.hmm ../data/msa/davide/tcoffee.fasta

# Perform a custom database search with hmmsearch (use HMMER command line tool)
# ../.binx/hmmer-3.3.1/src/hmmsearch --tblout ../data/msa/tcoffee-3.90.730.10.hmmer_tblout --domtblout ../data/msa/tcoffee-3.90.730.10.hmmer_domtblout ../data/msa/tcoffee-3.90.730.10.hmm ../data/msa/uniprot_sprot.fasta > ../data/msa/tcoffee-3.90.730.10.hmmer_align
