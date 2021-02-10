import pandas as pd
import parse_hmm as hmm

import requests as r
from Bio import SeqIO
from io import StringIO


def getSeq(ID):

    baseUrl = "http://www.uniprot.org/uniprot/"
    currentUrl = baseUrl+ID+".fasta"
    response = r.post(currentUrl)
    cData = ''.join(response.text)

    return list(SeqIO.parse(StringIO(cData), 'fasta'))[0].seq


def downloadAllSeqs(positives):
    # here we download all the sequences from the matches
    seqs = []
    for i in range(len(positives)):
        x = positives.iloc[i]
        accession = x['Accession ID']
        seqs.append([accession, getSeq(accession)])
        print("downloaded", accession)
    print("downloaded", len(positives), "sequences")
    df = pd.DataFrame(data=seqs, columns=[
        'ID', 'Sequence'])
    df.to_csv("../data/sequences.csv", index=False)

# 9.
# Evaluate the ability of matching domain position considering your ground truth, i.e. residues overlapping (and non overlapping) with Pfam domains.


def match(full_map):
    print("\n---------------")
    print("match\n")
    for x in full_map:
        print(full_map.get(x).get("model"),
              " <--> ", full_map.get(x).get("gt"))


# Find all the sequences that are present in the ground truth


def intersect():
    print("\n---------------")
    print("intersect\n")
    # hmm
    pp_map = hmm.build_map()
    # ground truth
    gt = pd.read_csv("../data/ground_truth.csv")
    # extract only the positives
    positives = gt.loc[gt['Annotated'] == True, ]
    # download all the sequences for the positives
    # downloadAllSeqs(positives) # ? only if not already downloaded
    # map of all the matches found and the positions of the residues
    full_map = {}
    # loop
    for i in range(len(positives)):
        x = positives.iloc[i]
        if pp_map.get(x['Accession ID']) is not None:
            full_map[x['Accession ID']] = {
                # position of my model
                'model': pp_map.get(x['Accession ID']),
                'gt': [x['start'], x['end']]  # position of the ground truth
            }
    print("found "+str(len(full_map))+"/"+str(len(pp_map))+" sequences")
    match(full_map)

# Calculate accuracy, precision, sensitivity, specificity, MCC, F-score, etc.


def stastistics():
    print("statistics")


if __name__ == "__main__":
    # getSeq("Q01638")
    intersect()
