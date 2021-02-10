import pandas as pd
import parse_hmm as hmm
import requests as r
from Bio import SeqIO
from io import StringIO
import os.path


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
        seq = full_map.get(x).get("seq")
        model = full_map.get(x).get("model")
        gt = full_map.get(x).get("gt")
        print(model, " - ", gt)
        model_cut = seq[model[0]:model[1]]
        gt_cut = seq[gt[0]:gt[1]]
        print(model_cut)
        print(gt_cut)
        shorter = len(gt_cut)
        if len(model_cut) < shorter:
            shorter = len(model_cut)
        count = 0
        for i in range(shorter):
            s1 = model_cut[i]
            s2 = gt_cut[i]
            if s1 == s2:
                count += 1
        print(count, "matches")

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
    if not os.path.isfile("../data/sequences.csv"):
        downloadAllSeqs(positives)
    # load sequences
    sequences = pd.read_csv("../data/sequences.csv")
    # create map
    seq_map = {}
    for i in range(len(sequences)):
        x = sequences.iloc[i]
        seq_map[x['ID']] = x['Sequence']
    # map of all the matches found and the positions of the residues
    full_map = {}
    # loop
    for i in range(len(positives)):
        x = positives.iloc[i]
        accession = x['Accession ID']
        if pp_map.get(accession) is not None:
            full_map[accession] = {
                'model': pp_map.get(accession),
                'gt': [x['start'], x['end']],
                'seq': seq_map.get(accession)
            }
    print("found "+str(len(full_map))+"/"+str(len(pp_map))+" sequences")
    match(full_map)

# Calculate accuracy, precision, sensitivity, specificity, MCC, F-score, etc.


def stastistics():
    print("statistics")


if __name__ == "__main__":
    # getSeq("Q01638")
    intersect()
