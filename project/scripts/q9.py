import pandas as pd
import parse_hmm as hmm
import parse_psiblast as psiblast
import requests as r
from Bio import SeqIO
from io import StringIO
import os.path
from q8 import statistics


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
    global_CM = [[0, 0], [0, 0]]  # every single CM sums into this one
    for x in full_map:
        seq = full_map.get(x).get("seq")
        model = full_map.get(x).get("model")
        gt = full_map.get(x).get("gt")
        # print(gt, " - ", model, " / ", len(seq))
        # for every residues, put 1 or 0 if the residues is present in gt or in model
        array_gt = []
        array_model = []
        for i in range(len(seq)):
            # ground truth
            if i >= gt[0] and i <= gt[1]:
                array_gt.append(True)
            else:
                array_gt.append(False)
            # model
            if i >= model[0] and i <= model[1]:
                array_model.append(True)
            else:
                array_model.append(False)

        # construct confusion matrix
        CM = [[0, 0], [0, 0]]
        for i in range(len(seq)):
            if array_model[i] == True and array_gt[i] == True:
                CM[0][0] += 1
            elif array_model[i] == True and array_gt[i] == False:
                CM[0][1] += 1
            elif array_model[i] == False and array_gt[i] == True:
                CM[1][0] += 1
            else:
                CM[1][1] += 1

        # update global CM
        global_CM[0][0] += CM[0][0]
        global_CM[0][1] += CM[0][1]
        global_CM[1][0] += CM[1][0]
        global_CM[1][1] += CM[1][1]

    print("GLOBAL CONFUSION MATRIX")
    for row in global_CM:
        print(row)
    statistics(global_CM)


def intersect(pp_map):
    print("\n---------------")
    print("intersect\n")
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
    # hmm
    print("\nHMM")
    pp_map_hmm = hmm.build_map()
    intersect(pp_map_hmm)
    print("\nPSI BLAST")
    # psi blast
    pp_map_psiblast = psiblast.build_map()
    intersect(pp_map_psiblast)
