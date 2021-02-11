import pandas as pd
import parse_hmm as hmm
import parse_psiblast as psiblast
import requests as r
from Bio import SeqIO
from io import StringIO
import os.path
from q8 import statistics

# total number of residues in SwissProt
total_residues_swissprot = 5000000


def getSeq(ID):

    baseUrl = "http://www.uniprot.org/uniprot/"
    currentUrl = baseUrl+ID+".fasta"
    response = r.post(currentUrl)
    cData = ''.join(response.text)

    return list(SeqIO.parse(StringIO(cData), 'fasta'))[0].seq


def downloadAllSeqs(data):
    # here we download all the sequences from the matches
    seqs = []
    for i in range(len(data)):
        x = data.iloc[i]
        accession = x['Accession ID']
        seqs.append([accession, getSeq(accession)])
        print("downloaded", accession)
    print("downloaded", len(data), "sequences")
    df = pd.DataFrame(data=seqs, columns=[
        'ID', 'Sequence'])
    df.to_csv("../data/sequences.csv", index=False)

# 9.
# Evaluate the ability of matching domain position considering your ground truth, i.e. residues overlapping (and non overlapping) with Pfam domains.


def intersect(pp_map, model_seqs):
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
    gt_seqs = pd.read_csv("../data/sequences.csv")

    # ? the correct solution is to create two separated maps
    # ? one for the sequences of the model
    # ? the other for the positives sequences of the ground truth
    # ? every map contains as keys the accession ID and as values the tuple:
    # ? (sequence,start,stop)
    model_map = {}
    gt_map = {}
    for i in range(len(model_seqs)):
        x = model_seqs.iloc[i]
        model_map[x['ID']] = {
            'seq': x['Sequence'],
            'start': pp_map.get(x['ID'])[0],
            'stop': pp_map.get(x['ID'])[1]
        }
    for i in range(len(gt_seqs)):
        x = gt_seqs.iloc[i]
        gt_map[x['ID']] = {
            'seq': x['Sequence'],
            'start': positives.iloc[i]['start'],
            'stop': positives.iloc[i]['end']
        }
    # compare
    evaluate(gt_map, model_map)


def evaluate(gt_map, model_map):
    print("\nevaluate\n")
    gt_set = set(gt_map.keys())
    model_set = set(model_map.keys())
    overlapping_set = gt_set & model_set  # ? need to compare residue by residue
    print("overlapping set", len(overlapping_set))
    overlapping_CM = evaluateOverlapping(overlapping_set, gt_map, model_map)
    # ? no need to compare
    false_negatives_set = gt_set - model_set  # ? all false negatives
    false_positives_set = model_set - gt_set  # ? all false positives
    print("false negatives set", len(false_negatives_set))
    print("false positives set", len(false_positives_set))
    FN = 0
    FP = 0
    for x in false_negatives_set:
        FN += len(gt_map.get(x)['seq'])
    for x in false_positives_set:
        FP += len(model_map.get(x)['seq'])

    # update values from overlapping CM
    TP = overlapping_CM[0][0]
    FN += overlapping_CM[1][0]
    FP += overlapping_CM[0][1]
    # ? the negatives are all that remains
    TN = total_residues_swissprot - TP - FN - FP
    # print("TP", TP)
    # print("FN", FN)
    # print("FP", FP)
    # print("TN", TN)
    CM = [[TP, FP], [FN, TN]]
    print("Confusion Matrix")
    for row in CM:
        print(row)
    statistics(CM)


def evaluateOverlapping(overlapping_set, gt_map, model_map):
    CM = [[0, 0], [0, 0]]
    for x in overlapping_set:
        seq = gt_map.get(x)['seq']
        array_gt = []
        array_model = []
        for i in range(len(seq)):
            # ground truth
            if i >= gt_map.get(x)['start'] and i <= gt_map.get(x)['stop']:
                array_gt.append(True)
            else:
                array_gt.append(False)
            # model
            if i >= model_map.get(x)['start'] and i <= model_map.get(x)['stop']:
                array_model.append(True)
            else:
                array_model.append(False)

        # update confusion matrix
        for i in range(len(seq)):
            if array_model[i] == True and array_gt[i] == True:
                CM[0][0] += 1
            elif array_model[i] == True and array_gt[i] == False:
                CM[0][1] += 1
            elif array_model[i] == False and array_gt[i] == True:
                CM[1][0] += 1
            else:
                CM[1][1] += 1
    # return the complete confusion matrix
    return CM


if __name__ == "__main__":
    # hmm
    print("\nHMM")
    pp_map_hmm = hmm.build_map()
    hmm_seqs = hmm.downloadSequences(pp_map_hmm)
    intersect(pp_map_hmm, hmm_seqs)

    # psi blast
    print("\nPSI BLAST")
    pp_map_psiblast = psiblast.build_map()
    psiblast_seqs = psiblast.downloadSequences(pp_map_psiblast)
    intersect(pp_map_psiblast, psiblast_seqs)
