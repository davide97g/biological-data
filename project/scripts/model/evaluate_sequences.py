import pandas as pd
import parse_psiblast as psiblast
import parse_hmm as hmm
import math

if __name__ == "__main__":
    # load GROUND TRUTH data

    df_gt = pd.read_csv("../../data/model/ground_truth.csv")

    # positives
    positives = df_gt.loc[df_gt['Annotated'] == True, ]
    P = len(positives)
    # negatives
    N = len(df_gt)-P


def buildCM(pp_map):

    TP = 0
    for i in range(P):
        x = positives.iloc[i]
        if pp_map.get(x['Accession ID']) is not None:
            TP += 1

    FN = P-TP  # from all the positives, remove the ones correctly classified
    # from the predicted positives, remove the true positives
    FP = len(pp_map)-TP
    TN = N - FN  # the true negatives are the ones not false negatives
    CM = [[TP, FP], [FN, TN]]
    print("\n---------------\nConfusion Matrix")
    for row in CM:
        print(row)
    return CM


def statistics(CM):
    TP = CM[0][0]
    FP = CM[0][1]
    FN = CM[1][0]
    TN = CM[1][1]
    P = TP+FN
    N = FP+TN
    print("\n---------------\nStatistics")

    # sensitivity
    sensitivity = TP/P
    # print("\tsensitivity", sensitivity)

    # specificity
    specificity = TN/N
    # print("\tspecificity", specificity)

    # precision
    precision = TP/(TP+FP)
    # print("\tprecision", precision)

    # accuracy
    accuracy = (TP+TN)/(P+N)
    # print("\taccuracy", accuracy)

    # balanced accuracy
    balanced_accuracy = (sensitivity+specificity)/2
    print("\tbalanced accuracy", round(balanced_accuracy, 3))

    # f1 score
    f1 = 2*TP/(2*TP+FP+FN)
    print("\tf1", round(f1, 3))

    # mcc
    mcc = (TP*TN-FP*FN) / math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    print("\tMCC", round(mcc, 3))


if __name__ == "__main__":
    # hmm
    pp_map_hmm = hmm.build_map()
    CM = buildCM(pp_map_hmm)
    statistics(CM)
    print("---------------\n")

    # psiblast
    pp_map_psiblast = psiblast.build_map()
    CM = buildCM(pp_map_psiblast)
    statistics(CM)
    print("---------------\n")
