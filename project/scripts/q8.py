import pandas as pd
import parse_psiblast as psiblast
import parse_hmm as hmm

# load GROUND TRUTH data

df_gt = pd.read_csv("../data/ground_truth.csv")

# positives
positives = df_gt.loc[df_gt['Annotated'] == True, ]
P = len(positives)
# negatives
N = len(df_gt)-P


def stastistics(pp_map):

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
    print("\n#####\nConfusion Matrix")
    for row in CM:
        print(row)

    print("\n#####\nStatistics")

    # sensitivity
    sensitivity = TP/P
    print("\tsensitivity", sensitivity)

    # specificity
    specificity = TN/N
    print("\tspecificity", specificity)

    # precision
    precision = TP/(TP+FP)
    print("\tprecision", precision)

    # accuracy
    accuracy = (TP+TN)/(P+N)
    print("\taccuracy", accuracy)

    # balanced accuracy
    balanced_accuracy = (sensitivity+specificity)/2
    print("\tbalanced accuracy", balanced_accuracy)

    # f1 score
    f1 = 2*TP/(2*TP+FP+FN)
    print("\tf1", f1)


# hmm
pp_map_hmm = hmm.build_map()
stastistics(pp_map_hmm)
print("-----------\n")

# psiblast
pp_map_psiblast = psiblast.build_map()
stastistics(pp_map_psiblast)
print("-----------\n")
