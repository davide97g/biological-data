import pandas as pd

version = "v"+input("version ")

# load GROUND TRUTH data

df_gt = pd.read_csv("../data/ground_truth.csv")

# load hmm hits

hmm_model_hits = pd.read_csv("../data/hmm/"+version+"/hmm_hits.csv")
pp_map = {}
for i in range(len(hmm_model_hits)):
    x = hmm_model_hits.iloc[i]
    pp_map[x['Target Accession']] = [x['start'], x['end']]

# print(pp_map.keys())

# positives
positives = df_gt.loc[df_gt['Annotated'] == True, ]
P = len(positives)
# print("P", P)
# negatives
N = len(df_gt)-P
# print("N", N)

TP = 0
for i in range(P):
    x = positives.iloc[i]
    if pp_map.get(x['Accession ID']) is not None:
        TP += 1

FN = P-TP  # from all the positives, remove the ones correctly classified
FP = len(pp_map)-TP  # from the predicted positives, remove the true positives
TN = N - FN  # the true negatives are the ones not false negatives
CM = [[TP, FP], [FN, TN]]
print("\nConfusion Matrix")
for row in CM:
    print(row)

print("\nStatistics")

# accuracy
accuracy = (TP+TN)/(P+N)
print("accuracy", accuracy)

# f1 score
f1 = 2*TP/(2*TP+FP+FN)
print("f1", f1)
