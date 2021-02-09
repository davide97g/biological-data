import pandas as pd

# hmm_search = pd.read_csv("../data/hmm_search.tsv", delimiter="\t")

# print(hmm_search.head())
# print(len(hmm_search))


# load GROUND TRUTH data

df_gt = pd.read_csv("../data/ground_truth.csv")

print(df_gt.head())
print(len(df_gt))


# ! better generalization of the mapping/conversion process for the ids
# load the mapping for the ids
# df_mapping_ids = pd.read_csv("../data/converted_ids.tsv", delimiter="\t")
# print(df_mapping_ids.head())

# mapping_ids = {}
# for i in range(len(df_mapping_ids)):
#     x = df_mapping_ids.iloc[i]
#     mapping_ids[x['Original']] = x['Converted']

# # create a map for the hmm search results
# pp_map = {}  # predicted positives map
# for i in range(len(hmm_search)):
#     x = hmm_search.iloc[i]
#     pp_map[mapping_ids.get(x['Target Accession'])] = [
#         x['Target Ali. Start'], x['Target Ali. End']]

# print(pp_map.keys())

# ? new version

hmm_model_hits = pd.read_csv("../data/hmm_model_hits.csv")
pp_map = {}
for i in range(len(hmm_model_hits)):
    x = hmm_model_hits.iloc[i]
    pp_map[x['Target Accession']] = [x['start'], x['end']]

print(pp_map.keys())

# positives
positives = df_gt.loc[df_gt['Annotated'] == True, ]
P = len(positives)
print("P", P)
# negatives
N = len(df_gt)-P
print("N", N)

TP = 0
for i in range(P):
    x = positives.iloc[i]
    # print(x['Accession ID'])
    if pp_map.get(x['Accession ID']) is not None:
        TP += 1

FP = P-TP
FN = len(pp_map)-TP
TN = N - FN
CM = [[TP, FP], [FN, TN]]
print("\nConfusion Matrix")
for row in CM:
    print(row)

print("\nStatistics")
# accuracy
accuracy = (TP+TN)/(P+N)
print("accuracy", accuracy)

f1 = 2*TP/(2*TP+FP+FN)
print("f1", f1)
