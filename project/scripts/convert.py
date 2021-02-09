import pandas as pd

hmm_search = pd.read_csv("../data/hmm_search_v2.1.tsv", delimiter="\t")


# for i in range(len(hmm_search)):
# x = hmm_search.iloc[i]
# print(x['Target Accession'])


converted = pd.read_csv("../data/mapping_table_v2.1.csv")
data = []

for i in range(len(hmm_search)):
    start = hmm_search.iloc[i]['Target Ali. Start']
    end = hmm_search.iloc[i]['Target Ali. End']
    accession = converted.iloc[i]['Entry']
    data.append(
        [accession, start, end])

print(data)
df = pd.DataFrame(data=data, columns=[
    'Target Accession', 'start', 'end'])

df.to_csv("../data/hmm_model_hits.csv", index=False)
