import pandas as pd


if __name__ == "__main__":
    version = "v"+input("hmm version: ")
    hmm_search = pd.read_csv("../data/hmm/"+version +
                             "/hmm_search.tsv", delimiter="\t")

    # todo: check if file "mapping" exists

    converted = pd.read_csv("../data/hmm/"+version+"/mapping.csv")
    data = []

    for i in range(len(hmm_search)):
        start = hmm_search.iloc[i]['Target Ali. Start']
        end = hmm_search.iloc[i]['Target Ali. End']
        accession = converted.iloc[i]['Entry']
        data.append(
            [accession, start, end])

    # print(data)
    df = pd.DataFrame(data=data, columns=[
        'Target Accession', 'start', 'end'])

    df.to_csv("../data/hmm/"+version+"/hmm_hits.csv", index=False)
    print("hmm hits saved")


def build_map():
    version = "v"+input("hmm version: ")
    hmm_model_hits = pd.read_csv("../data/hmm/"+version+"/hmm_hits.csv")
    pp_map = {}
    for i in range(len(hmm_model_hits)):
        x = hmm_model_hits.iloc[i]
        pp_map[x['Target Accession']] = [x['start'], x['end']]
    return pp_map
