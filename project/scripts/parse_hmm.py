import pandas as pd
import os.path

if __name__ == "__main__":
    version = "v"+input("hmm version: ")
    hmm_search = pd.read_csv("../data/hmm/"+version +
                             "/hmm_search.tsv", delimiter="\t")

    if os.path.isfile("../data/hmm/"+version+"/mapping.csv"):
        converted = pd.read_csv("../data/hmm/"+version+"/mapping.csv")
        data = []

        for i in range(len(hmm_search)):
            start_ali = hmm_search.iloc[i]['Target Ali. Start']
            end_ali = hmm_search.iloc[i]['Target Ali. End']
            start_env = hmm_search.iloc[i]['Target Env. Start']
            end_env = hmm_search.iloc[i]['Target Env. End']
            accession = converted.iloc[i]['Entry']
            data.append(
                [accession, start_ali, end_ali, start_env, end_env])

        # print(data)
        df = pd.DataFrame(data=data, columns=[
            'Target Accession', 'Start Ali.', 'End Ali.', 'Start Env.', 'End Env.'])

        df.to_csv("../data/hmm/"+version+"/hmm_hits.csv", index=False)
        print("hmm hits saved")
    elif os.path.isfile("../data/hmm/"+version+"/mapping.tsv"):

        df_mapping = pd.read_csv(
            "../data/hmm/"+version+"/mapping.tsv", delimiter="\t")

        data = []

        mapping = {}
        for i in range(len(df_mapping)):
            mapping[df_mapping.iloc[i]['Target Accession']
                    ] = df_mapping.iloc[i]['Entry']

        for i in range(len(hmm_search)):
            x = hmm_search.iloc[i]
            start_ali = x['Target Ali. Start']
            end_ali = x['Target Ali. End']
            start_env = x['Target Env. Start']
            end_env = x['Target Env. End']
            accession = mapping.get(x['Target Accession'])
            data.append(
                [accession, start_ali, end_ali, start_env, end_env])

        # print(data)
        df = pd.DataFrame(data=data, columns=[
            'Target Accession', 'Start Ali.', 'End Ali.', 'Start Env.', 'End Env.'])

        df.to_csv("../data/hmm/"+version+"/hmm_hits.csv", index=False)
        print("hmm hits saved")
    else:
        print("File mapping not found. Printing IDs that need to be mapped:")
        ids_to_map = hmm_search['Target Accession']
        print(ids_to_map)
        ids_to_map.to_csv("../data/hmm/"+version +
                          "/to_map.csv", header=False, index=False)


def build_map():
    version = "v"+input("hmm version: ")
    hmm_model_hits = pd.read_csv("../data/hmm/"+version+"/hmm_hits.csv")
    pp_map = {}
    for i in range(len(hmm_model_hits)):
        x = hmm_model_hits.iloc[i]
        pp_map[x['Target Accession']] = [x['Start Ali.'], x['End Ali.']]
        # pp_map[x['Target Accession']] = [x['Start Env.'], x['End Env.']]
    return pp_map
