import pandas as pd


# print(df.head())
# print(len(df))


def build_map():
    version = "v"+input("psi blast version: ")
    df = pd.read_csv("../data/psiblast/"+version +
                     "/psiblast.csv", delimiter="\t")

    pp_map = {}  # here we store the predicted positives and the positions where we found a hit
    for i in range(len(df)):
        x = df.iloc[i]
        accession = x['subject acc.ver'][:-2]
        pp_map[accession] = [x['s. start'], x['s. end']]

    # return the map
    return pp_map


if __name__ == "__main__":
    print(build_map())
