import pandas as pd
import os.path
import requests as r
from Bio import SeqIO
from io import StringIO
from progress.bar import ChargingBar


def getSeq(ID):

    baseUrl = "http://www.uniprot.org/uniprot/"
    currentUrl = baseUrl+ID+".fasta"
    response = r.post(currentUrl)
    cData = ''.join(response.text)

    return list(SeqIO.parse(StringIO(cData), 'fasta'))[0].seq


def downloadSequences(pp_map):
    version = "v"+input("sequences version: ")
    if os.path.isfile("../data/psiblast/"+version+"/sequences.csv"):
        print("Already downloaded sequences")
        df = pd.read_csv("../data/psiblast/"+version+"/sequences.csv")
        return df
    else:
        # here we download all the sequences from the matches
        seqs = []
        bar = ChargingBar('Downloading sequences', max=len(pp_map))
        for accession in pp_map:
            seqs.append([accession, getSeq(accession), pp_map.get(
                accession)[0], pp_map.get(accession)[1]])
            bar.next()
        bar.finish()
        df = pd.DataFrame(data=seqs, columns=[
            'ID', 'Sequence', 'Start', 'Stop'])
        df.to_csv("../data/psiblast/"+version+"/sequences.csv", index=False)
        return df


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
