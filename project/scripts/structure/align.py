import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from pandas import DataFrame
import os
import pandas as pd
from progress.bar import ChargingBar


tm_align = False  # change to perform all-vs-all pairwise structural alignment

path = "../../data/structure/"

pdb_df = pd.read_csv(path+"pdb_cut_positions.csv")

pdb_ids = list(pdb_df['PDB ID'])


if tm_align:
    if not os.path.isdir("..\\..\\data\\structure\\structural_alignments"):
        os.mkdir(
            "..\\..\\data\\structure\\structural_alignments")
    bar = ChargingBar(
        'All-vs-All pairwise structural alignment w/ TM-align', max=(len(pdb_ids)*len(pdb_ids)))
    # *** all-vs-all pairwise structural alignment using TM-align
    for pdb_i in pdb_ids:
        if not os.path.isdir("..\\..\\data\\structure\\structural_alignments\\{}".format(pdb_i)):
            # ? create dir for pdb id
            os.mkdir(
                "..\\..\\data\\structure\\structural_alignments\\{}".format(pdb_i))
        # ? loop over all other domains
        for pdb_j in pdb_ids:
            cmd = "..\\..\\..\\.binx\\TMalign ..\\..\\data\\structure\\pdb_cut\\pdb{}_cut.ent ..\\..\\data\\structure\\pdb_cut\\pdb{}_cut.ent -o ..\\..\\data\\structure\\structural_alignments\\{}\\{}_{}.sup > ..\\..\\data\\structure\\structural_alignments\\{}\\{}_{}.out".format(
                pdb_i, pdb_j, pdb_i, pdb_i, pdb_j, pdb_i, pdb_i, pdb_j)
            # execute the command
            os.system(cmd)
            bar.next()
    bar.finish()


def extract(data):
    data = [d for d in data if d != "\n"]
    data = data[:3]
    first_row = data[0].split(",")
    rmsd = float(first_row[1].replace("RMSD= ", ""))
    seq_id = float(first_row[2].replace(
        "Seq_ID=n_identical/n_aligned= ", "").replace("\n", ""))
    # TM score
    third_row = data[2]
    cut = third_row.replace("TM-score= ", "")
    tm_score = float(cut[:cut.index("(")])
    return (seq_id, rmsd, tm_score)


# ? results for every metric of interest
results = {}

for pdb_i in pdb_ids:
    results[pdb_i] = {}
    # ? loop over all other domains
    for pdb_j in pdb_ids:
        # if pdb_i != pdb_j:
        my_dir = "..\\..\\data\\structure\\structural_alignments\\{}\\{}_{}.out".format(
            pdb_i, pdb_i, pdb_j)
        f = open(my_dir, "r")
        count = 0
        data = []
        for line in f.readlines():
            if line == "\n":
                count += 1
            if count == 3:
                data.append(line)
        psi, rmsd, tm_score = extract(data)
        results[pdb_i][pdb_j] = {'psi': psi,
                                 'rmsd': rmsd, 'tm_score': tm_score}

# results_flat_psi = [["TM-Align | PSI"]]
# results_flat_rmsd = [["TM-Align | RMSD"]]
# results_flat_tm_score = [["TM-Align | TM-SCORE"]]

results_flat_psi = []
results_flat_rmsd = []
results_flat_tm_score = []
columns = []

# save the ids in the first row
for pdb_i in results:
    # results_flat_psi[0].append(pdb_i)
    # results_flat_rmsd[0].append(pdb_i)
    # results_flat_tm_score[0].append(pdb_i)
    columns.append(pdb_i)

ids = []
for pdb_i in results:
    ids.append(pdb_i)
    # array_psi_i = [pdb_i]
    # array_rmsd_i = [pdb_i]
    # array_tm_score_i = [pdb_i]
    array_psi_i = []
    array_rmsd_i = []
    array_tm_score_i = []
    for pdb_j in results:
        array_psi_i.append(results[pdb_i][pdb_j]['psi'])
        array_rmsd_i.append(results[pdb_i][pdb_j]['rmsd'])
        array_tm_score_i.append(results[pdb_i][pdb_j]['tm_score'])
    results_flat_psi.append(array_psi_i)
    results_flat_rmsd.append(array_rmsd_i)
    results_flat_tm_score.append(array_tm_score_i)

# print to console
print("\nPSI")
for row in results_flat_psi:
    print(row)

print("\nRSMD")
for row in results_flat_rmsd:
    print(row)

print("\nTM-SCORE")
for row in results_flat_rmsd:
    print(row)

# save files

df_psi = pd.DataFrame(data=results_flat_psi, columns=columns, index=ids)
df_psi.to_csv("../../data/structure/results/psi.csv", index=False)

df_rmsd = pd.DataFrame(data=results_flat_rmsd, columns=columns, index=ids)
df_rmsd.to_csv("../../data/structure/results/rmsd.csv", index=False)

df_tm_score = pd.DataFrame(data=results_flat_tm_score,
                           columns=columns, index=ids)
df_tm_score.to_csv("../../data/structure/results/tm_score.csv",
                   index=False)

# create heatmaps

sns.heatmap(df_psi, annot=False)
plt.savefig("../../data/structure/heatmaps/heatmap_psi.png")

plt.show()


sns.heatmap(df_rmsd, annot=False)
plt.savefig("../../data/structure/heatmaps/heatmap_rmsd.png")

plt.show()


sns.heatmap(df_tm_score, annot=False)
plt.savefig("../../data/structure/heatmaps/heatmap_tm_score.png")

plt.show()
