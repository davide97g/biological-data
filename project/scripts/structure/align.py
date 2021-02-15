import os
import pandas as pd
from progress.bar import ChargingBar


tm_align = False

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
            # print(pdb_i + " - vs - " + pdb_j)
            cmd = "..\\..\\..\\.binx\\TMalign ..\\..\\data\\structure\\pdb_cut\\pdb{}_cut.ent ..\\..\\data\\structure\\pdb_cut\\pdb{}_cut.ent -o ..\\..\\data\\structure\\structural_alignments\\{}\\{}_{}.sup > ..\\..\\data\\structure\\structural_alignments\\{}\\{}_{}.out".format(
                pdb_i, pdb_j, pdb_i, pdb_i, pdb_j, pdb_i, pdb_i, pdb_j)
            # execute the command
            os.system(cmd)
            bar.next()
        # print("---")
    bar.finish()


def extract(data: []):
    data = [d for d in data if d != "\n"]
    data = data[:3]
    first_row = data[0].split(",")
    # print(first_row)
    # al = int(first_row[0].replace("Aligned length= ", ""))
    rmsd = float(first_row[1].replace("RMSD= ", ""))
    seq_id = float(first_row[2].replace(
        "Seq_ID=n_identical/n_aligned= ", "").replace("\n", ""))
    # print("AL", al)
    # print("RMSD", rmsd)
    # print("SEQ ID", seq_id)
    return (seq_id, rmsd)


# *** read scores
show_psi = False
show_rmsd = True

flags = ['psi', 'rmsd']

mapping_psi_rmsd = {}
totals_psi = []
totals_rmsd = []
# full matrix
results = {}

for flag in flags:
    print(flag.upper())
    for pdb_i in pdb_ids:
        results[pdb_i] = {}
        output = pdb_i+"  "
        sum_psi = 0  # total psi over all domains
        sum_rmsd = 0  # total rmsd over all domains
        # ? loop over all other domains
        for pdb_j in pdb_ids:
            # if pdb_i != pdb_j:
            my_dir = "..\\..\\data\\structure\\structural_alignments\\{}\\{}_{}.out".format(
                pdb_i, pdb_i, pdb_j)
            f = open(my_dir, "r")
            # print(f.read())
            count = 0
            data = []
            for line in f.readlines():
                if line == "\n":
                    count += 1
                if count == 3:
                    # data lines
                    data.append(line)
            psi, rmsd = extract(data)
            results[pdb_i][pdb_j] = {'psi': psi, 'rmsd': rmsd}
            sum_psi += psi
            sum_rmsd += rmsd
            if flag == 'psi':
                output += str(psi)+"  "
            else:
                output += str(rmsd)+"  "
        print(output)
        mapping_psi_rmsd[pdb_i] = {
            'total_psi': sum_psi, 'total_rmsd': sum_rmsd}
        totals_psi.append(sum_psi)
        totals_rmsd.append(sum_rmsd)
    print("\n")

print("Flattening...\n")
results_flat_psi = [["TM-Align | PSI"]]
results_flat_rmsd = [["TM-Align | RMSD"]]

# save the ids in the first row
for pdb_i in results:
    results_flat_psi[0].append(pdb_i)
    results_flat_rmsd[0].append(pdb_i)

for pdb_i in results:
    array_psi_i = [pdb_i]
    array_rmsd_i = [pdb_i]
    for pdb_j in results:
        array_psi_i.append(results[pdb_i][pdb_j]['psi'])
        array_rmsd_i.append(results[pdb_i][pdb_j]['rmsd'])
    results_flat_psi.append(array_psi_i)
    results_flat_rmsd.append(array_rmsd_i)

print("PSI")
for row in results_flat_psi:
    print(row)

print("RSMD")
for row in results_flat_rmsd:
    print(row)

df_psi = pd.DataFrame(data=results_flat_psi, columns=results_flat_psi[0])
df_psi.to_csv("../../data/structure/psi.csv", index=False, header=False)

df_rmsd = pd.DataFrame(data=results_flat_rmsd, columns=results_flat_rmsd[0])
df_rmsd.to_csv("../../data/structure/rmsd.csv", index=False, header=False)
