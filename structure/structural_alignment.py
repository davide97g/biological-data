import os

user = input("user: ")
input_align = input("tm-align? ")
tm_align = False
if input_align == "true":
    tm_aling = True
else:
    tm_aling = False
pdb_ids = []

pdb_ids_davide = ['4npo', '3bm7', '1y0h', '1sqe',
                  '1tz0', '1iuj', '3hx9', '1vqs', '1lq9', '5k9f']
pdb_ids_darko = ['4co3', '4ush', '5d4n', '1uku',
                 '2nuh', '4e98', '1nza', '1p1l', '6gdx', '1naq']

if user == "davide":
    pdb_ids = pdb_ids_davide
else:
    pdb_ids = pdb_ids_darko


if tm_align:
    # *** all-vs-all pairwise structural alignment using TM-align
    for pdb_i in pdb_ids:
        # ? create dir for pdb id
        os.mkdir(
            "..\\data\\structure\\{}\\structural_alignments\\{}".format(user, pdb_i))
        # ? loop over all other domains
        for pdb_j in pdb_ids:
            print(pdb_i + " - vs - " + pdb_j)
            cmd = "..\\binx\\TMalign ..\\data\\structure\\{}\\pdb{}_cut.ent ..\\data\\structure\\{}\\pdb{}_cut.ent -o ..\\data\\structure\\{}\\structural_alignments\\{}\\{}_{}.sup > ..\\data\\structure\\{}\\structural_alignments\\{}\\{}_{}.out".format(user,
                                                                                                                                                                                                                                                            pdb_i, user, pdb_j, user, pdb_i, pdb_i, pdb_j, user, pdb_i, pdb_i, pdb_j)
            os.system(cmd)
        print("---")


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
# pdb_ids = pdb_ids[:2]  # only the first for now
show_psi = False
show_rmsd = True

flags = ['psi', 'rmsd']

mapping_psi_rmsd = {}
totals_psi = []
totals_rmsd = []
for flag in flags:

    if flag == 'psi':
        print("4.a\tPaste below a 10 x 10 matrix where cells represent the pairwise sequence identity obtained with the structural alignment (not sequence alignment).")
    else:
        print("4.b\tPaste below a 10 x 10 matrix where cells represent the pairwise RMSD.")

    for pdb_i in pdb_ids:
        output = pdb_i+"  "
        sum_psi = 0  # total psi over all domains
        sum_rmsd = 0  # total rmsd over all domains
        # ? loop over all other domains
        for pdb_j in pdb_ids:
            # if pdb_i != pdb_j:
            my_dir = "..\\data\\structure\\{}\\structural_alignments\\{}\\{}_{}.out".format(user,
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

print("4.c\tWhich is the domain more similar to all other domains looking at the sequence identity?")
for pdb in mapping_psi_rmsd:
    if mapping_psi_rmsd.get(pdb).get("total_psi") == max(totals_psi):
        print(pdb, max(totals_psi))
print("4.d\tWhich is the domain more similar to all other domains looking at the RMSD?")
for pdb in mapping_psi_rmsd:
    if mapping_psi_rmsd.get(pdb).get("total_rmsd") == min(totals_rmsd):
        print(pdb, min(totals_rmsd))
