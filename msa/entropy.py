"""
Bio.AlignIO module Wiki
https://biopython.org/wiki/AlignIO
"""

from Bio import AlignIO
import numpy as np
import scipy.stats
from collections import Counter

aa = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
      "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

user = input("user: ")

seqs = []
with open("../data/msa/{}/tcoffee.fasta".format(user)) as f:
    for record in AlignIO.read(f, "fasta"):
        # print(record.seq)
        seqs.append(list(record.seq))  # split seq string

entropies = []
map_index_entropy = {}
seqs = np.array(seqs, dtype="str")
for i, column in enumerate(seqs.T):
    count = Counter(column)  # count AA in column
    non_gap = np.count_nonzero(column != "-")  # count non gap AA
    if non_gap / column.size > 0.7:  # keep columns without gaps (almost)
        # AA probability in column
        probabilities = [count.get(k, 0.0) / column.size for k in aa]
        # Zero entropy = complete conservation
        entropy = scipy.stats.entropy(probabilities, base=20)
        entropies.append(entropy)
        map_index_entropy[i] = entropy
        # print("{} {:.2f} {} {}".format(i, entropy, non_gap, column))

#  python3 parse_msa.py | sort -k2 -f -r
threshold = np.percentile(entropies, 95)
print("5.b\tWhich are the most conserved columns looking at the column entropy?")
print("entropy threshold", threshold, "95% percentile")
for i in map_index_entropy:
    if map_index_entropy.get(i) >= threshold:
        print("column", i, "entropy", map_index_entropy.get(i))
