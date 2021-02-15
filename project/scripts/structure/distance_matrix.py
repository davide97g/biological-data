
import scipy.signal
import matplotlib.pyplot as plt
import matplotlib
from Bio.PDB import PDBList, NeighborSearch
from Bio.PDB.PDBParser import PDBParser
import numpy as np


def get_distance_matrix(residues, seq_sep=6):

    # Calculate the distance matrix
    distances = []
    for residue1 in residues:
        if residue1.id[0] == " ":  # Exclude hetero/water residues
            row = []
            for residue2 in residues:
                if residue2.id[0] == " ":  # Exclude hetero/water residues
                    if abs(residue1.id[1] - residue2.id[1]) >= seq_sep:
                        row.append(residue1["CA"] - residue2["CA"])
                    else:
                        row.append(None)
            distances.append(row)

    return np.array(distances, dtype=float)


#########################################


# Fetch a PDB file to the current dir
pdbl = PDBList()

# Input
pdb_id = '1ucd'

pdbl.retrieve_pdb_file(pdb_id, pdir='../data/structure/',
                       file_format='pdb')  # Will save to pdbXXXX.ent
structure = PDBParser(QUIET=True).get_structure(
    pdb_id, "../data/structure/pdb{}.ent".format(pdb_id))
selected_residues = structure[0]['A']  # select chain A of first model (0)

# Calculate the contact map using the NeighborSearch module (fast)
# https://biopython.org/docs/1.75/api/Bio.PDB.NeighborSearch.html
ns = NeighborSearch(
    [atom for residue in selected_residues for atom in residue.get_atoms()])
# level="R" returns pairs of residues in contact considering all atoms
for residue1, residue2 in ns.search_all(3.5, level="R"):
    if residue1.id[0] == " " and residue2.id[0] == " ":  # Exclude hetero/water residues
        if abs(residue1.id[1] - residue2.id[1]) > 6:  # Sequence separation > 6
            print(residue1.id, residue2.id)


sequence_separation = 12

dist_matrix = get_distance_matrix(selected_residues, sequence_separation)


# Plot distance matrix

# Alter the color map to highlight nan values (the diagonal)
current_cmap = matplotlib.cm.get_cmap()
current_cmap.set_bad(color='white')

fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(dist_matrix)
fig.colorbar(im, fraction=0.03, pad=0.05)
plt.savefig('../data/structure/ca_distances_{}.png'.format(pdb_id),
            bbox_inches='tight')


# Calculate distance matrix of another PDB
# Input
pdb_id = '1ioo'
pdbl.retrieve_pdb_file(pdb_id, pdir='../data/structure/',
                       file_format='pdb')  # Will save to pdbXXXX.ent
structure = PDBParser(QUIET=True).get_structure(
    pdb_id, "../data/structure/pdb{}.ent".format(pdb_id))
selected_residues = structure[0]['A']  # select chain A of first model (0)
dist_matrix2 = get_distance_matrix(selected_residues, sequence_separation)

fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(dist_matrix2)
fig.colorbar(im, fraction=0.03, pad=0.05)
plt.savefig('../data/structure/ca_distances_{}.png'.format(pdb_id),
            bbox_inches='tight')


# Align contact maps / distance matrices (cross-correlation)
# https://stackoverflow.com/questions/59651335/aligning-two-binary-matrices-for-maximum-overlap


# Full distance matrix

# Remove NaN from distance matrix in place (not necessary if sequence separation threshold is set to 0)
# Would be better removing diagonal NaN values
# np.nan_to_num(dist_matrix, copy=False, nan=0.0)
# np.nan_to_num(dist_matrix2, copy=False, nan=0.0)
# contact_matrix = dist_matrix
# contact_matrix2 = dist_matrix2


# Contact maps
# check distance matrix has diagonals, otherwise cross correlation fails
contact_matrix = (dist_matrix[:] < 8).astype(float)
contact_matrix2 = (dist_matrix2[:] < 8).astype(float)


print(contact_matrix)

# Shuffle one matrix
random_matrix = np.copy(contact_matrix)
np.random.shuffle([np.random.shuffle(c) for c in random_matrix])

corr = scipy.signal.correlate2d(
    contact_matrix2, contact_matrix, mode="same")  # cross correlation
corr_random = scipy.signal.correlate2d(
    contact_matrix2, random_matrix, mode="same")  # cross correlation random
# corr = scipy.signal.convolve2d(contact_matrix2, contact_matrix)  # convolution
print(corr)

y, x = np.unravel_index(np.argmax(corr), corr.shape)  # find the match
print(corr.max())
print(x, y)

# Plot cross-correlation
fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(corr)
fig.colorbar(im, fraction=0.03, pad=0.05)
plt.savefig('../data/structure/ca_distances_correlation.png'.format(pdb_id),
            bbox_inches='tight')
plt.close()

# Plot cross-correlation random matrix
fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(corr_random)
fig.colorbar(im, fraction=0.03, pad=0.05)
plt.savefig('../data/structure/ca_distances_correlation_random.png'.format(pdb_id),
            bbox_inches='tight')


# See correlation setting sequence separation to 0, 6 and 12 and compare distance matrix with contact maps
