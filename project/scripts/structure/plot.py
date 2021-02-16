from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
import pandas as pd


def create_dendrogram(df, name):
    Z = linkage(df.to_numpy())
    dendrogram(Z, labels=df.columns,
               above_threshold_color='#bcbddc',
               orientation='left')
    plt.savefig(f'../../data/structure/dendrograms/dendrogram_{name}.png')
    plt.show()


matrices = ['tm_score', 'rmsd', 'psi']

for m in matrices:
    df = pd.read_csv(f"../../data/structure/results/{m}.csv")
    create_dendrogram(df, m)
