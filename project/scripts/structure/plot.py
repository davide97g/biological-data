from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
import pandas as pd

df = pd.read_csv("../../data/structure/tm_score.csv")

df = df.drop(columns=[df.columns[0]])

Z = linkage(df.to_numpy())
dn = dendrogram(Z, labels=df.columns,
                above_threshold_color='#bcbddc',
                orientation='left')
plt.savefig('../../data/structure/dendrogram.png')
plt.show()
