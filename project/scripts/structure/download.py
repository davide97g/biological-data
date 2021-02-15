import urllib.request
import pandas as pd
from progress.bar import ChargingBar


path = "../../data/structure/pdb/"

pdb_df = pd.read_csv(path+"../pdb_cut_positions.csv")

domains = list(pdb_df['PDB ID'])

bar = ChargingBar('Downloading pdb files', max=len(domains))

for domain in domains:
    urllib.request.urlretrieve(
        f"http://files.rcsb.org/download/{domain}.pdb", f"{path}pdb{domain}.ent")
    bar.next()
bar.finish()
