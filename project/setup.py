# setup.py
import pandas as pd

df = pd.read_csv("data/input.csv")
print("Domain sequence:", df['Domain sequence'][0])
