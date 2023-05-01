import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("energy_energyEdepSpectrum.txt",sep=" ",header=9)

df = df.rename(columns={'0.0005': 'edep', '0.001': 'bin', '1588.3': 'freq'})


#print(df)

df = df[df["edep"] < 0.9]

print(df)
df.plot(x='edep',y='freq',kind='line')

#plt.xlim = [0, 0.8]
plt.show()