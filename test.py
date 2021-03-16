import matplotlib.pyplot as plt
import pandas as pd
df = pd.read_csv("exmoset/data/QM9_Data.csv")
data = df[["Dipole Moment","Isotropic Polarizability","Heat Capacity at 298.15K"]].to_numpy()

fig = plt.figure(figsize=(16,16))
ax = plt.subplot(projection="3d")

ax.scatter(data[:,0],data[:,1],data[:,2])
plt.show()
