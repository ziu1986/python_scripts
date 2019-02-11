import matplotlib.pyplot as plt
import numpy as np; np.random.seed(1)

x = np.random.poisson(size=(160))
y = np.random.poisson(size=(160))

fig, ax = plt.subplots()
ax.set_aspect("equal")
hist, xbins, ybins, im = ax.hist2d(x,y, bins=range(6))

for i in range(len(ybins)-1):
    for j in range(len(xbins)-1):
        ax.text(xbins[j]+0.5,ybins[i]+0.5, hist[i,j], 
                color="w", ha="center", va="center", fontweight="bold")

plt.show(block=False)
