import numpy as np
import numpy.random
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm as CM



f = open("results.txt")
lines = f.readlines()
f.close()

x = []
y = []

for l in lines:
    spl = l.split()
    x.append(float(spl[0]))
    y.append(float(spl[1]))

x = np.array(x)
y = np.array(y)

heatmap, xedges, yedges = np.histogram2d(x, y, bins=30)
#extent = [xedges[0], xedges[-1], yedges[-1], yedges[0]]
extent = [0, 40, 0, 6]

print np.average(x)
print np.std(x, ddof=1)

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.imshow(heatmap, extent=extent, cmap=CM.jet, origin='lower', interpolation='nearest')
cbar = fig.colorbar(cax)
#ax.scatter(x, y)

ax.set_aspect('auto')

plt.show()
