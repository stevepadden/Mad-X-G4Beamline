import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sea
dir = "/home/steve/Automatch_No_Spread"
filename = os.path.join(dir, "Handover.txt")  # Setting the file to this iteration
filename = filename.strip(
    ".txt")  # Removing the .txt for later formatting - if your file isnt a .txt change to the extension here.
pd.set_option('display.max_columns', None)  # set display columns pandas
cols = ["x", "y", "z", "Px", "Py", "Pz", "t", "PDGid", "EventID", "TrackID", "ParentID", "Weight"]
data = pd.read_csv(filename+".txt", sep=" ", skiprows=3, header=None, names=cols)  # read data with above vols
print(data)
title = "handover"

def findmax(x, y):
    x1 = max(x)
    x2 = abs(min(x))
    y1 = max(y)
    y2 = abs(min(y))
    if x2 > x1:
        maxx = x2
        minx = -x2
    else:
        maxx = x1
        minx = -x1

    if y2 > y1:
        maxy = y2
        miny = -y2
    else:
        maxy = y1
        miny = -y1
    endlim = (max(maxx, maxy))
    return minx, maxx, miny, maxy, endlim

nreplicas= 100
bsfrac = 0.3

rmsx_a = []
rmsy_a = []
rmspx_a = []
rmspy_a = []

for i in range(nreplicas):
    dat = data.sample(frac=bsfrac,replace=True)
    rmsx = np.sqrt(np.sum(np.power(dat["x"], 2)) / len(dat["x"]))
    rmsy = np.sqrt(np.sum(np.power(dat["y"], 2)) / len(dat["y"]))
    rmspx = np.sqrt(np.sum(np.power(dat["Px"], 2)) / len(dat["Px"]))
    rmspy = np.sqrt(np.sum(np.power(dat["Py"], 2)) / len(dat["Py"]))

    rmsx_a.append(rmsx)
    rmsy_a.append(rmsy)
    rmspx_a.append(rmspx)
    rmspy_a.append(rmspy)

rmsx_a = np.array(rmsx_a)
rmsy_a = np.array(rmsy_a)
rmspx_a = np.array(rmspx_a)
rmspy_a = np.array(rmspy_a)

N = len(rmsx_a)
rxa = rmsx_a.mean()
rxe = rmsx_a.std()/np.sqrt(N)
rya = rmsy_a.mean()
rye = rmsy_a.std()/np.sqrt(N)

rpxa = rmspx_a.mean()
rpxe = rmspx_a.std()/np.sqrt(N)
rpya = rmspy_a.mean()
rpye = rmspy_a.std()/np.sqrt(N)


print("RMSX : %.2f pm %.4f : RMSY : %.2f pm %.4f" %(rxa,rxe,rya,rye))
print("RMS PX : %.6f pm %.2e : RMS PY : %.6fe pm %.2e " % (rpxa,rpxe,rpya,rpye))

print(rmspx_a)






# PLOTTING X Y HEATMAP
cmap = "viridis"
fig, ax = plt.subplots()
heatmap, xedges, yedges = np.histogram2d(data["x"], data["y"], bins=50)
rmsx = np.sqrt(np.sum(np.power(data["x"], 2)) / len(data["x"]))
print("RMSX")
print(rmsx)
rmsy = np.sqrt(np.sum(np.power(data["y"], 2)) / len(data["y"]))
print("RMSY")
print(rmsy)
# print(heatmap)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
# sea.heatmap(heatmap,square=True)
# fig.show()
# plt.clf()
ax.set_title("Spatial profile at "+ title)
ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")
minx, maxx, miny, maxy, endlim = findmax(data["x"], data["y"])
# plt.xlim((minx,maxx))
# plt.ylim((miny,maxy))
ax.set_xlim((-endlim, endlim))
ax.set_ylim((-endlim, endlim))
tt = ax.imshow(heatmap.T, extent=extent, origin='lower', cmap=cmap, aspect="equal")
ax.set_facecolor('#440154')
str = "$\sigma_{RMSX} =$ %.2f (mm) \n$\sigma_{RMSY} =$ %.2f (mm)" % (rmsx, rmsy)
ax.plot([], [], ' ', label=str)
ax.legend(frameon=False, facecolor='#440154', labelcolor="w")
cbar = plt.colorbar(tt, ax=ax)

fig2, ax2 = plt.subplots()
heatmap, xedges, yedges = np.histogram2d(data["Px"], data["Py"], bins=50)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
rmspx = np.sqrt(np.sum(np.power(data["Px"], 2)) / len(data["Px"]))
rmspy = np.sqrt(np.sum(np.power(data["Py"], 2)) / len(data["Py"]))

# endlim = 10
# plt.clf()
ax2.set_title("Momentum profile at " + title)
ax2.set_xlabel("Px (MeV/c)")
ax2.set_ylabel("Py (MeV/c)")
minx, max2x, miny, max2y, endlim = findmax(data["Px"], data["Py"])
str =  "$\sigma_{RMS PX} =$ %.5f (MeV/c) \n$\sigma_{RMS PY} =$ %.5f (MeV/c)" % (rmspx, rmspy)
ax2.plot([], [], ' ', label=str)
ax2.legend(frameon=False, facecolor='#440154', labelcolor="w")
# plt.xlim((minx,max2x))
# plt.ylim((miny,max2y))
# endlim = 0.1
ax2.set_xlim((-endlim, endlim))
ax2.set_ylim((-endlim, endlim))
tt = ax2.imshow(heatmap.T, extent=extent, origin='lower', cmap=cmap, aspect="equal")
ax2.set_facecolor('#440154')
fig2.colorbar(tt, ax=ax2)

plt.show()