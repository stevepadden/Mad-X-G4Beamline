#This file uses the backend saved from minimisation to produce quads, it is useful for checking things without having to run a full minimisation routine!


#Imports
import emcee
import pymadx
import numpy as np
import pandas as pd
import os
from matplotlib.pyplot import plot
import matplotlib.pyplot as plt
import corner
sampler = emcee.backends.HDFBackend("/home/steve/Automatch_No_Spread/Alpha_Emittance_Matched_19th_April.h5", read_only=True)	#The location of the backend file produced via minimisation!
flatchain = sampler.get_chain(flat=True)	#Getting the flat chain
from Twiss_Read import Twiss_read
from Make_Objects import Objects  # Importing the class that handles object making
path = "/home/steve/Automatch/"	#Path
File = path+"New_Lines/twiss_centre_lne00_lne01_lne03_lne04_nom.tfs" #Mad-X file location
read = Twiss_read(File)  # Reading in twiss file to new class object read of type Twiss_Read
quads = read.getquads()  # loading quads
iniquads=[]
for i in quads:
    iniquads.append(i["grad"])

ndim = len(quads)
all = read.getall()  # loading all
tfs_data = pymadx.Data.Tfs(File) #Getting TFS data
madx = pymadx.Plot._GetOpticalDataFromTfs(tfs_data)	#Optical data and subsequent extractions
madxbx = madx["betx"]
madxbx = np.float64(np.sqrt(madxbx))
madxby = madx["bety"]
madxby = np.float64(np.sqrt(madxby))
madxs = madx["s"]
madxbxerror = 0.05*np.average(madxbx)
madxbyerror =0.05*np.average(madxby)
global zspacing
zspacing =1
for i in all:
    if i["name"] == "LNE04$END":  #Adjust to the correct end name point! ~ADJUST HERE!!!!!!!
        zspacing = (float(i["S"])) / (
                    len(madxbx) - 1)  # Matching the number of events in our profile to those of madx for loss computing

#Extracting flat, probability and position values from the sampler.
flat = sampler.get_chain(flat=True)
prob = sampler.get_log_prob()
pos = sampler.get_last_sample()

#The maximum probabilty extraction routine - not recomended.
argmax = np.argmax(prob)
#print(argmax)
quads = flat[argmax]
flat_samples = sampler.get_chain(flat=True)
#The percentile extraction routine - recomended
quads=[]
for i in range(ndim):
    mcmc = np.percentile(flat_samples[:,i],[15,50,84])
    #print(mcmc)
    quads.append(mcmc[1])
print(len(quads))

#From here down, some advanced plotting routines are constructed, you should run these, it will show you alot about how the walkers explore the space.

quadnums = [0,1,2,3]	#Selecting some quadrupoles to analyse, here using the first four.

quadnames  = []
print(quadnames)
vals = []
ij = 0
iniquadvals = []
print(flat.shape)
walkers = []
#print(all)
for i in all:
    if i["type"] == "QUADRUPOLE":
        if ij in quadnums:
            quadnames.append(i["name"])
            vals.append(quads[ij])
            iniquadvals.append(iniquads[ij])
            walkers.append(flat[:,ij])
        ij = ij +1
quadnames = [x.replace("LNE.","") for x in quadnames]

fig,axes = plt.subplots(len(quadnums),figsize=(15,7),sharex=True)
for i in range(len(quadnums)):
    ax = axes[i]
    #ax.plot(flat[:, i], "k", alpha=0.6)
    ax.plot(walkers[i],"k",alpha=0.6)
    ax.axhline(y=vals[i], color='#0505FF', linestyle='-')
    ax.axhline(y=iniquadvals[i],color='#FF1515',linestyle='-')
    ax.set_xlim(0, len(flat))
    string = "%s\n(%.4e)" % (quadnames[i],vals[i])
    ax.set_ylabel(string)
    #ax.set_ylabel(quadnames[i] + "\n" +str(vals[i]))
    ax.yaxis.set_label_coords(-0.05, 0.5)
axes[-1].set_xlabel("step number")
from matplotlib.lines import Line2D
#plt.rcParams.update({'font.size': 8})
custom_lines = [Line2D([0], [0], color='#0050FF', lw=4),
                Line2D([0], [0], color='#FF1515', lw=4),
                Line2D([0], [0], color='#000000',alpha=0.4, lw=4)]

fig.legend(custom_lines,['Minimised Value','Initial Value','Walker Value'],loc='upper left',frameon=False,ncol=len(custom_lines))
fig.suptitle("Trace of walker values for four randomly selected quadrupoles")
#fig.suptitle("Trace of walker values over each step")
fig.show()

from matplotlib.lines import Line2D

custom_lines = [Line2D([0], [0], color='#0050FF', lw=4),
                Line2D([0], [0], color='#FF1515', lw=4),
                Line2D([0], [0], color='#000000', lw=4)]


fig1,ax = plt.subplots(len(quadnums),len(quadnums),figsize=(11.69,8.27))
flat_samples = sampler.get_chain(flat=True)
print(flat_samples.shape)
flat_samples = flat_samples[:,quadnums]
print("FLAT SAMPLES FOR CORNER")
print(flat_samples)
print(flat_samples.shape)
fig = corner.corner(
    #flat_samples,labels=quadnames,bins=20,fig=fig1,quantiles=[0.16, 0.5, 0.84])
    flat_samples, labels=quadnames, bins=20, fig=fig1,show_titles=True,title_fmt='.2e',labelpad=0.2,title_kwargs={'fontsize':10},quantiles=[0.16,0.84])

#labels=quadnames,truths = vals, truth_color = '#0505FF'
corner.overplot_lines(fig,vals,color="#0050FF") #Used values
value2 = iniquadvals
corner.overplot_lines(fig, value2, color="#FF1515") #Initial Values
#corner.overplot_points(fig, value2[None], marker="s", color="#FF1515")
fig.legend(custom_lines,['Minimised Value','Initial Value'],loc='upper right',frameon=False)
fig.suptitle("Posterior probability of the four initial quadrupoles")
fig.set_tight_layout('tight')
#fig.set_fontsize(20)
fig.show()



g4bllocation = "/home/steve/G4beamline-3.06/bin"

def readbeta():
    os.system("bash rung4bl.sh %s %s > /dev/null" % (g4bllocation, new_file)) # > /dev/null"
    g4bl_data = pd.read_csv(profile, skiprows=1, sep=" ")
    g4bl_data = g4bl_data.drop(index=0)
    g4bl_data = g4bl_data.apply(pd.to_numeric)
    betax = np.float64(np.sqrt(g4bl_data["betaX"]/1000))
    betay = np.float64(np.sqrt(g4bl_data["betaY"]/1000))
    return (betax, betay)

def model(quads,zpsacing=zspacing,nevents=200):
    #print("In Model")
    j=0
    for i in all:
        if i["type"] == "QUADRUPOLE":
            i["grad"] = quads[j]
            j = j+1
    objects = Objects(file=new_file, nevents=nevents, zspacing=zspacing,profile=profile)
    objects.Open_File()
    objects.Place_Objects(all)
    #t1error = totalerror()
    bx,by=readbeta()
    return(bx,by)

def matcharraylength(A, B):
    #print(len(B))
    #print(len(A))
    lendiff = len(B) - len(A)
    return np.pad(A, pad_width=(0, lendiff), mode='constant',constant_values=-np.inf)


def lnlike(quads,madxxmeasured,madxxerror,madxymeasured,madxyerror,zspacing=zspacing):
    #print("In Lnlike")
    bx,by = model(quads,zpsacing=zspacing)
    print(bx)
    print(len(bx))
    print(len(madxbx))
    bx = matcharraylength(bx,madxbx)
    by = matcharraylength(by,madxby)
    error1 = -0.5*np.sum(np.power (np.subtract(bx,madxxmeasured) ,2) )
    error2 = -0.5*np.sum(np.power(np.subtract(by, madxymeasured), 2) )
    #error1 = -0.5*np.sum(((bx-madxxmeasured)**2/madxxerror))
    #error2 = -0.5*np.sum(((by-madxymeasured)**2/madxyerror))
    #error2 = 0
    #print(error2)
    #print(error1)
    error = error1 + error2
    print("Current error:%s " %error)
    if abs(error)<10:
        print(quads)
        from Analyse import plotanalyse
        plotanalyse(tfsfile=File,profile=profile)
    if error == np.nan:
        print("NAN error? Setting to massive error.")
        error = -np.inf
    if np.isnan(error):
        print("NAN error? Setting to massive error.")
        error = -np.inf
    return error


zspacing=zspacing
nevents=200
new_file="backend_plot.g4bl"
profile =path+"profile_No_Spread_Fixed.txt"
jj = 0
for i in all:
    if i["type"] == "QUADRUPOLE":
        i["grad"] = quads[jj]
        jj = jj+1
        print(i["grad"])
objects = Objects(file=new_file, nevents=nevents, zspacing=zspacing)
objects.Open_File()
objects.Place_Objects(all)

model(quads)
from Analyse import plotanalyse
plotanalyse(tfsfile=File,profile=profile,name=r"Error minimised $\beta$ matching")


print(iniquads)
print(madxbx)
print(madxbxerror)
print(madxby)
print(madxbyerror)

inierror = lnlike(iniquads,madxbx,madxbxerror,madxby,madxbyerror,zspacing=zspacing)
plotanalyse(tfsfile=File,profile=profile,name=r"Automatic G4Beamline $\beta$ vs MAD-X $\beta$")
error = lnlike(quads,madxbx,madxbxerror,madxby,madxbyerror,zspacing=zspacing)

print("Initial error : %f \n Error after adjustment: %f" % (inierror,error))
fig2,ax2 = plt.subplots()
errorbar = np.sqrt(100)
errs = sampler.get_log_prob()
thiserr = errs.std(axis=1)/errorbar
error_per_iteration = errs.mean(axis=1)
error_per_iteration = -1*error_per_iteration
errs_iterations_axis = np.arange(len(error_per_iteration))


plt.errorbar(errs_iterations_axis,error_per_iteration,thiserr)
plt.fill_between(errs_iterations_axis,error_per_iteration-thiserr,error_per_iteration+thiserr,alpha=0.2)

ax2.set_xlabel("Iteration")
ax2.set_ylabel("Average Error per Iteration")#\n $-\sqrt{(G4BL_{x}-MADX_{x})^2+(G4BL_{y}-MADX_{y})^2}$")
ax2.set_title("Error per Iteration")
fig2.show()


np.savetxt("Adjusted_Quads.csv",quads)	#This file saves the new quad values into a csv file, so you can have an easy reccord of them!
