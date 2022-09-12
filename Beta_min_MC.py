#This file uses emcee to minimise the beta functions between G4Beamline and MAD-X, the logic is not explained within this script, if the user wishes for a more indepth analysis behind the logic, you are encouraged to read my thesis entitled "End to End simulations of Antiproton Transport and Degredagtion" - no knowledge of antimatter physics is required to understand the impacts to this script and minimisation functions.

#Included in this file are "total error" and "calc error", These are left in as a result of legacy code, you can ignore them.

#Note any errors due to array differences are likely due to line 36 not handling the end coord of the line correctly

#imports
import emcee 
import pymadx
import numpy as np
import pandas as pd
import os
from matplotlib.pyplot import plot
#load in initial values from madx - adjust all these baths to reflect your unique setup.

path = "/home/steve/Automatch_No_Spread/"	#folder
initial_file = path + "Automatch_No_Spread.g4bl"	#file name of G4BL file that has already been automatically constructed
new_file = path + "automatched_fixed_emit.g4bl"		#The name of the new G4BL file that will be made
g4bllocation = "/home/steve/G4beamline-3.06/bin"	#G4Beamline bin location 
File = path+"New_Lines/twiss_centre_lne00_lne01_lne03_lne04_nom.tfs"	#MAD-X TFS file used
profile =path+"profile_No_Spread_Fixed.txt"	#filename for the profile made by this matching routine
iniprofile = "profile_no_spread.txt"	#Initial profile name
tfs_data = pymadx.Data.Tfs(File)# Extracting data from the MAD-X file
madx = pymadx.Plot._GetOpticalDataFromTfs(tfs_data) #Grabbing optical data
madxbx = madx["betx"]
madxbx = np.float64(np.sqrt(madxbx))
madxby = madx["bety"]
madxby = np.float64(np.sqrt(madxby))
madxs = madx["s"]
nevents=300
zspacing = 30

#load in initial G4BL - not running them yet
from Twiss_Read import Twiss_read #Class that handled object reading
from Make_Objects import Objects  # Importing the class that handles object making

read = Twiss_read(File)  # Reading in twiss file to new class object read of type Twiss_Read
quads = read.getquads()  # loading quads
all = read.getall()  # loading all
#print(all)
#load quadstrenghts into own array
for i in all:
    if i["name"] == "LNE04$END":
        zspacing = (float(i["S"])) / (
                    len(madxbx) - 1)  # Matching the number of events in our profile to those of madx for loss computing


quadstrengths = []

for i in quads:
    quadstrengths.append(i["grad"])
numquads = len(quadstrengths)
#Here we set the maximum upper and lower limits on quad strengths, smaller limits can dramatically speed up convergence times, but may not be optimised, similarly too large a quad strength may be unphysical for the quadrupole set up - be careful here
maximum = 2*(max(quadstrengths))      
minumum = 2*(min(quadstrengths))
print("MIN QUADS = ")
print(min(quadstrengths))
print("MINIMUM USED =")
print(minumum)
print("MAX QUADS = ")
print(max(quadstrengths))
print("MAX USED =")
print(maximum)


def matcharraylength(A, B): #Matches the array length
    lendiff = len(B) - len(A)
    return np.pad(A, pad_width=(0, lendiff), mode='constant',constant_values=-np.inf)


def calcerror(madxBx, g4blBx, madxBy, g4blBy): #An error calculation  - adjust to reflect your specific error calcs!
    g4blBx = matcharraylength(g4blBx, madxbx)
    g4blBx = np.nan_to_num(g4blBx)
    g4blBy = matcharraylength(g4blBy, madxby)
    g4blBy = np.nan_to_num(g4blBy)
    total1 = np.sum(np.power((np.subtract(g4blBx, madxBx)), 2))
    total2 = np.sum(np.power((np.subtract(g4blBy, madxBy)), 2))
    total = total1 + total2
    return (total)

def totalerror():	#Semi-middleman class, enables fast switching between which error you would like to use
    bx, by = readbeta()	#Reads in the beta values
    total = calcerror(madxbx, bx, madxby, by)	#Calls the error calculation
    return total



def readbeta():	#This is called "readbeta" but in reality it is much more than that - it runs G4Beamline for this specific iteration, then reads all the beta values. If you were minimising on, say, beam dimenstions, this is where you would have to change things significantly!
    os.system("bash rung4bl.sh %s %s > /dev/null" % (g4bllocation, new_file)) # > /dev/null"
    g4bl_data = pd.read_csv(profile, skiprows=1, sep=" ")
    g4bl_data = g4bl_data.drop(index=0)
    g4bl_data = g4bl_data.apply(pd.to_numeric)
    betax = np.float64(np.sqrt(g4bl_data["betaX"]/1000))
    betay = np.float64(np.sqrt(g4bl_data["betaY"]/1000))
    return (betax, betay)


def model(quads):	#This models the quads, it is used because we have to essentially change quadrupole values each iteration!
    #print("In Model")
    j=0
    for i in all:
        if i["type"] == "QUADRUPOLE":
            i["grad"] = quads[j]
            j = j+1
    objects = Objects(file=new_file, nevents=nevents, zspacing=zspacing,profile=profile)
    objects.Open_File()
    objects.Place_Objects(all)
    bx,by=readbeta()
    return(bx,by)


def lnlike(quads,madxxmeasured,madxxerror,madxymeasured,madxyerror):	#This class handles the "error" as used within EMCEE itself, the line "error" is what you again change to reflect your error calculations! For example you could minimised spot size at handover at this point
    bx,by = model(quads)	#Extracting Beta x and Beta y , adjust to extract your own minimisation parameters
    bx = matcharraylength(bx,madxbx)	#Match the array lengths between Mad-X and G4Beamline - essential.
    by = matcharraylength(by,madxby)

    error = -0.5*np.sqrt(np.sum((bx**2 + by**2) ))	#The actual error function, here I minimise the entire line, note the error is negative, this is always required by EMCEE, ensure you do the same. 
    print("Current error:%s " %error)	#Error for this iteration, used for monitoring by the user
    if abs(error)<0:
        print(quads)
        from Analyse import plotanalyse
        plotanalyse(tfsfile=File,profile=profile)
    if error == np.nan:
        print("NAN error? Setting to massive error.")
        error = -np.inf	#If the G4Beamline fails to converge (IE the beam hits the wall for example and no data is returned) we set the error to be huge, this ensures "failed" runs are heavily penalised!
    if np.isnan(error):
        print("NAN error? Setting to massive error.")
        error = -np.inf
    return error


def lnprior(quads):	#Prior setup, if the walker explores a space that is outside the pre set boundary conditions, it doesnt run anything and redraws.

    t=False
    for i in quads:
        if i > maximum or i< minumum:
            t=True
    if t == True:
        return -np.inf
    else:
        return 0.0


madxbxerror = 0.05*np.average(madxbx)	#Generic 5 % error on MAD-X values, adjust to your desired error.
madxbyerror =0.05*np.average(madxby)


def lnprob(quads,madxxmeasured,madxxerror,madxymeasured,madxbyerror):	#Brings together the Lnprior (initial quad limits) and actual error calculation (LnLike) into one function over which EMCEE will minimise.

    lp = lnprior(quads)
    if not np.isfinite(lp):
        return -np.inf
    xx= lp+lnlike(quads,madxxmeasured,madxxerror,madxymeasured,madxbyerror)
    if xx == np.nan:
        xx = -np.inf
    return xx

ndim,nwalkers = len(quads),100	#Definining the number of walkers, here I use 100 walkers per quadrupole in the line
niter = 75	#Total number of iterations used for training the model
burnin = 12	#Number of burn in iterations, used for exploring space to improve convergence times
initial = np.array(quadstrengths)	#The initial quadrupole strengths
inierror = lnlike(quadstrengths,madxbx,madxbxerror,madxby,madxbyerror)	#The initial error
from Analyse import plotanalyse	#Some initial analysis
plotanalyse(tfsfile=File,profile=iniprofile,name="initial")	#The plot of initial matching between G4Beamline and MAD-X
p0 = [initial+2e-4*np.random.randn(ndim) for i in range(nwalkers)]`#The starting values of our quads used within G4Beamline, initiated in some small gaussian distribution around their actual value, ensures that if the beam is already perfectly optimised, that we explore the space surrounding it

def main(p0,nwalkers,niter,ndim,lnprob):	#The meat of this code, this runs the EMCEE minimisation routines! 
    backend.reset(nwalkers, ndim)	#Clearing the backend
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(madxbx,madxbxerror,madxby,madxbyerror),backend=backend)	#Spinning up the sampler, initated with specific errors and walker dimensions etc, also saves to the backend 

    print ("Burning in")	
    p0, _, _ = sampler.run_mcmc(p0,burnin)	#running burn in
    res=plot(sampler.chain[:,:,0].T,"-",color="k",alpha=0.3)	#Viewing residuals
    sampler.reset()	#Clearing out the "burn in " values
    print("Running Production")	#Performing minimisation
    backend.reset(nwalkers, ndim)	#resetting backend before minimisation
    pos,prob,state = sampler.run_mcmc(p0,niter)	#minimising! this runs for some time (days!). Good things come to those who wait! Suggest using low iteration numbers for initial testing, results will be rubbish but you can check the code actually works!


    return sampler,pos,prob,state	#Returning all the sampler data!


filename = "Alpha_Emittance_Matched.h5" #Name of the .h5 file used to store all minimisation data, huge file but strongly reccomended as if minimisation breaks you can resume from this!
backend = emcee.backends.HDFBackend(filename)	


sampler,pos,prob,state = main(p0,nwalkers,niter,ndim,lnprob)	#Setting up the sampler by calling the function "Main" 

print(sampler)	#Printing the sampler
samples = sampler.flatchain	#extrcting the sampled chain!


#Printing some things returned from the sampler, this is not necessary but is good to see to ensure its working as intended
print("INCOMING FLATCHAIN")
flat = sampler.get_chain(flat=True)
print(flat)
print("INCOMING POS")
print(pos)
print("INCOMING PROB")
print(prob)

#READ THE NEXT LINES!
#quad vals 1!
#Key lines! We use the maximum value from the probability (ie the error),
#Set the quads equal to the pos index value for that run.

#argmax = np.argmax(prob)
#quads = pos[argmax]
#print(quads)


#Quad vals 2 - preffered!

#There is an alternative "quads" value which i reccomend, rather than the "perfect" value for that quad, use the 50th percentile (ie middle of the gaussian) the sampler explored, it is more stable and minimises better!
quads = []
for i in range(ndim):
	mcmc = np.percentile(flat[:,i],[15,50,84])	#Extracting the 15th, 50th and 84th percentiles (middle gaussian and 1 sigma either side)
	#print(mcmc)
	quads.append(mcmc[1])	#Selecting the 50th percentile as the optimised quad values

#Constructs a simulation which has the optimised quadrupole values
jj = 0
for i in all:
    if i["type"] == "QUADRUPOLE":
        i["grad"] = quads[jj]
        jj = jj+1
        print(i["grad"])
objects = Objects(file=new_file, nevents=nevents, zspacing=zspacing)
objects.Open_File()
objects.Place_Objects(all)
err = totalerror()
print(err)

error = lnlike(quads,madxbx,madxbxerror,madxby,madxbyerror)	#Finding the error of the optimised G4Beamline code

from Analyse import plotanalyse
plotanalyse(tfsfile=File,profile=profile,name="Adjusted")	#Analysing the optimised (adjusted) G4Beamline code

print("Initial error : %f \n Error after adjustment: %f" % (inierror,error))	#Comparring to initial errors!





import matplotlib.pyplot as plt
quadnums = [0,6,19,22]	#Selecting some quadrupoles to look at for our analysis! Pick and random numbers here!
#quadnums = [0,1,2,3,4,5]
quadnames  = []
vals = []
ij = 0
s=[]
inivals = []

for i in all:
    if i["type"] == "QUADRUPOLE":
        if ij in quadnums:
            quadnames.append(i["name"].replace("LNE.",""))
            vals.append(quads[ij])
            inivals.append(initial[ij])
        ij = ij +1


#Everything beneath here is plotting routines to analyse how the beta minimisation has worked.

fig,axes = plt.subplots(len(quadnums),figsize=(15,7),sharex=True)
for i in range(len(quadnums)):
    ax = axes[i]
    ax.plot(flat[:, i], "k", alpha=0.8)
    ax.axhline(y=vals[i], color='r', linestyle='-')
    ax.axhline(y=inivals[i],color='royalblue',linestyle='-')
    ax.set_xlim(0, len(flat))
    string = "%s\n(%.3e)" % (quadnames[i],vals[i])
    ax.set_ylabel(string,fontsize=9)
    #ax.set_ylabel(quadnames[i] + "\n" +str(vals[i]))
    ax.yaxis.set_label_coords(-0.05, 0.5)
    ax.tick_params(axis="y",labelsize=7)
axes[-1].set_xlabel("step number")
fig.show()

#Get correlation time & thin by roughly half of that
#discard is number of iterations of N walkers to throw away
#Whilst thin takes only every "thin" number of steps, ie thin=8 == every 8 steps
# tau = sampler.get_autocorr_time()
# print(tau)
# flat_samples = sampler.get_chain(discard=10, thin=8, flat=True)
# print(flat_samples.shape)
import corner
fig1,ax = plt.subplots(len(quadnums),len(quadnums),figsize=(15,15))
flat_samples = sampler.get_chain(flat=True)
print(flat_samples.shape)
flat_samples = flat_samples[:,quadnums]
print(flat_samples.shape)
fig = corner.corner(
    flat_samples, labels=quadnames, truths=vals,bins=15,fig=fig1,quantiles=[0.16, 0.5, 0.84])
fig.show()

fig2,ax2 = plt.subplots()
errorbar = np.sqrt(nwalkers)
errs = sampler.get_log_prob()
error_per_iteration = errs.mean(axis=1)
error_per_iteration = -1*error_per_iteration
errs_iterations_axis = np.arange(len(error_per_iteration))
#plt.scatter(errs_iterations_axis,error_per_iteration,marker="+")
plt.errorbar(errs_iterations_axis,error_per_iteration,errorbar)
plt.fill_between(errs_iterations_axis,error_per_iteration-errorbar,error_per_iteration+errorbar,alpha=0.5)
ax2.set_xlabel("Iteration")
ax2.set_ylabel("Error per Iteration")#\n $-\sqrt{(G4BL_{x}-MADX_{x})^2+(G4BL_{y}-MADX_{y})^2}$")
ax2.set_title("Error per Iteration")
fig2.show()

fig3,ax3 = plt.subplots()
tt = sampler.get_log_prob(flat=True)
it = np.arange(len(tt))
aa = np.argmax(tt)
print(flat[aa]) #potential improvement?
plt.scatter(it,tt)

