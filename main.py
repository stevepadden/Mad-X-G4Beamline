#Steps
#1 Read in the Twiss tables requested
#2 Convert the tables into objects
#3 Place objects in a text file
import pymadx
import os

#1)
File="twiss_centre_lne00_lne01_lne03_lne04_nom.tfs"#Point towards the TFS file location returned from MAD-X twiss output"
g4bl_fname = "Automatch_No_Spread.g4bl"	#The file name that the python script creates - read into G4BL after
g4bl_profile_fname = "profile_no_spread.txt"	#Name of the beam profile file - used for analysis via the user 
g4bl_location = "/home/steve/G4beamline-3.06/bin" #Locaion of the G4Beamline bin directory - used for sourcing etc within linux
from Twiss_Read import Twiss_read		#Importing the ability to read objects from the file Twiss read - this needs editing and is where each object is defined.
beamstart = "LNE00LNE01LNE03LNE04$START"	#Where in the TFS file the beam construction is started - this is a FLAG in the tfs file	
tbeam = "beam gaussian particle=anti_proton nEvents=300 x=%s y=%s "	#Setting up a basic beam 
ebtoggle = 1 # set to 1 for equads, 0 for magnetic - If using electric, G4BL will need to be adjusted before compile. Without adjustments change Egrad to just Grad in twiss_read.

read = Twiss_read(File,ebtoggle=ebtoggle) #Reading in twiss file to new class object read of type Twiss_Read
quads=read.getquads()   #Extracting the quads from the object called read
hkickers=read.getHkickers() #Getting all horizontal kickers
vkickers=read.getVkickers() #Getting all veritcal kickers
bends=read.getbends()	#Getting all bending objects
beam_start = read.getBeamstart(beamstart)	#Finding the beam start
#print(beam_start)
all = read.getall()		#Gets ALL objects within the transfer line

#The below lines can be used to create a "generic" beam each run, the code currently does not use this but it is included for user ease.
tfs_header = pymadx.Data.Tfs(File)	
print(tfs_header.header)
EX = tfs_header.header["EX"]    #emittance X
EY = tfs_header.header["EY"]
relgamma = tfs_header.header["GAMMA"]  #Relatvisitic gamma
relbeta = tfs_header.header["BETA"]
import numpy as np
start = tfs_header[beamstart]
betx_start = start["BETX"]
bety_start = start["BETY"]
sigx = np.sqrt(EX*betx_start)
sigy = np.sqrt(EY*bety_start)

#The next line is CRITICAL for matching routines between MAD-X and G4Beamline, this ensures that each point is matched correctly to a point from the other program
for i in all:
    print(i)
    if i["name"] == "LNE04$END":	#This is the END of the line, again a FLAG from g4beamline
        tfs_data = pymadx.Data.Tfs(File)
        madx = pymadx.Plot._GetOpticalDataFromTfs(tfs_data)
        madxbx = madx["betx"]
        zspacing = (float(i["S"])) / (len(madxbx) - 1)  # Matching the number of events in our profile to those of madx for loss computing

from Make_Objects import Objects    #Importing the class that handles object making - Needs heavily editing to build each object unique to the transfer line, providing they are modular.
objects = Objects(zspacing=zspacing,profile=g4bl_profile_fname,file=g4bl_fname,ebtoggle=ebtoggle,fringe=True,nevents=10000) #Making a new version of it

objects.Open_File()     #This writes some stuff to the file, mainly the opening exchanges without handling the objects themselves
objects.Place_Objects(all)	#Places all objects as edited in the file Make_Objects.py

os.system("bash rung4bl.sh %s %s " % (g4bl_location, g4bl_fname))  # > /dev/null" # This is how G4Beamline is ran! Edit the shell script included in this repository as required!, include > /dev/null if you want no output from running G4Beamline.

from Analyse import plotanalyse
#plotanalyse(tfsfile=File,profile=g4bl_profile_fname,name=r"Automatic G4Beamline $\beta$ vs MAD-X $\beta$")
plotanalyse(tfsfile=File,profile=g4bl_profile_fname,name=r"Initial G4Beamline $\beta$ vs MAD-X $\beta$ matching") #Prepares some simple plots to check between each program.
