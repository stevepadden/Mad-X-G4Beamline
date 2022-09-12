#This file generates a gaussian beam profile which is later imported into G4Beamline, it reads all data directly from the MAD-X file so the beam is approximately matched
#Edit beamstart etc to reflect your specific beam, it should produce a file called "matched_beam.G4Beam". 

import numpy as np
from random import gauss
import pymadx

File="/home/steve/Automatch_No_Spread/New_Lines/twiss_centre_lne00_lne01_lne03_lne04_nom.tfs"
madx_beam = "inrays_lne00.madx"
g4beam_beam ="matched_beam.G4Beam"
tfs_header = pymadx.Data.Tfs(File)
beamstart = "LNE00LNE01LNE03LNE04$START"
EX = tfs_header.header["EX"]/1000    #emittance X
EY = tfs_header.header["EY"]/1000
relgamma = tfs_header.header["GAMMA"]  #Relatvisitic gamma
relbeta = tfs_header.header["BETA"]
import numpy as np

start = tfs_header[beamstart]
betx = start["BETX"]
bety= start["BETY"]
sigx = np.sqrt(EX*betx)
sigy = np.sqrt(EY*bety)
dx = start["DX"]*relbeta
dy = start["DY"] *relbeta
dpx = start["DPX"] *relbeta
dpy =start["DPY"] *relbeta
alfx = start["ALFX"]
alfy = start["ALFY"]
PDGid = -2212  # particle ID. Change this to reflect which particle is propegated through the line

print("ALFX %.2e \n BETX %.2e \n dx %.2e \n dpx %.2e \n EX %.2e \n" %(alfx,betx,dx,dpx,EX))
print("ALFY %.2e \n BETY %.2e \n dy %.2e \n dpy %.2e \n EY %.2e \n" %(alfy,bety,dy,dpy,EY))

Npart = 10000 #Number of Particles #set to 1 for dispersion calc


pi= 3.14

ex = EX
ey = EY
dp = 0   #Elena beam dispersion del p /p <- check this

mass_pbar = 938.272  # mass antiproton
beam_energy_mev = 0.1  # beam energy in MeV

beam_mom_mev_c = np.sqrt((mass_pbar + beam_energy_mev) ** 2 - (mass_pbar) ** 2)  # momentum in MeV/c
total_energy = mass_pbar + beam_energy_mev
relative_beta = beam_mom_mev_c / total_energy
dp_for_madx = dp * relative_beta

dispx_inj = dx * relative_beta  # madx output * relative_beta
dispxprime_inj = dpx * relative_beta  # madx output * relative_beta
dispy_inj = 0 * relative_beta  # madx output * relative_beta
dispyprime_inj = 0 * relative_beta  # madx output * relative_beta
print(dispx_inj, dispxprime_inj)
print("beam mom mev/c", beam_mom_mev_c)
print("total energy", total_energy, "mev")
print("beta_rel", relative_beta)
print("ex = ", ex, ", ey = ", ey, ", dp = ", dp, "dp_for_madx = ", dp_for_madx)

filej = open(madx_beam, "w")
fileg = open(g4beam_beam, "w")

print("#BLTrackFile inrays_dp beam", file=fileg)
print("#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight", file=fileg)
print("#mm mm mm MeV/c MeV/c MeV/c ns - - - - -", file=fileg)

zpos = 0.0  # beam starting point in z
t = 0  # 0.375 
ParentID = 0
Weight = 1.0
TrackID = 1
sigma_x, sigma__momentum_x, sigma_y, sigma__momentum_y, se = np.sqrt(ex * betx), np.sqrt(ex / betx), np.sqrt(ey * bety), np.sqrt(ey / bety), dp_for_madx
for i in range(Npart):
    momentum_spread = gauss(0, se)
    xpos = gauss(0, sigma_x) + momentum_spread * (dispx_inj / relative_beta)
    momentum_x = (gauss(0, sigma__momentum_x) - (xpos - momentum_spread * (dispx_inj / relative_beta)) * alfx / betx + momentum_spread * (
                dispxprime_inj / relative_beta)) * beam_mom_mev_c
    ypos = gauss(0, sigma_y) + momentum_spread * (dispy_inj / relative_beta)
    momentum_y = (gauss(0, sigma__momentum_y) - (ypos - momentum_spread * (dispy_inj / relative_beta)) * alfy / bety + momentum_spread * (
                dispyprime_inj / relative_beta)) * beam_mom_mev_c
    momentum_z = momentum_spread * beam_mom_mev_c + beam_mom_mev_c
    # Calculate energy taking into account momentum_spread (momentum spread)
    E = np.sqrt((beam_mom_mev_c * (1 + momentum_spread)) ** 2 + mass_pbar ** 2)
    EventID = i
    print("ptc_start, x= %.15f , px= %.15f, y= %.15f, py= %.15f, pt= %.15f ;" \
          % (xpos, momentum_x, ypos, momentum_y, momentum_spread), file=filej)
    print("%.15f %.15f %.15f %.15f %.15f %.15f %.3f %4.f %1.f %1.f %1.f %1.f" \
          % (xpos * 1000, ypos * 1000, zpos * 1000, momentum_x, momentum_y, momentum_z, t, PDGid, EventID, TrackID, ParentID, Weight), file=fileg)
    # filej.close()

filej.close()
fileg.close()

#import main
#import Emittance_w_dispersion
# myfile.close()

