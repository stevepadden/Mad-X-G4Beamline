# Mad-X-G4Beamline
A selection of scripts that convert between MAD-X and G4Beamline automatically, aswell as an MCMC minimisation routine

main runs the automatic construction required by conversion from MAD-X to G4Beamline, Twiss_Read extracts details about each object, 
whilst Make_Objects handles all construction.

It is Make_Objects that the user must heavily edit, this is where the deftinitions of each modular object is defined, make sure this is properly constructed!

Analyse simply has some basic plotting routines.

Beam_Prep constructs a beam that is not changing run to run, it is strongly suggested to use this module before anything else!
It will make your life significantly easier - who doesnt want that!

The shell script "rung4bl.sh" will have to be edited, it is essential for running G4BL from python, which interacts nicely with the shell script but less so
an executable. Just point the script at your G4BL setup, nothing else is needed

Emittance_With_Dispersion calculates emittance profiles for each BPM aswell as emittance and dispersion along the line,
just point it at the folders where your BPM's are and let it work. Edit the momentum of your particles in this file!

Handover_Analyse does some simple analysis for just the virtual BPM automatically constructed at the end of the simulated line

Finally two extra, and quite large codes are included.

Beta_Min_MC presents a novel approach to acceletator physics, by utilising a markov-chain-monte-carlo approach to adjusting quadrupole values.
This file is heavily commented and explains more within the text, should you wish to use G4Beamline in conjunction with MCMC, this file may be of use to you!

Backend_Percentile is used to anlayse the minimisation routine, it will provide extra infomation useful for seeing how the minimisation routine is exploring the space!
