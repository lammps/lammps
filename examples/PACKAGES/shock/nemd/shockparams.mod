# matches conditions of Holian et al. (2010)
# https://doi.org/10.1063/1.3486088
# J. Chem. Phys. 133, 114502 (2010)

variable 	nequil index 1000               # equilbration run length
variable 	nsteps index 3000               # shock run length
variable 	rho index 0.7071067811865475    # initial density
variable 	temp index 1.251                # initial temperature
variable        nx index 30                     # fcc units cells in axial direction
variable        ny index 5                      # fcc units cells in laterial direction
variable        nz index 5                      # fcc units cells in laterial direction
variable	nthermo equal 100               # human output interval
variable        nevery equal 1                  # computer sampling interval
variable	eps index 1.0e-5                # offset of first atom from (0,0,0)
variable	dt index 0.0005                 # timestep
variable	up index 22.4                   # shock particle velocity
variable	mass index 1.0                  # LJ mass
variable        leftwindow index 50             # number of time samples before shockfront
variable        rightwindow index 150           # number of time samples after shockfront
variable        leftbulkpad index 5.0           # ignore particles within this distance from left surface 
variable        rightbulkpad index 30.0         # ignore particles within this distance from right surface 
