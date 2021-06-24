This directory has 6 scripts that compute the viscosity (eta) of fluid
using 6 different methods. 5 of them are for a Lennard-Jones fluid
and the last one is for SPC/E water model. See the discussion in
Section 6.21 of the manual for an overview of the methods and pointers
to doc pages for the commands which implement them.  Citations for the
various methods can also be found in the manual.

These scripts are provided for illustration purposes.  No guarantee is
made that the systems are fully equilibrated or that the runs are long
enough to generate good statistics and highly accurate results.

-------------

These are the 5 methods for computing viscosity of a LJ fluid.  The first 3 are
non-equilibrium methods; the last 2 are equilibrium methods.

in.wall = move a wall to shear the fluid between two walls
in.nemd = use fix deform and fix nvt/sllod to perform a NEMD shear simulation
in.mp = use fix viscosity and the Muller-Plathe method
in.gk = use the Green-Kubo method
in.einstein = use the Einstein version of Green-Kubo method

All the systems have around 800 atoms.  The NEMD methods run for short
times; the G-K and Einstein systems need to run longer to generate good statistics.

The scripts were all run on a single processor.  They all run in a
minute or so and produce the accompanying log files and profile files
(for velocity or momentum flux).

See the Movies page of the LAMMPS web site
(https://www.lammps.org/movies.html), for animations of the NEMD
scripts, created using the dump image command.

The state point of the LJ fluid is rho* = 0.6, T* = 1.0, and Rcut =
2.5 sigma.  This system should have a shear viscosity of about 1.0.

-------------

Here is how to extract viscosity from the log file output for each
method.

The NEMD methods use the formula eta = - dM / Velocity-gradient, where
dM = momentum flux in the y-direction, and Vel gradient = dVelX / dY =
the change in x-velocity over a distance dY in the y-direction.

(1) in.wall.2d

mom flux = pxy
dVelX = Srate = 2.7
dY = Y box length = 41.99

eta = 0.946 = running average output as last log file column

(2) in.nemd.2d

mom flux = pxy
dVelX = velocity of top box edge = Srate = 2.7
dY = Y box length = 36.51

eta = 1.18 = running average output as last log file column

(3) in.mp.2d

mom flux = dMom in Y / 2 / Area-perp-to-Y / dTime
dMom = -1370.2 from log file, tallied by MP
       factor of 2 since system is periodic and dMom goes 2 ways
       Area for 2d = lx
       dTime = elapsed time in tau for accumulating dMom
dVelX = 4th column of log output, from fix ave/spatial
dY = 1/2 of Y box length

eta = 0.997 = running average output as last log file column

(4) in.gk.2d

eta is computed directly within the script, by performing a time
integration of the formula discussed in Section 6.21 of the manual,
analogous to the formula for thermal conductivity given on the compute
heat/flux doc page - the resulting value prints at the end of the run
and is in the log file

eta = 1.07

(5) in.einstein.2d

eta is computed directly within the script, by performing a time
integration of the formula discussed in Section 6.21 of the manual,
analogous to the formula for thermal conductivity given on the compute
heat/flux doc page - the resulting value prints at the end of the run
and is in the log file

eta = 1.07

-------------

in.cos.1000SPCE is an example script of using cosine periodic perturbation method
to calculate the viscosity of SPC/E water model.

The reciprocal of eta is computed within the script, and printed out as v_invVis
in thermo_style command. The result will converge after hundreds of picoseconds.
Then eta is obtained from the reciprocal of time average of v_invVis.

eta = 0.75 mPa*s

Note that the calculated viscosity by this method decreases with increased acceleration.
It is therefore generally necessary to perform calculation at different accelerations
and extrapolate the viscosity to zero shear.
