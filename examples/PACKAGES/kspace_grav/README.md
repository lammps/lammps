# Examples for kspace_style ewald/grav and pair_style grav/long

Newtonian gravity in a fully periodic box, computed with the Ewald method of
Hernquist, Bouchet & Suto, ApJS 75, 231 (1991).  The point masses move in a
uniform negative background of the same total mass, as in comoving
cosmological simulations.  All examples use `units lj` with G = 1.

## in.jeans: Jeans instability of a cold medium

Particles on a lattice are an (unstable) equilibrium.  A plane-wave
displacement along x grows exactly as exp(Gamma t), with
Gamma = sqrt(4 pi G rho), independent of the wavelength.  The thermo
columns `v_mode` and `v_exact` are the measured and exact mode amplitude.

## in.fewbody: a periodic three-body orbit in a periodic box

The orbit "Hill 3 months/yr" of Vanderbei (three unit masses: one body on a
near-circle around a tight binary that turns three "months" per "year") in a box
of side L, which perturbs it through the periodic images and the background.

R. J. Vanderbei, New Orbits for the n-Body Problem, Ann. N.Y. Acad. Sci.
1017, 422-433 (2004), https://doi.org/10.1196/annals.1311.024
