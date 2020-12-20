The input file in2ddipole.brownian demonstrates how to run a
2d simulation of particles undergoing overdamped brownian motion
in both translational and rotational degrees of freedom. Dipole
updating is turned on, so the package DIPOLE must be built
for this example to work.

The input file in3d_virial_on.brownian demonstrates how to run
a similar simulation but in 3d. In this case, the virial
contribution of the brownian dynamics (the sum
sum_i <r_i dot sqrt{2D_t}W>/(3*volume) where W is
a random variable with mean 0 and variance dt) is
calculated via the fix_modify command. For long
enough times, this will be equal to rho*D_t*gamma_t
(the ideal gas term in equilibrium systems). Note that
no dipole updating is performed.

To confirm rotational diffusion is working correctly,
run the above simulations with dump files on and
measure \sum_i<e_i(0) dot e_i(t)> where e_i is the
dipole vector of particle i, and one should
find that this decays as an exponential with
timescale 1/((d-1)*D_r).

Note that both of the simulations above are not long
enough to get good statistics on e.g. ideal gas
pressure, rotational diffusion, or translational diffusion.
