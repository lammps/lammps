This package contains a bias potential that is used to study solid-liquid transitions with the interface pinning method.
An interface between a solid and a liquid is simulated by applying a field that bias the system towards two-phase configurations.
This is done by adding a harmonic potential to the Hamiltonian. The bias field couple to an order-parameter of crystallinity Q:

  U_bias = 0.5*k*(Q-a)^2

Here, We user long-range order for "crystallinity". Q=rho_k wher rho_k is the collective density field.

# References
The main reference for the method is
 [Ulf R. Pedersen, J. Chem. Phys. 139, 104102 (2013)] 

Please visit 
  urp.dk/interface_pinning.htm 
for a detailed bibliography. 

# Build
Remember to include the following command when building LAMMPS
  make yes-user-pinning

# Use

   fix [name] [groupID] rhok [nx] [ny] [nz] [kappa] [anchor-point]

where the parameters set the harmonic bias potential U=0.5*kappa*(|rho_k|-anchor-point)^2 
with the wave-vector elements of rho_k to k_x = (2 pi / L_x) * n_x, k_y = (2 pi / L_y) * n_y and k_z = (2 pi / L_z) * n_z.

# Usage example
In the following we will apply use the interface pinning method for the Lennard-Jones system (trunctaed at 2.5)
at temperature 0.8 and pressure 2.185. This happens to be a coexistence state-point, but we will later show how interface pinning
can be used to determine this. The present directory contains input files, that we will use.

## Density of crystal
First we will determine the density of the crystal with the following LAMMPS input file
{crystal.lmp}
from the output we get that the average density is 0.9731. We need this density to ensure hydrostatic pressure
when in the crystal slap of a two-phase simulation.

## Setup two-phase configuration
Next, setup a two-phase configuration using the density determined in the previous step.
{setup.lmp}


## Setup two-phase configuration
Finally, we run simulation with the bias field applied.
{pinning.lmp}

# Contact
  Ulf R. Pedersen
  http://www.urp.dk
  ulf AT urp.dk

# Cite
Please cite 
  [Ulf R. Pedersen, J. Chem. Phys. 139, 104102 (2013)] 
when using the package for a publication.
