This example demonstrate using a bias potential that can be used to study solid-liquid transitions 
with the interface pinning method. This is done by adding a harmonic potential to the Hamiltonian 
that bias the system towards two-phase configurations. 

  U_bias = 0.5*K*(Q-a)^2

The bias field couple to an order-parameter of crystallinity Q. The implementation use long-range order:

  Q=|rho_k|, 

where rho_k is the collective density field of the wave-vector k. 
For future reference we note that the structure factor S(k) is given by the variance of the collective density field: 

  S(k)=|rho_k|^2.

## Reference

It is recommended to get familiar with the interface pinning method by reading:

  [Ulf R. Pedersen, J. Chem. Phys. 139, 104102 (2013)] 

A detailed bibliography is provided at

  <http://urp.dk/interface_pinning.htm>

## Use of rhok fix

For this example we will be using the rhok fix.

   fix [name] [groupID] rhok [nx] [ny] [nz] [K] [a]

This fix include a harmonic bias potential U_bias=0.5*K*(|rho_k|-a)^2 to the force calculation.
The elements of the wave-vector k is given by the nx, ny and nz input: 

  k_x = (2 pi / L_x) * n_x, k_y = (2 pi / L_y) * n_y and k_z = (2 pi / L_z) * n_z. 

We will use a k vector that correspond to a Bragg peak.

# The Interface Pinning method for studying melting transitions of the Lennard-Jones (LJ) system

We will use the interface pinning method to study melting of the LJ system
at temperature 0.8 and pressure 2.185. This is a coexistence state-point, and the method
can be used to show this. The present directory contains the input files that we will use:

  in.crystal
  in.setup
  in.pinning

1. First we will determine the density of the crystal with the LAMMPS input file in.crystal.
  From the output we get that the average density after equilibration is 0.9731. 
  We need this density to ensure hydrostatic pressure when in a two-phase simulation.

2. Next, we setup a two-phase configuration using in.setup.

3. Finally, we run a two-phase simulation with the bias-field applied using in.pinning.
  The last column in the output show |rho_k|. We note that after a equilibration period
  the value fluctuates around the anchor point (a) -- showing that this is indeed a coexistence
  state point.

The reference [J. Chem. Phys. 139, 104102 (2013)] gives details on using the method to find coexistence state points,
and the reference [J. Chem. Phys. 142, 044104 (2015)] show how the crystal growth rate can be computed from fluctuations.
That method have been experienced to be most effective in the slightly super-heated regime above the melting temperature.

# Contact

  Ulf R. Pedersen
  http://www.urp.dk
  ulf AT urp.dk
