This example demonstrate using a bias potential that can be used to study solid-liquid transitions with the interface pinning method. This is done by adding a harmonic potential to the Hamiltonian that bias the system towards two-phase configurations. 

  U_bias = 0.5*k*(Q-a)^2

The bias field couple to an order-parameter of crystallinity Q 
This implimentation use long-range order: Q=|rho_k|, where rho_k is the collective density field of the wave-vector k.

# Reference
[Ulf R. Pedersen, J. Chem. Phys. 139, 104102 (2013)] 

Please visit 
  urp.dk/interface_pinning.htm 
for a detailed bibliography. 

# Use

   fix [name] [groupID] rhok [nx] [ny] [nz] [spring-constant] [anchor-point]

include a harmonic bias potential U_bias=0.5*k*(|rho_k|-a)^2 to the force calculation.
The elements of the wave-vector rho_k is k_x = (2 pi / L_x) * n_x, k_y = (2 pi / L_y) * n_y and k_z = (2 pi / L_z) * n_z.

# The Interface Pinning method for studying melting transitions
We will use the interface pinning method to study melting of the Lennard-Jones system
at temperature 0.8 and pressure 2.185. This is a coexistence state-point, and the method
can be used to show this. The present directory contains the input files:

  crystal.lmp
  setup.lmp 
  pinning.lmp

1.  First we will determine the density of the crystal with the LAMMPS input file crystal.lmp.
    From the output we get that the average density after equbriliation is 0.9731. 
    We need this density to ensure hydrostatic pressure when in a two-phase simulation.

2.  Next, we setup a two-phase configuration using setup.lmp.

3.  Finally, we run a two-phase simulation with the bias-field applied using pinning.lmp.
    The last coulmn in the output show |rho_k|. We note that after a equbriliation period
	the value fluctuates aroung the anchor point (a) -- showing that this is indeed a coexistence
	state point.

The reference [J. Chem. Phys. 139, 104102 (2013)] gives details on using the method to find coexitence state points,
and the referecee [J. Chem. Phys. 142, 044104 (2015)] show how the crystal growth rate can be computed.
That method have been experienced to be most effective in the slightly super-heated regime above the melting temperature.

# Contact
  Ulf R. Pedersen
  http://www.urp.dk
  ulf AT urp.dk
