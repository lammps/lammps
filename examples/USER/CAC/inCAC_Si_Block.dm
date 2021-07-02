units       metal

dimension    3
boundary     s s s

#two CAC element styles currently available
#one has charge and the other does not
#the first number specified is the maximum number
#of nodes an element can have in your input model
#the second is the maximum number of internal 
#degrees of freedom an element can have in the model
#atom_style     cac/charge 8 8
atom_style     cac 8 8

#use tiled style for load balancing non uniform
#element scales
comm_style cac

#read model data file using CAC format
read_data     Si_Model_CAC2aCG.txt
#read_restart   restart.30000



#example pair styles and pair_coeff commands
#in addition to the arguments of the non-CAC version
#of this pair style you may specify the one keyword
#to simplify the quadrature scheme

#pair_style   cac/sw
#pair_coeff     * * Si.sw Si

pair_style   cac/tersoff one
pair_coeff * * Si.tersoff Si
	
timestep     0.001

restart      10000 restart.*

#min style for cac
#min_style cac/cg
#minimize 1.0e-4 1.0e-6 100 1000

#computes the kinetic energy using nodal velocities
compute Ntemp all cac/nodal/temp

#thermo_style custom step c_Ntemp 
thermo_style custom step temp 
thermo	     100

#required to weight the elements and load balance multiresolution models
compute Eweight all cac/quad/count
variable Eweights atom c_Eweight
fix comm all balance 1000 1.00 rcb weight var Eweights

#usual lammps output routine; outputs temp and load balancing 
#info
#fix callertest all ave/time 100 1 100 c_Ntemp f_comm[1] f_comm[2] f_comm[3] f_comm file Node_temps.txt
fix callertest all ave/time 100 1 100 c_Ntemp file atom_temps.txt

#rescales the temperature in the same way the lammps version #does; this one just uses nodal velocities instead
#fix tempre all cac/temp/rescale 5 0.01 300 0.02 1.0
#fix 3 all cac/temp/rescale 10 50.0 50.0 1.0 1.0

#region loading block EDGE EDGE 7750 EDGE EDGE EDGE units box
#group displace region loading

#Verlet integrator for nodal quantities
fix NVE all cac/nve

#a viscous dampener for nodal velocities and atoms
#fix visc all cac/viscous 0.1

#outputs all the nodal positions 
#dump	      1 all cac/nodal/positions 100 CACmesh.txt

#outputs all the nodal kinetic energy; using square of nodal
#velocities
#dump	      2 all cac/kinetic/energy 1000 CACmeshkin.txt

#outputs atom equivalent of the CAC model in xyz format
#dump	      3 notch_region1 cac/xyz 2000 atoms.txt

#can be used to fix nodal forces in the specified group to a #value for now
#fix	loading displace cac/setvelocity 0 1 0
#fix	BCs BD cac/setforce 0 0 0

run 5000
