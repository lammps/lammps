
kspace_style  	pppm 1e-6

neighbor        1.0 bin
neigh_modify    delay 0 every 1 check yes

group    	type1 type 1
compute   	charge1 type1 property/atom q
compute   	q1 type1 reduce ave c_charge1
group    	type2 type 2
compute   	charge2 type2 property/atom q
compute   	q2 type2 reduce ave c_charge2
variable   	qtot equal count(type1)*c_q1+count(type2)*c_q2
variable        nqeq equal f_2

thermo_style  	custom step pe c_q1 c_q2 v_qtot v_nqeq
thermo 		10
thermo_modify   format line "%4d %12.9g %12.8f %12.8f %16.12f %6.0f"

timestep 	0.0001

fix  		1 all nve
