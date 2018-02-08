# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# Choose potentials (spin or spin-lattice)
pair_style hybrid/overlay eam/alloy pair/spin/exchange 4.0 
pair_coeff * * eam/alloy ../examples/SPIN/cobalt/Co_PurjaPun_2012.eam.alloy Co
air_coeff * * pair/spin/exchange exchange 4.0 0.0446928 0.003496 1.4885

# Setup neighbor style
neighbor 0.1 bin
neigh_modify every 10 check yes delay 20

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic

# Setup output
compute out_mag    all compute/spin
compute out_temp   all temp

variable magnorm   equal c_out_mag[5]
variable tmag      equal c_out_mag[7]

thermo_style    custom step time v_magnorm v_tmag temp
thermo          1


#thermo		1
#thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
#thermo_modify norm no
