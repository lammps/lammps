# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

variable Pu string H
print "potential chosen ${Pu}"
# Choose potential
pair_style meam/ms
print		"we just executed"

pair_coeff      * * library.msmeam ${Pu} Ga4  HGa.meam ${Pu} Ga4
# Setup neighbor style
neighbor 1.0 bin
neigh_modify once no every 1 delay 0 check yes

# Setup minimization style
variable dmax equal 1.0e-2
min_style	     cg
min_modify	     dmax ${dmax} line quadratic
compute eng all pe/atom
compute eatoms all reduce sum c_eng

# Setup output
thermo		100
thermo_style custom step temp etotal  press pxx pyy pzz pxy pxz pyz lx ly lz vol c_eatoms
thermo_modify norm yes
