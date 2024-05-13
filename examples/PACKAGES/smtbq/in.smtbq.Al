# Al2O3 crystal, qeq on, minimizes, then calculates elastic constants

variable T_depart equal 300
variable dt equal 0.0002

#Constante
variable rac3 equal sqrt(3.0)
variable rac1_2 equal sqrt(0.5)
variable rac3_2 equal sqrt(1.5)

#Structure
variable a equal 4.05

variable nx equal 10 
variable ny equal 7
variable nz equal 4

variable bx equal ${a}*${nx}*${rac1_2}
variable by equal ${a}*${ny}*${rac3_2}
variable bz equal ${a}*${nz}*${rac3}
# =======================================================================

units		metal
atom_style	charge
dimension       3
boundary        p p p


lattice		sc 1.0
region          box_vide prism 0 ${bx} 0 ${by} 0 ${bz} 0.0 0.0 0.0
create_box      1 box_vide

# Aluminium atoms z = [111]
lattice custom ${a} a1 ${rac1_2}  0.0      0.0 &
                    a2   0.0   ${rac3_2}   0.0 &
                    a3   0.0      0.0    ${rac3} &
                        basis 0.0 0.0 0.0 &
                        basis 0.5 0.5 0.0 &
                        basis 0.5 0.166666667 0.33333 &
                        basis 0.0 0.666666667 0.33333 &
                        basis 0.0 0.333333333 0.66667 &
                        basis 0.5 0.833333333 0.66667

create_atoms 	1 region box_vide

mass		1 26.98

velocity	all create ${T_depart} 277387

pair_style	smtbq
pair_coeff	* * ffield.smtbq.Al Al

neighbor	0.5 bin
neigh_modify	every 20 delay 0 check yes

timestep        ${dt}

thermo_style    custom step temp press pe ke etotal lx ly lz vol 
thermo_modify	flush yes
thermo		1


#dump            5 all custom 1 box_Al.lammpstrj id type q x y z  

fix		3 all nve
run		10

unfix		3
#thermo		15
fix             1 all box/relax tri 0.0 vmax 0.001
minimize        1.0e-8 1.0e-10 1000 10000

unfix		1
thermo		1
fix             3 all nve
run             10

