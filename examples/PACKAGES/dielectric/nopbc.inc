# Interface
newton off
units       lj
atom_style  dielectric
atom_modify map array
dimension   3
boundary    f f f

variable    method index gmres  # gmres = BEM/GMRES
                                # icc   = BEM/ICC*
                                # dof   = Direct optimization of the functional
                                # none
 
variable    data index data.sphere

read_data   ${data}

group interface type 1
group ions type 2 3

pair_style      lj/cut/coul/cut/dielectric 1.122 20.0
pair_coeff      * * 1.0 1.0
pair_coeff      1 1 0.0 1.0

neigh_modify    one 5000

#compute         ef all efield/atom
#dump            1 all custom 100 all.dump id mol type q x y z #c_ef[1] c_ef[2] c_ef[3]
#
#dump            2 interface custom 100 interface.dump id mol type q x y z #c_ef[1] c_ef[2] c_ef[3]
#dump_modify     1 sort id

fix             1 ions nve

if "${method} == gmres" then &
  "fix             3 interface polarize/bem/gmres 1 1.0e-4" &
elif "${method} == icc"&
  "fix             3 interface polarize/bem/icc 1 1.0e-4" &
elif "${method} == dof" &
  "fix             3 interface polarize/functional 1 1.0e-4" &
else &
  "print 'Unsupported method for polarization' "

thermo          1000
thermo_style    custom step evdwl ecoul elong epair
thermo_modify   flush yes

run             0



