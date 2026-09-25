# 2d Vicsek-type collective motion: overdamped active Brownian particles
# with a polar alignment torque toward the orientations of their neighbors

variable        temp equal 0.1       # D_t = temp/gamma_t, D_r = temp/gamma_r
variable        gamma_t equal 10.0
variable        gamma_r equal 1.0
variable        fp equal 10.0        # self-propulsion speed v0 = fp/gamma_t = 1.0
variable        kalign equal 1.0     # alignment torque prefactor, alignment rate = kalign/gamma_r
variable        rcut equal 1.5       # alignment cutoff
variable        seed equal 12345

units           lj
dimension       2
atom_style      hybrid sphere dipole
newton          off

lattice         sq 0.4
region          box block 0 20 0 20 -0.5 0.5
create_box      1 box
create_atoms    1 box
set             group all diameter 1.0 density 1.0
set             group all dipole/random ${seed} 1.0

# WCA potential (purely repulsive)
pair_style      lj/cut 1.122462
pair_coeff      * * 1.0 1.0
pair_modify     shift yes
neigh_modify    every 1 delay 0 check yes

# the alignment cutoff exceeds the pair cutoff: extend the communication cutoff
comm_modify     cutoff $(v_rcut+0.5)

fix             step all brownian/sphere ${temp} ${seed} gamma_t ${gamma_t} gamma_r ${gamma_r}
fix             prop all propel/self dipole ${fp}
fix             align all align/neighbor dipole ${kalign} ${rcut} symmetry polar
fix             plane all enforce2d

# polar and nematic order parameters
compute         mu all property/atom mux muy
compute         pol all reduce ave c_mu[1] c_mu[2]
variable        polar equal sqrt(c_pol[1]^2+c_pol[2]^2)
variable        c2 atom 2.0*c_mu[1]^2-1.0
variable        s2 atom 2.0*c_mu[1]*c_mu[2]
compute         nem all reduce ave v_c2 v_s2
variable        nematic equal sqrt(c_nem[1]^2+c_nem[2]^2)
compute         press all pressure NULL virial

thermo_style    custom step time pe v_polar v_nematic c_press
thermo_modify   norm no
thermo          500
timestep        0.002

run             10000
