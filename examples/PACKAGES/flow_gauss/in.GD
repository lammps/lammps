#LAMMPS input script
#in.GD
#see README for details

###############################################################################
#initialize variables
clear

#frequency for outputting info (timesteps)
variable        dump_rate       equal 50
variable        thermo_rate     equal 10

#equilibration time (timesteps)
variable        equil           equal 1000

#stabilization time (timesteps to reach steady-state)
variable        stabil          equal 1000

#data collection time (timesteps)
variable        run             equal 2000

#length of pipe
variable        L               equal 30

#width of pipe
variable        d               equal 20

#flux (mass/sigma*tau)
variable        J               equal 0.1

#simulation box dimensions
variable        Lx              equal 100
variable        Ly              equal 40

#bulk fluid density
variable        dens            equal 0.8

#lattice spacing for wall atoms
variable        aWall           equal 1.0 #1.7472

#timestep
variable        ts              equal 0.001

#temperature
variable        T               equal 2.0

#thermostat damping constant
variable        tdamp           equal ${ts}*100

units           lj
dimension       2
atom_style      atomic


###############################################################################
#create box

#create lattice with the spacing aWall
variable        rhoWall         equal ${aWall}^(-2)
lattice         sq ${rhoWall}

#modify input dimensions to be multiples of aWall
variable        L1 equal round($L/${aWall})*${aWall}
variable        d1 equal round($d/${aWall})*${aWall}
variable        Ly1 equal round(${Ly}/${aWall})*${aWall}
variable        Lx1 equal round(${Lx}/${aWall})*${aWall}

#create simulation box
variable        lx2 equal ${Lx1}/2
variable        ly2 equal ${Ly1}/2
region          simbox block -${lx2} ${lx2} -${ly2} ${ly2} -0.1 0.1 units box
create_box      2 simbox

#####################################################################
#set up potential

mass            1 1.0           #fluid atoms
mass            2 1.0           #wall atoms

pair_style      lj/cut 2.5
pair_modify     shift yes
pair_coeff      1 1 1.0 1.0 2.5
pair_coeff      1 2 1.0 1.0 1.12246
pair_coeff      2 2 0.0 0.0

neigh_modify  exclude type 2 2

timestep        ${ts}

#####################################################################
#create atoms

#create wall atoms everywhere
create_atoms    2 box

#define region which is "walled off"
variable        dhalf equal ${d1}/2
variable        Lhalf equal ${L1}/2
region          walltop block -${Lhalf} ${Lhalf} ${dhalf} EDGE -0.1 0.1 &
                units box
region          wallbot block -${Lhalf} ${Lhalf} EDGE -${dhalf} -0.1 0.1 &
                units box
region          outsidewall union 2 walltop wallbot side out

#remove wall atoms outside wall region
group           outside region outsidewall
delete_atoms    group outside

#remove wall atoms that aren't on edge of wall region
variable        x1 equal ${Lhalf}-${aWall}
variable        y1 equal ${dhalf}+${aWall}
region          insideTop block -${x1} ${x1} ${y1} EDGE -0.1 0.1 units box
region          insideBot block -${x1} ${x1} EDGE -${y1} -0.1 0.1 units box
region          insideWall union 2 insideTop insideBot
group           insideWall region insideWall
delete_atoms    group insideWall

#define new lattice, to give correct fluid density
#y lattice const must be a multiple of aWall
variable        atrue equal ${dens}^(-1/2)
variable        ay equal round(${atrue}/${aWall})*${aWall}

#choose x lattice const to give correct density
variable        ax equal (${ay}*${dens})^(-1)

#change Lx to be multiple of ax
variable        Lx1 equal round(${Lx}/${ax})*${ax}
variable        lx2 equal ${Lx1}/2
change_box      all x final -${lx2} ${lx2} units box

#define new lattice
lattice         custom ${dens} &
                a1 ${ax} 0.0 0.0 a2 0.0 ${ay} 0.0 a3 0.0 0.0 1.0 &
                basis 0.0 0.0 0.0

#fill in rest of box with bulk particles
variable        delta equal 0.001
variable        Ldelt equal ${Lhalf}+${delta}
variable        dDelt equal ${dhalf}-${delta}
region          left block EDGE -${Ldelt} EDGE EDGE -0.1 0.1 units box
region          right block ${Ldelt} EDGE EDGE EDGE -0.1 0.1 units box
region          pipe block -${Ldelt} ${Ldelt} -${dDelt} ${dDelt} -0.1 0.1 &
                units box

region          bulk union 3 left pipe right
create_atoms    1 region bulk

group           bulk type 1
group           wall type 2

#remove atoms that are too close to wall
delete_atoms    overlap 0.9 bulk wall

neighbor        0.3 bin
neigh_modify    delay 0 every 1 check yes
neigh_modify    exclude group wall wall

velocity        bulk create $T 78915 dist gaussian rot yes mom yes loop geom

#####################################################################
#set up PUT
#see Evans and Morriss, Phys. Rev. Lett. 56(20) 1986, p. 2172

#average number of particles per box, Evans and Morriss used 2.0
variable        NperBox equal 8.0

#calculate box sizes
variable        boxSide equal sqrt(${NperBox}/${dens})
variable        nX equal round(lx/${boxSide})
variable        nY equal round(ly/${boxSide})
variable        dX equal lx/${nX}
variable        dY equal ly/${nY}

#temperature of fluid (excluding wall)
compute         myT bulk temp

#profile-unbiased temperature of fluid
compute         myTp bulk temp/profile 1 1 0 xy ${nX} ${nY}

#thermo setup
thermo          ${thermo_rate}
thermo_style    custom step c_myT c_myTp etotal press

#dump initial configuration
# dump            55 all custom 1 all.init.lammpstrj id type x y z vx vy vz
# dump            56 wall custom 1 wall.init.lammpstrj id type x y z
# dump_modify     55 sort id
# dump_modify     56 sort id
run             0
# undump          55
# undump          56

#####################################################################
#equilibrate without GD

fix             nvt bulk nvt temp $T $T ${tdamp}
fix_modify      nvt temp myTp
fix             2 bulk enforce2d

run             ${equil}

#####################################################################
#initialize the COM velocity and run to achieve steady-state

#calculate velocity to add: V=J/rho_total
variable        Vadd            equal $J*lx*ly/count(bulk)

#first remove any COM velocity, then add back the streaming velocity
velocity        bulk zero linear
velocity        bulk set ${Vadd} 0.0 0.0 units box sum yes mom no

fix             GD bulk flow/gauss 1 0 0 #energy yes
#fix_modify     GD energy yes

run             ${stabil}

#####################################################################
#collect data

#print the applied force and total flux to ensure conservation of Jx
variable        Fapp equal f_GD[1]
compute         vxBulk bulk reduce sum vx
compute         vyBulk bulk reduce sum vy
variable        invVol equal 1.0/(lx*ly)
variable        jx equal c_vxBulk*${invVol}
variable        jy equal c_vyBulk*${invVol}
variable        curr_step equal step
variable        p_Fapp format Fapp %.3f
variable        p_jx format jx %.5g
variable        p_jy format jy %.5g
fix             print_vCOM all print ${dump_rate} &
                "${curr_step} ${p_Fapp} ${p_jx} ${p_jy}" file GD.out screen no &
                title "timestep Fapp Jx Jy"

#compute IK1 pressure profile
#see Todd, Evans, and Davis, Phys. Rev. E 52(2) 1995, p. 1627
#use profile-unbiased temperature to remove the streaming velocity
#from the kinetic part of the pressure
compute         spa bulk stress/atom myTp

#for the pressure profile, use the same grid as the PUT
compute         chunkX bulk chunk/atom bin/1d x lower ${dX} units box

#output pressure profile and other profiles
#the pressure profile is (-1/2V)*(c_spa[1] + c_spa[2]), where
#V is the volume of a slice
fix             profiles bulk ave/chunk 1 1 ${dump_rate} chunkX &
                vx density/mass  c_spa[1] c_spa[2] &
                file x_profiles ave running overwrite

#compute velocity profile across the pipe with a finer grid
variable        dYnew equal ${dY}/10
compute         chunkY bulk chunk/atom bin/1d y center ${dYnew} units box &
                region pipe
fix             velYprof bulk ave/chunk 1 1 ${dump_rate} chunkY &
                vx file Vy_profile ave running overwrite

#full trajectory
# dump          7 bulk custom ${dump_rate} bulk.lammpstrj id type x y z
# dump_modify     7 sort id

run             ${run}
