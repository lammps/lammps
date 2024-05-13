# Monte Carlo relaxation of perturbed 2d hex lattice

# set these parameters
# make sure neigh skin > 2*deltamove

variable iter loop 3000            # number of Monte Carlo moves
variable deltaperturb equal 0.2    # max size of initial perturbation per dim
variable deltamove equal 0.1       # max size of MC move in one dimension
variable density equal 1.0         # reduced LJ density of atoms on lattice
variable kT equal 0.05             # effective T in Boltzmann factor
variable seed equal 582783         # RNG seed

# problem setup

units		lj
atom_style	atomic
atom_modify     map array sort 0 0.0

dimension       2

lattice		hex ${density}
region		box block 0 10 0 5 -0.5 0.5

create_box	1 box
create_atoms	1 box
mass		1 1.0

pair_style	lj/cut 2.5
pair_coeff	1 1 1.0 1.0 2.5
pair_modify     shift yes

neighbor	0.3 bin
neigh_modify	delay 0 every 1 check yes

variable        e equal pe

# run 0 to get energy of perfect lattice
# emin = minimum energy

run             0
variable        emin equal $e

# disorder the system
# estart = initial energy

variable        x atom x+v_deltaperturb*random(-1.0,1.0,${seed})
variable        y atom y+v_deltaperturb*random(-1.0,1.0,${seed})

set             group all x v_x
set             group all y v_y

#dump		1 all atom 25 dump.mc

#dump		2 all image 25 image.*.jpg type type &
#		zoom 1.6 adiam 1.0
#dump_modify	2 pad 5

#dump		3 all movie 25 movie.mpg type type &
#		zoom 1.6 adiam 1.0
#dump_modify	3 pad 5

variable        elast equal $e
thermo_style    custom step v_emin v_elast pe

run             0

variable        estart equal $e
variable        elast equal $e

# loop over Monte Carlo moves

variable        naccept equal 0
variable        increment equal v_naccept+1
variable        irandom equal floor(atoms*random(0.0,1.0,${seed})+1)
variable        rn equal random(0.0,1.0,${seed})
variable        boltzfactor equal "exp(atoms*(v_elast - v_e) / v_kT)"
variable        xnew equal x[v_i]+v_deltamove*random(-1.0,1.0,${seed})
variable        ynew equal y[v_i]+v_deltamove*random(-1.0,1.0,${seed})
variable        xi equal x[v_i]
variable        yi equal y[v_i]

label           loop

variable        i equal ${irandom}

variable        x0 equal ${xi}
variable        y0 equal ${yi}

set             atom $i x ${xnew}
set             atom $i y ${ynew}

run             1 pre no post no

if              "$e <= ${elast}" then &
                  "variable elast equal $e" &
                  "variable naccept equal ${increment}" &
                elif "${rn} <= ${boltzfactor}" &
                  "variable elast equal $e" &
                  "variable naccept equal ${increment}" &
                else &
                  "set atom $i x ${x0}" &
                  "set atom $i y ${y0}"

next            iter
jump            SELF loop

# final energy and stats

variable       nb equal nbuild
variable       nbuild equal ${nb}

run             0

print           "MC stats:"
print           "  starting energy = ${estart}"
print           "  final energy = $e"
print           "  minimum energy of perfect lattice = ${emin}"
print           "  accepted MC moves = ${naccept}"
print           "  neighbor list rebuilds = ${nbuild}"
