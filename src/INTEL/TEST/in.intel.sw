# bulk Si via Stillinger-Weber

variable        N index on      # Newton Setting
variable        w index 10      # Warmup Timesteps
variable        t index 6200	# Main Run Timesteps
variable        m index 1       # Main Run Timestep Multiplier
variable        n index 0       # Use NUMA Mapping for Multi-Node
variable        p index 0       # Use Power Measurement

variable	x index 2
variable	y index 2
variable	z index 4

variable	xx equal 20*$x
variable	yy equal 20*$y
variable	zz equal 10*$z
variable        rr equal floor($t*$m)
variable        root getenv LMP_ROOT

newton          $N
if "$n > 0"     then "processors * * * grid numa"

units		metal
atom_style	atomic

lattice		diamond 5.431
region		box block 0 ${xx} 0 ${yy} 0 ${zz}
create_box	1 box
create_atoms	1 box

pair_style	sw
pair_coeff	* * ${root}/bench/POTENTIALS/Si.sw Si
mass            1 28.06

velocity	all create 1000.0 376847 loop geom

neighbor	1.0 bin
neigh_modify    delay 5 every 1

fix		1 all nve
thermo		1000
timestep	0.001

if "$p > 0"     then "run_style verlet/power"

if "$w > 0"     then "run $w"
run             ${rr}
