# set boundary in main script because ffield is periodic
units real
# distribute electrode atoms among all processors:
if "$(extract_setting(world_size) % 2) == 0" then "processors * * 2"

atom_style full
pair_style lj/cut/coul/long 16
bond_style harmonic
angle_style harmonic
kspace_style pppm/electrode 1e-7
# kspace_modify in main script because ffield is periodic

read_data "data.graph-il"

# replicate 4 4 1 # test different sys sizes

variable zpos atom "z > 0"
group zpos variable zpos
group ele type 5
group top intersect ele zpos
group bot subtract ele top

group bmi type 1 2 3
group electrolyte type 1 2 3 4

fix nvt electrolyte nvt temp 500.0 500.0 100
fix shake bmi shake 1e-4 20 0 b 1 2 a 1

variable q atom q
compute qtop top reduce sum v_q
compute qbot bot reduce sum v_q
compute ctemp electrolyte temp
