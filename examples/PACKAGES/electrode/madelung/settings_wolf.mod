# set boundary in main script because ffield is periodic
units real
# distribute electrode atoms among all processors:
if "$(extract_setting(world_size) % 2) == 0" then "processors * * 2"

atom_style full
pair_style lj/cut/coul/wolf/gauss 0.05 25
neigh_modify one 10000

read_data "data.au-elyt"

pair_coeff 1 1 0 0 2
pair_coeff 2 2 0 0 2
pair_coeff 3 3 0 0 NULL

group bot type 1
group top type 2

# get electrode charges
variable q atom q
compute qbot bot reduce sum v_q
compute qtop top reduce sum v_q

compute compute_pe all pe
variable vpe equal c_compute_pe
variable charge equal c_qtop
compute press all pressure NULL virial
variable p3 equal c_press[3]
fix fxprint all print 1 "${vpe}, ${charge}, ${p3}" file "out.csv"

dump dump_forces all custom 1 forces.lammpstrj id fx fy fz
