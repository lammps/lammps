# set boundary in main script because ffield is periodic
units real
# distribute electrode atoms among all processors:
if "$(extract_setting(world_size) % 2) == 0" then "processors * * 2"

atom_style full
pair_style lj/cut/coul/long 12

read_data "data.au-elyt"

group bot type 1
group top type 2

# get electrode charges
variable q atom q
compute qbot bot reduce sum v_q
compute qtop top reduce sum v_q

compute compute_pe all pe
variable vpe equal c_compute_pe
variable charge equal c_qtop
fix fxprint all print 1 "${vpe}, ${charge}" file "out.csv"
