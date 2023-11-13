
# distribute electrode atoms among all processors:
if "$(extract_setting(world_size) % 2) == 0" then "processors * * 2"

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

thermo_style custom step pe c_qbot c_qtop

