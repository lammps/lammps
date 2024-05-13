# set boundary in main script because ffield is periodic
units real
# distribute electrode atoms among all processors:
if "$(extract_setting(world_size) % 2) == 0" then "processors * * 2"

atom_style full
pair_style lj/cut/coul/long 14
# kspace_style and _modify in main script to test different approaches

read_data "data.au-vac"

group bot molecule 1
group top molecule 2

# ramping potential difference
variable v equal ramp(0,2)

# get electrode charges
variable q atom q
compute qbot bot reduce sum v_q
compute qtop top reduce sum v_q

# get theoretical charges:
# calculate distance dz between electrodes
compute zbot bot reduce max z
compute ztop top reduce min z
variable dz equal c_ztop-c_zbot

# calculate theoretical capacitance as eps0 * area / dz
variable eps0 equal 55.26349406/10000 # epsilon zero
variable capac equal "v_eps0 * lx * ly / v_dz"

# calculate theoretical charges and deviation of constant potential charges from theory
variable qtheory equal "v_capac * v_v"
variable percdev equal "100 * (c_qtop - v_qtheory) / (v_qtheory + 1e-9)" # avoid divide-by-zero
