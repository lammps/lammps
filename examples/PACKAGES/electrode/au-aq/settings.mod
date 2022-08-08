units real
# distribute electrode atoms among all processors:
if "$(extract_setting(world_size) % 2) == 0" then "processors * * 2"

atom_style full
pair_style lj/cut/coul/long 15
bond_style harmonic
angle_style harmonic
kspace_style pppm/electrode 1e-7

read_data "data.au-aq"

group bot type 6
group top type 7

group SPC type 1 2 3
group electrolyte type 1 2 3 4 5

fix nvt electrolyte nvt temp 298.0 298.0 241
fix shake SPC shake 1e-4 20 0 b 1 2 a 1

variable q atom q
variable qz atom q*(z-lz/2)
compute qtop top reduce sum v_q
compute qbot bot reduce sum v_q
compute qztop top reduce sum v_qz
compute qzbot bot reduce sum v_qz
compute ctemp electrolyte temp

