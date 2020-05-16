variable  newton_pair     index  on
variable  newton_bond     index  on
variable  units           index  metal

atom_style       atomic
atom_modify      map array
neigh_modify     delay 2 every 2 check no
timestep         0.0001
units            ${units}
newton           ${newton_pair} ${newton_bond}

pair_style       zero 8.0
read_data        data.metal
