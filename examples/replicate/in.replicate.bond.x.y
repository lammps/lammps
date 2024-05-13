# test of replicating system with periodic bonds in both x and y

dimension       2
atom_style      molecular

read_data       data.bond.x.y

#replicate       3 3 1
replicate       3 3 1 bond/periodic

mass            1 1.0
velocity        all create 0.02 87287 loop geom

pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0

bond_style      harmonic
bond_coeff      1 50.0 1.0

special_bonds   fene

fix             1 all nve

write_data      tmp.data.x.y

dump		1 all image 100 tmp.image.x.y.*.ppm type type &
                adiam 0.2 bond type 0.1 zoom 1.6
dump_modify	1 pad 5

#dump		2 all movie 100 tmp.movie.x.y.mpg type type &
#                adiam 0.2 bond type 0.1 zoom 1.6
#dump_modify	2 pad 5

run             5000
