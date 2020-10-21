#!/bin/sh

if [ ! -f rubylammps.so ]
then \
    echo "Need to compile 'rubylammps.so' first for this script to work"
    exit 1
fi

cat > example.rb <<EOF
require 'rubylammps'

osinfo = Rubylammps.lammps_get_os_info(512)
puts osinfo

lmp = Rubylammps.lammps_open_no_mpi(0,nil,nil)
ver = Rubylammps.lammps_version(lmp)

npair_styles = Rubylammps.lammps_style_count(lmp, 'pair')
puts "LAMMPS includes #{npair_styles} pair styles"
(0..9).each { |i|
   res = Rubylammps.lammps_style_name(lmp, 'pair', i, 128)
   if res[0] == 1
      puts "Pair style #{i}: #{res[1]}\n"
   end
}

Rubylammps.lammps_command(lmp, 'units real')
Rubylammps.lammps_command(lmp, 'lattice fcc 2.5')
Rubylammps.lammps_command(lmp, 'region box block -5 5 -5 5 -5 5')
Rubylammps.lammps_command(lmp, 'create_box 1 box')
Rubylammps.lammps_command(lmp, 'create_atoms 1 box')

boxlo_p = Rubylammps.new_double_1d(3)
boxhi_p = Rubylammps.new_double_1d(3)
xy_p = Rubylammps.new_double_p()
yz_p = Rubylammps.new_double_p()
xz_p = Rubylammps.new_double_p()
pflags_p = Rubylammps.new_int_1d(3)
boxflag_p = Rubylammps.new_int_p()

Rubylammps.lammps_extract_box(lmp, boxlo_p, boxhi_p, xy_p, yz_p, xz_p, pflags_p, boxflag_p)

print('boxlo:    ', Rubylammps.double_1d_getitem(boxlo_p, 0), ' ', Rubylammps.double_1d_getitem(boxlo_p, 1), ' ', Rubylammps.double_1d_getitem(boxlo_p, 2), "\n")
print('boxhi:    ', Rubylammps.double_1d_getitem(boxhi_p, 0), ' ', Rubylammps.double_1d_getitem(boxhi_p, 1), ' ', Rubylammps.double_1d_getitem(boxhi_p, 2), "\n")
print('xy/yz/xz: ', Rubylammps.double_p_value(xy_p), ' ', Rubylammps.double_p_value(yz_p), ' ', Rubylammps.double_p_value(xz_p), "\n")
print('periodicity: ', Rubylammps.int_1d_getitem(pflags_p, 0), ' ', Rubylammps.int_1d_getitem(pflags_p, 1), ' ', Rubylammps.int_1d_getitem(pflags_p, 2), "\n")
print('boxflag:  ', Rubylammps.int_p_value(boxflag_p), "\n")
Rubylammps.delete_double_1d(boxlo_p)
Rubylammps.delete_double_1d(boxhi_p)
Rubylammps.delete_int_1d(pflags_p)
Rubylammps.delete_double_p(xy_p)
Rubylammps.delete_double_p(yz_p)
Rubylammps.delete_double_p(xz_p)
Rubylammps.delete_int_p(boxflag_p)

puts "LAMMPS version #{ver}"
print('Number of created atoms: ', Rubylammps.lammps_get_natoms(lmp), "\n")
print('Current size of timestep: ', Rubylammps.double_p_value(Rubylammps.void_p_to_double_p(Rubylammps.lammps_extract_global(lmp,'dt'))), "\n")
Rubylammps.lammps_close(lmp)
EOF

export RUBYLIB=${PWD}:${RUBYLIB-${PWD}}
ruby example.rb
