#!/bin/sh

if [ ! -f lualammps.so ]
then \
    echo "Need to compile 'lualammps.so' first for this script to work"
    exit 1
fi

cat > example.lua <<EOF
require("lualammps")

-- Lua does not seem to support %cstring_output_maxsize so this
-- currently does not work and complains about a missing argument.
--
-- osinfo = lualammps.lammps_get_os_info(512)
-- print(osinfo)

lmp = lualammps.lammps_open_no_mpi(0,nil,nil)
ver = lualammps.lammps_version(lmp)

npair_styles = lualammps.lammps_style_count(lmp, "pair")
print("LAMMPS includes ", npair_styles, " pair styles")

-- this also does not work (see above)
-- for i = 0, 10, 1 do
--    res = lualammps.lammps_style_name(lmp, 'pair', i, 128)
--    print("Pair style ",i, ': ', res[1])
-- end

lualammps.lammps_command(lmp, "units real")
lualammps.lammps_command(lmp, "lattice fcc 2.5")
lualammps.lammps_command(lmp, "region box block -5 5 -5 5 -5 5")
lualammps.lammps_command(lmp, "create_box 1 box")
lualammps.lammps_command(lmp, "create_atoms 1 box")

boxlo_p = lualammps.new_double_1d(3)
boxhi_p = lualammps.new_double_1d(3)
xy_p = lualammps.new_double_p()
yz_p = lualammps.new_double_p()
xz_p = lualammps.new_double_p()
pflags_p = lualammps.new_int_1d(3)
boxflag_p = lualammps.new_int_p()

lualammps.lammps_extract_box(lmp, boxlo_p, boxhi_p, xy_p, yz_p, xz_p, pflags_p, boxflag_p)

print('boxlo:    ', lualammps.double_1d_getitem(boxlo_p, 0), ' ', lualammps.double_1d_getitem(boxlo_p, 1), ' ', lualammps.double_1d_getitem(boxlo_p, 2))
print('boxhi:    ', lualammps.double_1d_getitem(boxhi_p, 0), ' ', lualammps.double_1d_getitem(boxhi_p, 1), ' ', lualammps.double_1d_getitem(boxhi_p, 2))
print('xy/yz/xz: ', lualammps.double_p_value(xy_p), ' ', lualammps.double_p_value(yz_p), ' ', lualammps.double_p_value(xz_p))
print('periodicity: ', lualammps.int_1d_getitem(pflags_p, 0), ' ', lualammps.int_1d_getitem(pflags_p, 1), ' ', lualammps.int_1d_getitem(pflags_p, 2))
print('boxflag:  ', lualammps.int_p_value(boxflag_p))
lualammps.delete_double_1d(boxlo_p)
lualammps.delete_double_1d(boxhi_p)
lualammps.delete_int_1d(pflags_p)
lualammps.delete_double_p(xy_p)
lualammps.delete_double_p(yz_p)
lualammps.delete_double_p(xz_p)
lualammps.delete_int_p(boxflag_p)

print("LAMMPS version ", ver)
print("Number of created atoms: ", lualammps.lammps_get_natoms(lmp))
print("Current size of timestep: ", lualammps.double_p_value(lualammps.void_p_to_double_p(lualammps.lammps_extract_global(lmp,"dt"))))
lualammps.lammps_close(lmp)
EOF

lua example.lua
