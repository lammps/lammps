# Install/unInstall package classes in LAMMPS

if (test $1 = 1) then

    cp -p pair_reax_c.cpp ..
    cp -p fix_qeq_reax.cpp ..
    cp -p fix_reax_c.cpp ..
    cp -p fix_reaxc_bonds.cpp ..

    cp -p pair_reax_c.h ..
    cp -p fix_qeq_reax.h ..
    cp -p fix_reax_c.h ..
    cp -p fix_reaxc_bonds.h ..

    cp -p reaxc_allocate.cpp ..
    cp -p reaxc_basic_comm.cpp ..
    cp -p reaxc_bond_orders.cpp ..
    cp -p reaxc_bonds.cpp ..
    cp -p reaxc_control.cpp ..
    cp -p reaxc_ffield.cpp ..
    cp -p reaxc_forces.cpp ..
    cp -p reaxc_hydrogen_bonds.cpp ..
    cp -p reaxc_init_md.cpp ..
    cp -p reaxc_io_tools.cpp ..
    cp -p reaxc_list.cpp ..
    cp -p reaxc_lookup.cpp ..
    cp -p reaxc_multi_body.cpp ..
    cp -p reaxc_nonbonded.cpp ..
    cp -p reaxc_reset_tools.cpp ..
    cp -p reaxc_system_props.cpp ..
    cp -p reaxc_tool_box.cpp ..
    cp -p reaxc_torsion_angles.cpp ..
    cp -p reaxc_traj.cpp ..
    cp -p reaxc_valence_angles.cpp ..
    cp -p reaxc_vector.cpp ..

    cp -p reaxc_allocate.h ..
    cp -p reaxc_basic_comm.h ..
    cp -p reaxc_bond_orders.h ..
    cp -p reaxc_bonds.h ..
    cp -p reaxc_control.h ..
    cp -p reaxc_defs.h ..
    cp -p reaxc_ffield.h ..
    cp -p reaxc_forces.h ..
    cp -p reaxc_hydrogen_bonds.h ..
    cp -p reaxc_init_md.h ..
    cp -p reaxc_io_tools.h ..
    cp -p reaxc_list.h ..
    cp -p reaxc_lookup.h ..
    cp -p reaxc_multi_body.h ..
    cp -p reaxc_nonbonded.h ..
    cp -p reaxc_reset_tools.h ..
    cp -p reaxc_system_props.h ..
    cp -p reaxc_tool_box.h ..
    cp -p reaxc_torsion_angles.h ..
    cp -p reaxc_traj.h ..
    cp -p reaxc_types.h ..
    cp -p reaxc_valence_angles.h ..
    cp -p reaxc_vector.h ..

elif (test $1 = 0) then

    rm -f ../pair_reax_c.cpp
    rm -f ../fix_qeq_reax.cpp
    rm -f ../fix_reax_c.cpp
    rm -f ../fix_reaxc_bonds.cpp 

    rm -f ../pair_reax_c.h
    rm -f ../fix_qeq_reax.h
    rm -f ../fix_reax_c.h
    rm -f ../fix_reaxc_bonds.h

    rm -f ../reaxc_allocate.cpp
    rm -f ../reaxc_basic_comm.cpp
    rm -f ../reaxc_bond_orders.cpp
    rm -f ../reaxc_bonds.cpp
    rm -f ../reaxc_control.cpp
    rm -f ../reaxc_ffield.cpp
    rm -f ../reaxc_forces.cpp
    rm -f ../reaxc_hydrogen_bonds.cpp
    rm -f ../reaxc_init_md.cpp
    rm -f ../reaxc_io_tools.cpp
    rm -f ../reaxc_list.cpp
    rm -f ../reaxc_lookup.cpp
    rm -f ../reaxc_multi_body.cpp
    rm -f ../reaxc_nonbonded.cpp
    rm -f ../reaxc_reset_tools.cpp
    rm -f ../reaxc_system_props.cpp
    rm -f ../reaxc_tool_box.cpp
    rm -f ../reaxc_torsion_angles.cpp
    rm -f ../reaxc_traj.cpp
    rm -f ../reaxc_valence_angles.cpp
    rm -f ../reaxc_vector.cpp

    rm -f ../reaxc_allocate.h
    rm -f ../reaxc_basic_comm.h
    rm -f ../reaxc_bond_orders.h
    rm -f ../reaxc_bonds.h
    rm -f ../reaxc_control.h
    rm -f ../reaxc_defs.h
    rm -f ../reaxc_ffield.h
    rm -f ../reaxc_forces.h
    rm -f ../reaxc_hydrogen_bonds.h
    rm -f ../reaxc_init_md.h
    rm -f ../reaxc_io_tools.h
    rm -f ../reaxc_list.h
    rm -f ../reaxc_lookup.h
    rm -f ../reaxc_multi_body.h
    rm -f ../reaxc_nonbonded.h
    rm -f ../reaxc_reset_tools.h
    rm -f ../reaxc_system_props.h
    rm -f ../reaxc_tool_box.h
    rm -f ../reaxc_torsion_angles.h
    rm -f ../reaxc_traj.h
    rm -f ../reaxc_types.h
    rm -f ../reaxc_valence_angles.h
    rm -f ../reaxc_vector.h

fi
