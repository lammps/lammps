# Install/unInstall package classes in LAMMPS

if (test $1 = 1) then

    cp -p pair_reax_c.cpp ..
    cp -p fix_qeq_reax.cpp ..
    cp -p fix_reax_c.cpp ..

    cp -p pair_reax_c.h ..
    cp -p fix_qeq_reax.h ..
    cp -p fix_reax_c.h ..

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

    rm ../pair_reax_c.cpp
    rm ../fix_qeq_reax.cpp
    rm ../fix_reax_c.cpp

    rm ../pair_reax_c.h
    rm ../fix_qeq_reax.h
    rm ../fix_reax_c.h

    rm ../reaxc_allocate.cpp
    rm ../reaxc_basic_comm.cpp
    rm ../reaxc_bond_orders.cpp
    rm ../reaxc_bonds.cpp
    rm ../reaxc_control.cpp
    rm ../reaxc_ffield.cpp
    rm ../reaxc_forces.cpp
    rm ../reaxc_hydrogen_bonds.cpp
    rm ../reaxc_init_md.cpp
    rm ../reaxc_io_tools.cpp
    rm ../reaxc_list.cpp
    rm ../reaxc_lookup.cpp
    rm ../reaxc_multi_body.cpp
    rm ../reaxc_nonbonded.cpp
    rm ../reaxc_reset_tools.cpp
    rm ../reaxc_system_props.cpp
    rm ../reaxc_tool_box.cpp
    rm ../reaxc_torsion_angles.cpp
    rm ../reaxc_traj.cpp
    rm ../reaxc_valence_angles.cpp
    rm ../reaxc_vector.cpp

    rm ../reaxc_allocate.h
    rm ../reaxc_basic_comm.h
    rm ../reaxc_bond_orders.h
    rm ../reaxc_bonds.h
    rm ../reaxc_control.h
    rm ../reaxc_defs.h
    rm ../reaxc_ffield.h
    rm ../reaxc_forces.h
    rm ../reaxc_hydrogen_bonds.h
    rm ../reaxc_init_md.h
    rm ../reaxc_io_tools.h
    rm ../reaxc_list.h
    rm ../reaxc_lookup.h
    rm ../reaxc_multi_body.h
    rm ../reaxc_nonbonded.h
    rm ../reaxc_reset_tools.h
    rm ../reaxc_system_props.h
    rm ../reaxc_tool_box.h
    rm ../reaxc_torsion_angles.h
    rm ../reaxc_traj.h
    rm ../reaxc_types.h
    rm ../reaxc_valence_angles.h
    rm ../reaxc_vector.h

fi
