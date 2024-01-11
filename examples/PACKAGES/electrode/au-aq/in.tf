# electrodes with constant potential using finite field
# for z-periodic gold-saline electrochemical cell
# using Thomas-Fermi metallicity model: electrode q and qz will be
# smaller because of more delocalized charge

boundary p p p # ffield uses periodic z-boundary and no slab
include settings.mod # styles, groups, computes and fixes

fix conp bot electrode/conp -1.0 1.805132 couple top 1.0 symm on ffield yes etypes on
fix_modify conp tf 6 1.0 18.1715745
fix_modify conp tf 7 1.0 18.1715745

thermo 50
thermo_style custom step temp c_ctemp epair etotal c_qtop c_qbot c_qztop c_qzbot
run 500
