import lammps

infile = "in.grid.pair"
faclist = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]

print("# Timestep DeltaE DeltaE/Timestep^2")
for fac in faclist:
    cmdlist = ["-screen","none","-var","dtfac", "%g" % fac]
    lmp = lammps.lammps(cmdargs = cmdlist)
    lmp.file(infile)
    dt = lmp.extract_global("dt", lammps.LAMMPS_DOUBLE)
    de = lmp.extract_fix("avede", lammps.LMP_STYLE_GLOBAL, lammps.LMP_TYPE_SCALAR)
    dedt2 = de/dt**2
    print(f"{dt} {de} {dedt2}")
    
    
