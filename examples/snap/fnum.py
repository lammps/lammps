import lammps

def sqerr(a,b):
    return (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2

infile = "in.fnum"

fdeltalist = [1.0e-2,1.0e-3,1.0e-4,1.0e-5,1.0e-6,1.0e-7,1.0e-8,1.0e-9,1.0e-10]

print("Fdelta RMSE")
for fdelta in fdeltalist:
    cmdlist = ["-screen","none","-var","fdelta",f'{fdelta}']
    lmp = lammps.lammps(cmdargs = cmdlist)
    lmp.file(infile)
    natoms = lmp.get_natoms()

    f = lmp.extract_atom("f")
    fnum = lmp.extract_fix("fnum", lammps.LMP_STYLE_ATOM, lammps.LMP_TYPE_ARRAY)

    sumsq = 0
    for i in range(nlocal):
        sumsq += sqerr(fnum[i],f[i])
    rmse = (sumsq/natoms)**0.5
    print(f"{fdelta} {rmse}")
