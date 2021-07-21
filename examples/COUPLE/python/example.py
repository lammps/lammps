# this example requires the LAMMPS Python package (python/lammps) to be installed
# and LAMMPS to be loadable as shared library in LD_LIBRARY_PATH

import lammps

def callback(caller, ntimestep, nlocal, tag, x, fext):
    """
    This callback receives a caller object that was setup when registering the callback

    In addition to timestep and number of local atoms, the tag and x arrays are passed as
    NumPy arrays. The fext array is a force array allocated for fix external, which
    can be used to apply forces to all atoms. Simply update the value in the array,
    it will be directly written into the LAMMPS C arrays
    """
    print("Data passed by caller (optional)", caller)
    print("Timestep:", ntimestep)
    print("Number of Atoms:", nlocal)
    print("Atom Tags:", tag)
    print("Atom Positions:", x)
    print("Force Additions:", fext)
    fext.fill(1.0)
    print("Force additions after update:", fext)
    print("="*40)

L = lammps.lammps()
L.file("in.fix_external")

# you can pass an arbitrary Python object to the callback every time it is called
# this can be useful if you need more state information such as the LAMMPS ptr to
# make additional library calls
custom_object = ["Some data", L]

L.set_fix_external_callback("2", callback, custom_object)
L.command("run 100")


