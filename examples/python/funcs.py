# Python function that implements a loop of short runs
# calls back to LAMMPS via "lmp" instance
# lammps() must be called with ptr=lmpptr for this to work
from __future__ import print_function

def loop(N,cut0,thresh,lmpptr):
  print("LOOP ARGS",N,cut0,thresh,lmpptr)
  from lammps import lammps
  lmp = lammps(ptr=lmpptr)
  natoms = lmp.get_natoms()

  try:
    for i in range(N):
      cut = cut0 + i*0.1

      lmp.set_variable("cut",cut)                 # set a variable in LAMMPS

      lmp.command("pair_style lj/cut ${cut}")     # LAMMPS command
      #lmp.command("pair_style lj/cut %d" % cut)  # LAMMPS command option

      lmp.command("pair_coeff * * 1.0 1.0")       # ditto
      lmp.command("run 10")                       # ditto
      pe = lmp.extract_compute("thermo_pe",0,0)   # extract total PE from LAMMPS
      print("PE",pe/natoms,thresh)
      if pe/natoms < thresh: return
  except Exception as e:
    print("LOOP error:", e)

