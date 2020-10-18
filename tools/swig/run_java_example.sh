#!/bin/sh

cat > example.java <<EOF
public class example {
  static {
    System.loadLibrary("javalammps");
  }

  public static void main(String argv[]) {
    SWIGTYPE_p_void lmp = javalammps.lammps_open_no_mpi(0, null, null);
    int ver = javalammps.lammps_version(lmp);

    javalammps.lammps_command(lmp, "units real");
    javalammps.lammps_command(lmp, "lattice fcc 2.5");
    javalammps.lammps_command(lmp, "region box block -5 5 -5 5 -5 5");
    javalammps.lammps_command(lmp, "create_box 1 box");
    javalammps.lammps_command(lmp, "create_atoms 1 box");

    System.out.println("LAMMPS version " + ver);
    System.out.println("Number of created atoms: " + javalammps.lammps_get_natoms(lmp));
    javalammps.lammps_close(lmp);
  }
}
EOF

CLASSPATH=$PWD:${CLASSPATH-${PWD}}
LD_LIBRARY_PATH=$PWD:${LD_LIBARARY_PATH-${PWD}}

export CLASSPATH LD_LIBRARY_PATH

javac *.java
java example
