# Example system
The example contains two Fe3O4 magnetic nanoparticles, each functionalized by
19 oleic acid ligands.

Properties of the particles:
* NP saturation magnetization 480 KAmpere/Meter 
* NP volume 1767.15 Angstrom^3
* NP size 1.6 nm

## Specifying dipole moments of the particles
First create a new atom property using LAMMPS `property/atom` command.
Next identify two atoms per nanoparticle that will be used to create a
diameter vector to define the magnetic dipole moment.

```
fix     qm      all property/atom d_qm
```
To calculate the magnetic charge value, use `qm = mu / d` formula, where
`mu` is the dipole moment value in real units, and `d` is the distance
between the two chosen atoms. If the chosen two atoms are located at the 
nanoparticle surface, this will be same as the particle diameter.
To convert SI to real units, as an example, use
`mu = 3.0E-21 A m^-2 / 1.6022E-21 = 1.8724 Ke * Angstrom^2 /fs`

A dipole moment will be defined from the `-qm` atom towards the `+qm` atom.
Currently, each NP can carry only 1 magnetic dipole moment, specified by one
pair of atoms carrying non-zero `qm` values.

```
set     atom 79      d_qm  0.0353   # NP 1
set     atom 140     d_qm -0.0353   # NP 1
set     atom 1347    d_qm  0.0353   # NP 2
set     atom 1283    d_qm -0.0353   # NP 2
```

## External field direction
Direction of the external field can be specified using the field vector
components. In the example, a 1 T uniform external field is applied in
the +z direction using the `bfield 0.0 0.0 1.0` keyword.
