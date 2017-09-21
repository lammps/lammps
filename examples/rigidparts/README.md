# Rigid Parts Example

--------------------------------------------------------------------------------
Description
--------------------------------------------------------------------------------

In this example, [Playmol] and [LAMMPS] are used to simulate a single sucrose
molecule solvated in TIP3P water. The solute is modeled as a pair of rigid
bodies (the residues of glucose and fructose) connected together by chemical
bonds to an oxygen atom (a glycosidic linkage).

**NOTE**: This example benefits from ability of [Playmol] to deal with rigid
bodies. An angle, dihedral, or improper is deleted if all involved atoms
belong to the same body. Bonds in the same situation are preserved, but have
`bond_style zero` assigned to them in the generated [LAMMPS] data file.
[Playmol] also includes a section `BodyTags` to this data file, which can be
read in [LAMMPS] using the following commands:

    fix         bodies all property/atom i_bodytag
    read_data   sucrose_in_water.data fix bodies NULL BodyTags

This example also makes use of the following [LAMMPS] commands:

    neigh_modify   exclude custom/match bodytag all
    fix            1 all rigid/small custom bodytag

--------------------------------------------------------------------------------
Requirements
--------------------------------------------------------------------------------

[Playmol] and [LAMMPS] are required for executing this examples.

--------------------------------------------------------------------------------
Preparation (Done!)
--------------------------------------------------------------------------------

File `GLYCAM_06j.playmol` was created from `GLYCAM_06j.dat` using playmoltools:

    playmoltools -f amber -i GLYCAM_06j.dat -o GLYCAM_06j.playmol

File `1.pdb` containing a single sucrose molecule was obtained from Glycam-Web's
[Carbohydrate Builder] using the following condensed code:

    DFrufb2-1DGlcpa

File `sucrose.playmol` was generated from `1.pdb` using playmoltools:

    playmoltools -f pdb -p GLYCAM_06j-1.prep -i 1.pdb -o sucrose.playmol

File `sucrose_in_water.data` was generated from `sucrose_in_water.playmol`
using [Playmol], which creates a simulation box with one sucrose molecule and
800 water molecules:

    playmol sucrose_in_water.playmol

--------------------------------------------------------------------------------
Execution
--------------------------------------------------------------------------------

Execute [LAMMPS] to simulate the system:

    $LAMMPS_DIR/lmp_mpi -in in.sucrose_in_water

--------------------------------------------------------------------------------

[AMBER Tools]:			http://ambermd.org/#AmberTools
[Carbohydrate Builder]:		http://glycam.org/cb
[GLYCAM]:			http://glycam.org/docs/forcefield
[Playmol]:			https://github.com/atoms-ufrj/playmol
[LAMMPS]:			http://lammps.sandia.gov
