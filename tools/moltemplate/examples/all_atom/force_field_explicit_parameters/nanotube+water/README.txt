This is a small version of a carbon-nanotube, water capillary system.
It was inspired by this paper:

    Laurent Joly, J. Chem. Phys. 135(21):214705 (2011)

-------- Requirements: -------
To run this system at constant pressure, it might help to compile LAMMPS
with the optional RIGID package, and use "fix rigid" on the carbon.
(The use of fix rigid is controversial.)  Running at NVT does not require this.
------------------------------

Note: To investigate the behavior from that paper, it might be a good
      idea to increase the size of the water reservoir, the spacing between
      the walls, and the size of the system in the X and Y directions.

Note: Explicit carbon-carbon bonds:
      In the graphene and nanotube structures, I did not try to connect the
      carbon atoms together with bonds.  Instead we will hold these structures
      rigid by not integrating their equations of motion.
      (If you want to simulate movement of the carbon atoms at high
       temperatures or tension, LAMMPS has 3-body/many-body LAMMPS force-fields
       available for simulating the behaviour of carbon in graphite. I know
       that you don't need to specify bonds to use these force fields.  I do
       not know know if these force fields work for nanotubes or graphene.)

Note: Other modeling tools:
      If you need explicit bonds between carbon atoms, then you must add them
      yourself or use a different tool. Currently (2012-10-20), moltemplate does
      not generate bonds automatically.  The "Nanotube Builder" and "topotools"
      plugins for for VMD can generate a nanotube with bonds in LAMMPS data
      format.  You can then convert this data file to .LT format using the
      ltemplify.py utility and then import it into another .LT file and play
      with it later.  (In the "cnad-cnt" example, the carbon nanotube was built
      using "Nanotube Builder" and topotools, and processed with ltemplify.py)


# WARNING: THIS IS NOT A REALISTIC MODEL OF A GRAPHENE-NANOTUBE JUNCTION.
#           A real junction would be curved and deformed near the boundary,
#           (not 90 degrees) and it would not be built entirely from hexagons.
#           (This is not a problem in this example because the carbon atoms
#           are immobilized.)  If you want to simulate the behavior of
#           real graphene or nanotube junctions, you must be more careful.
#          To solve this problem:
#           Moltemplate allows you to move, customize or delete individual
#           atoms near the boundary.  You can move atoms by overwriting their
#           coordinates using additional write("Data Atoms") statements (after
#           the walls and tube are created).  You can also change their charge.
#          Alternately, you could start with the structure provided here, and
#           relax/minimize the coordinates of the carbon atoms using LAMMPS
#           before using it in other simulations.
#          Or you could do both (customization & minimization).
