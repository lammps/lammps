2d simulations
==============

Use the :doc:`dimension <dimension>` command to specify a 2d simulation.

Make the simulation box periodic in z via the :doc:`boundary <boundary>`
command.  This is the default.

If using the :doc:`create_box <create_box>` command to define a
simulation box, set the z dimensions narrow, but finite, so that the
:doc:`create_atoms <create_atoms>` command will fill the 3d simulation
box with a single z plane of atoms - e.g.

.. code-block:: LAMMPS

   create box 1 -10 10 -10 10 -0.25 0.25

If using the :doc:`read data <read_data>` command to read in a file of
atom coordinates, set the "zlo zhi" values to be finite but narrow,
similar to the create_box command settings just described.  For each
atom in the file, assign a z coordinate so it falls inside the
z-boundaries of the box - e.g. 0.0.

Use the :doc:`fix enforce2d <fix_enforce2d>` command as the last
defined fix to insure that the z-components of velocities and forces
are zeroed out every timestep.  The reason to make it the last fix is
so that any forces induced by other fixes will be zeroed out.

Many of the example input scripts included in the LAMMPS distribution
are for 2d models.

.. note::

   Some models in LAMMPS treat particles as finite-size spheres, as
   opposed to point particles.  See the :doc:`atom_style sphere <atom_style>` and :doc:`fix nve/sphere <fix_nve_sphere>`
   commands for details.  By default, for 2d simulations, such particles
   will still be modeled as 3d spheres, not 2d discs (circles), meaning
   their moment of inertia will be that of a sphere.  If you wish to
   model them as 2d discs, see the :doc:`set density/disc <set>` command
   and the *disc* option for the :doc:`fix nve/sphere <fix_nve_sphere>`,
   :doc:`fix nvt/sphere <fix_nvt_sphere>`, :doc:`fix nph/sphere <fix_nph_sphere>`, :doc:`fix npt/sphere <fix_npt_sphere>`
   commands.
