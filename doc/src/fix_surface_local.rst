.. index:: fix surface/local

fix surface/local command
===============

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID surface/local source

* ID, group-ID are documented in :doc:`fix <fix>` command
* surface/local = style name of this fix command
* source = NULL or molecule template ID or STL filename

Examples
""""""""

.. code-block:: LAMMPS

   read_data data.tris
   fix 1 all surface/local NULL

   molecule tris surf.tri
   fix 1 all surface/local tris

   fix 1 all surface/local surf.tri.stl

Description
"""""""""""

Enable granular surfaces to be used as boundary conditions on
particles in a granular simulation.  Granular surfaces are defined as
a set of triangles (for 3d models) or a set of line segments (for 2d
models).

The :doc:`Howto granular surfaces <Howto_granular_surfaces>` doc page
gives an overview of granular surfaces of two types, *global* and
*local*, and gives guidelines for how they should be defined.

This command must be used for models with *local* surfaces.
The :doc:`fix surface/global <fix_surface_global>` command must be
used for models with *global* surfaces.  As explained on the
:doc:`Howto granular surfaces <Howto_granular_surfaces>` doc page,
*local* surfaces should typically use triangles/lines whose size is no
more than a few times larger than the spherical particles used in a
granular model.  There can be as many of them as needed to describe
the physical surfaces at high resolution.

*Local* triangles or line segments are distributed across processors
in the same manner as particles, based on which processor's sub-domain
the center of the triangle or line segment is inside of.

*Local* triangles/lines can be defined in one of 3 ways, which
correpond to the 3 options listed above for the *source* argument:

* via a data file, read by the :doc:`read_data <read_data>` command
* via a molecule file(s), read by the :doc:`molecule <molecule>` command
* via an STL file, read by this commmand

If triangles/lines were previously read in by the :doc:`read_data
<read_data>` command, then the *source* argument is specified as NULL.

If triangles/lines were previously read in by the :doc:`molecule
<molecule>` command, then the *source* argument is specified as the
molecule template ID used with the :doc:`molecule <molecule>` command.

STL (stereolithography) files define a set of triangles.  For use with
this command, the *source* argument is specified as the name of the
STL file.  The file can be in text or binary format; this command
auto-detects the format.  Note that STL files cannot be used for 2d
simulations.

This `Wikepedia page
<https://en.wikipedia.org/wiki/STL_(file_format)>`_ describes the
format of both text and binary STL files.  Binary STL files can be
converted to ASCII for editing with the stl_bin2txt tool in the
lammps/tools directory.  Examples of text STL files with the suffix
".stl" are included in the examples/gransurf directory.

If the *source* argument is NULL, a set of distributed triangles or
lines already exist.  As explained on the :doc:`Howto granular
surfaces <Howto_granular_surfaces>` doc page, these are "particles" as
defined by the :doc:`atom_style tri or line <atom_style>` command,
typically as a sub-style of the :doc:`atom_style hybrid <atom_style>`
command.

For the other two *source* options, this command creates a new
triangle or line particle from the information in the molecule
template or STL file.  This is one triangle or line particle for each
triangle or line in the molecule template.  Or one triangle particle
for each triangle in the STL file.

Once all the distributed triangle/line particles are defined, this
command calculates the connectivity of the set of triangles/lines and
stores that information with each triangle/line particle.  Two
triangles are "connected" if they have the same corner point in
common, or the same edge in common (2 corner points).  Two line
segments are "connected" if the they have the same end point in
common.  More technical details on connectivity and its significance
for granular simulations with surfaces is given on :doc:`Howto
granular surfaces <Howto_granular_surfaces>` doc page.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.  No global or per-atom quantities are stored by
this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`fix surface/global <fix_surface_global>`

Default
"""""""

none
