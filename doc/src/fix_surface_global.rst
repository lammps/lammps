.. index:: fix surface/global

fix surface/global command
===============

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID surface/global source fstyle fstyle_params

* ID, group-ID are documented in :doc:`fix <fix>` command
* surface/global = style name of this fix command
* source = molecule template ID or STL filename
* fstyle = style of force interactions between particles and surfaces

  .. parsed-literal::

       possible choices: hooke, hooke/history, hertz/history, granular

* fstyle_params = parameters associated with force interaction style

  .. parsed-literal::

       For *hooke*, *hooke/history*, and *hertz/history*, *fstyle_params* are:
             Kn = elastic constant for normal particle repulsion (force/distance units or pressure units - see discussion below)
             Kt = elastic constant for tangential contact (force/distance units or pressure units - see discussion below)
             gamma_n = damping coefficient for collisions in normal direction (1/time units or 1/time-distance units - see discussion below)
             gamma_t = damping coefficient for collisions in tangential direction (1/time units or 1/time-distance units - see discussion below)
             xmu = static yield criterion (unitless value between 0.0 and 1.0e4)
             dampflag = 0 or 1 if tangential damping force is excluded or included
             optional keyword = *limit_damping*, limit damping to prevent attractive interaction

  .. parsed-literal::

       For *granular*, *fstyle_params* are set using the same syntax as for the *pair_coeff* command of :doc:`pair_style granular <pair_granular>`

  
Examples
""""""""

.. code-block:: LAMMPS

   molecule lines surf.line
   fix 1 all surface/global lines hooke 4000.0 NULL 100.0 NULL 0.5 1

   fix 1 all surface/global surf.tri.sgl hooke/history 4000.0 NULL 100.0 NULL 0.5 1

Description
"""""""""""

Enable granular surfaces to be used as boundary conditions on
particles in a granular simulation.  Granular surfaces are defined as
a set of triangles (for 3d models) or a set of line segments (for 2d
models).

The :doc:`Howto granular surfaces <Howto_granular_surfaces>` doc page
gives an overview of granular surfaces of two types, *global* and
*local*, and gives guidelines for how they should be defined.

This command must be used for models with *global* surfaces.  The
:doc:`fix surface/local <fix_surface_local>` command must be used for
models with *local* surfaces.  As explained on the :doc:`Howto
granular surfaces <Howto_granular_surfaces>` doc page, *global*
surfaces are most appropriate when there is a modest number of them.
Each triangle/line can be of any size, even as large as a dimension of
the simulation box.

A copy of *global* triangles or line segments are stored on all
processors.

*Global* triangles/lines can be defined in one of 2 ways, which
correpond to the 2 options listed above for the *source* argument:

* via a molecule file(s), read by the :doc:`molecule <molecule>` command
* via an STL file, read by this commmand

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

Once the *global* triangles/lines are defined, this command calculates
the connectivity of the set of triangles/lines and stores that
information as well.  Two triangles are "connected" if they have the
same corner point in common, or the same edge in common (2 corner
points).  Two line segments are "connected" if the they have the same
end point in common.  More technical details on connectivity and its
significance for granular simulations with surfaces is given on
:doc:`Howto granular surfaces <Howto_granular_surfaces>` doc page.

Note that all particles in the system interact with the surface when
they are close enough to touch any individual triangle/line.

NOTE: change this to be just particles in the specified group ?

The nature of the surface/particle interactions are determined by the
*fstyle* setting.  It can be any of the styles defined by the
:doc:`pair_style gran/\* <pair_gran>` or the more general
:doc:`pair_style granular <pair_granular>` commands.  Currently the
options are *hooke*, *hooke/history*, or *hertz/history* for the
former, and *granular* with all the possible options of the associated
*pair_coeff* command for the latter.  The equation for the force
between a triangle/line and a particle touching it is the same as the
corresponding equation on the :doc:`pair_style gran/\* <pair_gran>`
and :doc:`pair_style granular <pair_granular>` doc pages, in the limit
of one of the two particles going to infinite radius and mass (flat
surface).  Specifically, delta = radius - r = overlap of particle with
triangle/line, m_eff = mass of particle, and the effective radius of
contact = RiRj/Ri+Rj is set to the radius of the particle.

The parameters *Kn*, *Kt*, *gamma_n*, *gamma_t*, *xmu*, *dampflag*,
and the optional keyword *limit_damping* have the same meaning and
units as those specified with the :doc:`pair_style gran/\*
<pair_gran>` commands.  This means a NULL can be used for either *Kt*
or *gamma_t* as described on that page.  If a NULL is used for *Kt*,
then a default value is used where *Kt* = 2/7 *Kn*\ .  If a NULL is
used for *gamma_t*, then a default value is used where *gamma_t* = 1/2
*gamma_n*.

NOTE: Include more info here from the fix wall/gran doc page on the
topic of surface/particle interactions ?

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

The :doc:`fix_modify move <fix_modify>` option is supported by this
fix.  Its syntax is as follows:

.. code-block:: LAMMPS

   move mode args

* mode = *rotate* (other modes of surface motion will be added later)

  .. parsed-literal::

       *mode* args = Px Py Pz Rx Ry Rz period
         Px,Py,Pz = origin point of axis of rotation (distance units)
         Rx,Ry,Rz = axis of rotation vector
         period = period of rotation (time units)

As an example:

.. code-block:: LAMMPS

fix_modify 1 move rotate 0 0 0 0 0 1 25

No global or per-atom quantities are stored by this fix for access by
various :doc:`output commands <Howto_output>`.  No parameter of this
fix can be used with the *start/stop* keywords of the :doc:`run <run>`
command.  This fix is not invoked during :doc:`energy minimization
<minimize>`.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`fix surface/local <fix_surface_local>`

Default
"""""""

none
