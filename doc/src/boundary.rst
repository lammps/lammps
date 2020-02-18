.. index:: boundary

boundary command
================

Syntax
""""""


.. parsed-literal::

   boundary x y z

* x,y,z = *p* or *s* or *f* or *m*\ , one or two letters
  
  .. parsed-literal::
  
       *p* is periodic
       *f* is non-periodic and fixed
       *s* is non-periodic and shrink-wrapped
       *m* is non-periodic and shrink-wrapped with a minimum value



Examples
""""""""


.. parsed-literal::

   boundary p p f
   boundary p fs p
   boundary s f fm

Description
"""""""""""

Set the style of boundaries for the global simulation box in each
dimension.  A single letter assigns the same style to both the lower
and upper face of the box.  Two letters assigns the first style to the
lower face and the second style to the upper face.  The initial size
of the simulation box is set by the :doc:`read_data <read_data>`,
:doc:`read_restart <read_restart>`, or :doc:`create_box <create_box>`
commands.

The style *p* means the box is periodic, so that particles interact
across the boundary, and they can exit one end of the box and re-enter
the other end.  A periodic dimension can change in size due to
constant pressure boundary conditions or box deformation (see the :doc:`fix npt <fix_nh>` and :doc:`fix deform <fix_deform>` commands).  The *p*
style must be applied to both faces of a dimension.

The styles *f*\ , *s*\ , and *m* mean the box is non-periodic, so that
particles do not interact across the boundary and do not move from one
side of the box to the other.

For style *f*\ , the position of the face is fixed.  If an atom moves
outside the face it will be deleted on the next timestep that
reneighboring occurs.  This will typically generate an error unless
you have set the :doc:`thermo_modify lost <thermo_modify>` option to
allow for lost atoms.

For style *s*\ , the position of the face is set so as to encompass the
atoms in that dimension (shrink-wrapping), no matter how far they
move. Note that when the difference between the current box dimensions
and the shrink-wrap box dimensions is large, this can lead to lost
atoms at the beginning of a run when running in parallel. This is due
to the large change in the (global) box dimensions also causing
significant changes in the individual sub-domain sizes. If these
changes are farther than the communication cutoff, atoms will be lost.
This is best addressed by setting initial box dimensions to match the
shrink-wrapped dimensions more closely, by using *m* style boundaries
(see below).

For style *m*\ , shrink-wrapping occurs, but is bounded by the value
specified in the data or restart file or set by the
:doc:`create_box <create_box>` command.  For example, if the upper z
face has a value of 50.0 in the data file, the face will always be
positioned at 50.0 or above, even if the maximum z-extent of all the
atoms becomes less than 50.0.  This can be useful if you start a
simulation with an empty box or if you wish to leave room on one side
of the box, e.g. for atoms to evaporate from a surface.

For triclinic (non-orthogonal) simulation boxes, if the 2nd dimension
of a tilt factor (e.g. y for xy) is periodic, then the periodicity is
enforced with the tilt factor offset.  If the 1st dimension is
shrink-wrapped, then the shrink wrapping is applied to the tilted box
face, to encompass the atoms.  E.g. for a positive xy tilt, the xlo
and xhi faces of the box are planes tilting in the +y direction as y
increases.  These tilted planes are shrink-wrapped around the atoms to
determine the x extent of the box.

See the :doc:`Howto triclinic <Howto_triclinic>` doc page for a
geometric description of triclinic boxes, as defined by LAMMPS, and
how to transform these parameters to and from other commonly used
triclinic representations.

Restrictions
""""""""""""


This command cannot be used after the simulation box is defined by a
:doc:`read_data <read_data>` or :doc:`create_box <create_box>` command or
:doc:`read_restart <read_restart>` command.  See the
:doc:`change_box <change_box>` command for how to change the simulation
box boundaries after it has been defined.

For 2d simulations, the z dimension must be periodic.

Related commands
""""""""""""""""

See the :doc:`thermo_modify <thermo_modify>` command for a discussion
of lost atoms.

Default
"""""""


.. parsed-literal::

   boundary p p p
