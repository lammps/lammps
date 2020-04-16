.. index:: change_box

change_box command
==================

Syntax
""""""

.. code-block:: LAMMPS

   change_box group-ID parameter args ... keyword args ...

* group-ID = ID of group of atoms to (optionally) displace
* one or more parameter/arg pairs may be appended

  .. parsed-literal::

     parameter = *x* or *y* or *z* or *xy* or *xz* or *yz* or *boundary* or *ortho* or *triclinic* or *set* or *remap*
       *x*\ , *y*\ , *z* args = style value(s)
         style = *final* or *delta* or *scale* or *volume*
           *final* values = lo hi
             lo hi = box boundaries after displacement (distance units)
           *delta* values = dlo dhi
             dlo dhi = change in box boundaries after displacement (distance units)
           *scale* values = factor
             factor = multiplicative factor for change in box length after displacement
           *volume* value = none = adjust this dim to preserve volume of system
       *xy*\ , *xz*\ , *yz* args = style value
         style = *final* or *delta*
           *final* value = tilt
             tilt = tilt factor after displacement (distance units)
           *delta* value = dtilt
             dtilt = change in tilt factor after displacement (distance units)
       *boundary* args = x y z
         x,y,z = *p* or *s* or *f* or *m*\ , one or two letters
         *p* is periodic
         *f* is non-periodic and fixed
         *s* is non-periodic and shrink-wrapped
         *m* is non-periodic and shrink-wrapped with a minimum value
       *ortho* args = none = change box to orthogonal
       *triclinic* args = none = change box to triclinic
       *set* args = none = store state of current box
       *remap* args = none = remap atom coords from last saved state to current box

* zero or more keyword/value pairs may be appended
* keyword = *units*

  .. parsed-literal::

       *units* value = *lattice* or *box*
         lattice = distances are defined in lattice units
         box = distances are defined in simulation box units

Examples
""""""""

.. code-block:: LAMMPS

   change_box all xy final -2.0 z final 0.0 5.0 boundary p p f remap units box
   change_box all x scale 1.1 y volume z volume remap

Description
"""""""""""

Change the volume and/or shape and/or boundary conditions for the
simulation box.  Orthogonal simulation boxes have 3 adjustable size
parameters (x,y,z).  Triclinic (non-orthogonal) simulation boxes have
6 adjustable size/shape parameters (x,y,z,xy,xz,yz).  Any or all of
them can be adjusted independently by this command.  Thus it can be
used to expand or contract a box, or to apply a shear strain to a
non-orthogonal box.  It can also be used to change the boundary
conditions for the simulation box, similar to the
:doc:`boundary <boundary>` command.

The size and shape of the initial simulation box are specified by the
:doc:`create_box <create_box>` or :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` command used to setup the simulation.
The size and shape may be altered by subsequent runs, e.g. by use of
the :doc:`fix npt <fix_nh>` or :doc:`fix deform <fix_deform>` commands.
The :doc:`create_box <create_box>`, :doc:`read data <read_data>`, and
:doc:`read_restart <read_restart>` commands also determine whether the
simulation box is orthogonal or triclinic and their doc pages explain
the meaning of the xy,xz,yz tilt factors.

See the :doc:`Howto triclinic <Howto_triclinic>` doc page for a
geometric description of triclinic boxes, as defined by LAMMPS, and
how to transform these parameters to and from other commonly used
triclinic representations.

The keywords used in this command are applied sequentially to the
simulation box and the atoms in it, in the order specified.

Before the sequence of keywords are invoked, the current box
size/shape is stored, in case a *remap* keyword is used to map the
atom coordinates from a previously stored box size/shape to the
current one.

After all the keywords have been processed, any shrink-wrap boundary
conditions are invoked (see the :doc:`boundary <boundary>` command)
which may change simulation box boundaries, and atoms are migrated to
new owning processors.

.. note::

   This means that you cannot use the change_box command to enlarge
   a shrink-wrapped box, e.g. to make room to insert more atoms via the
   :doc:`create_atoms <create_atoms>` command, because the simulation box
   will be re-shrink-wrapped before the change_box command completes.
   Instead you could do something like this, assuming the simulation box
   is non-periodic and atoms extend from 0 to 20 in all dimensions:

.. code-block:: LAMMPS

   change_box all x final -10 20
   create_atoms 1 single -5 5 5       # this will fail to insert an atom

   change_box all x final -10 20 boundary f s s
   create_atoms 1 single -5 5 5
   change_box all boundary s s s      # this will work

.. note::

   Unlike the earlier "displace_box" version of this command, atom
   remapping is NOT performed by default.  This command allows remapping
   to be done in a more general way, exactly when you specify it (zero or
   more times) in the sequence of transformations.  Thus if you do not
   use the *remap* keyword, atom coordinates will not be changed even if
   the box size/shape changes.  If a uniformly strained state is desired,
   the *remap* keyword should be specified.

.. note::

   It is possible to lose atoms with this command.  E.g. by
   changing the box without remapping the atoms, and having atoms end up
   outside of non-periodic boundaries.  It is also possible to alter
   bonds between atoms straddling a boundary in bad ways.  E.g. by
   converting a boundary from periodic to non-periodic.  It is also
   possible when remapping atoms to put them (nearly) on top of each
   other.  E.g. by converting a boundary from non-periodic to periodic.
   All of these will typically lead to bad dynamics and/or generate error
   messages.

.. note::

   The simulation box size/shape can be changed by arbitrarily large
   amounts by this command.  This is not a problem, except that the
   mapping of processors to the simulation box is not changed from its
   initial 3d configuration; see the :doc:`processors <processors>`
   command.  Thus, if the box size/shape changes dramatically, the
   mapping of processors to the simulation box may not end up as
   optimal as the initial mapping attempted to be.  You may wish to
   re-balance the atoms by using the :doc:`balance <balance>` command
   if that is the case.

.. note::

   You cannot use this command after reading a restart file (and
   before a run is performed) if the restart file stored per-atom
   information from a fix and any of the specified keywords change the
   box size or shape or boundary conditions.  This is because atoms
   may be moved to new processors and the restart info will not
   migrate with them.  LAMMPS will generate an error if this could
   happen.  Only the *ortho* and *triclinic* keywords do not trigger
   this error.  One solution is to perform a "run 0" command before
   using the change_box command.  This clears the per-atom restart
   data, whether it has been re-assigned to a new fix or not.

.. note::

   Because the keywords used in this command are applied one at a time
   to the simulation box and the atoms in it, care must be taken with
   triclinic cells to avoid exceeding the limits on skew after each
   transformation in the sequence.  If skew is exceeded before the
   final transformation this can be avoided by changing the order of
   the sequence, or breaking the transformation into two or more
   smaller transformations.  For more information on the allowed
   limits for box skew see the discussion on triclinic boxes on
   :doc:`Howto triclinic <Howto_triclinic>` doc page.

----------

For the *x*\ , *y*\ , and *z* parameters, this is the meaning of their
styles and values.

For style *final*\ , the final lo and hi box boundaries of a dimension
are specified.  The values can be in lattice or box distance units.
See the discussion of the units keyword below.

For style *delta*\ , plus or minus changes in the lo/hi box boundaries
of a dimension are specified.  The values can be in lattice or box
distance units.  See the discussion of the units keyword below.

For style *scale*\ , a multiplicative factor to apply to the box length
of a dimension is specified.  For example, if the initial box length
is 10, and the factor is 1.1, then the final box length will be 11.  A
factor less than 1.0 means compression.

The *volume* style changes the specified dimension in such a way that
the overall box volume remains constant with respect to the operation
performed by the preceding keyword.  The *volume* style can only be
used following a keyword that changed the volume, which is any of the
*x*\ , *y*\ , *z* keywords.  If the preceding keyword "key" had a *volume*
style, then both it and the current keyword apply to the keyword
preceding "key".  I.e. this sequence of keywords is allowed:

.. code-block:: LAMMPS

   change_box all x scale 1.1 y volume z volume

The *volume* style changes the associated dimension so that the
overall box volume is unchanged relative to its value before the
preceding keyword was invoked.

If the following command is used, then the z box length will shrink by
the same 1.1 factor the x box length was increased by:

.. code-block:: LAMMPS

   change_box all x scale 1.1 z volume

If the following command is used, then the y,z box lengths will each
shrink by sqrt(1.1) to keep the volume constant.  In this case, the
y,z box lengths shrink so as to keep their relative aspect ratio
constant:

.. code-block:: LAMMPS

   change_box all x scale 1.1 y volume z volume

If the following command is used, then the final box will be a factor
of 10% larger in x and y, and a factor of 21% smaller in z, so as to
keep the volume constant:

.. code-block:: LAMMPS

   change_box all x scale 1.1 z volume y scale 1.1 z volume

.. note::

   For solids or liquids, when one dimension of the box is
   expanded, it may be physically undesirable to hold the other 2 box
   lengths constant since that implies a density change.  For solids,
   adjusting the other dimensions via the *volume* style may make
   physical sense (just as for a liquid), but may not be correct for
   materials and potentials whose Poisson ratio is not 0.5.

For the *scale* and *volume* styles, the box length is expanded or
compressed around its mid point.

----------

For the *xy*\ , *xz*\ , and *yz* parameters, this is the meaning of their
styles and values.  Note that changing the tilt factors of a triclinic
box does not change its volume.

For style *final*\ , the final tilt factor is specified.  The value
can be in lattice or box distance units.  See the discussion of the
units keyword below.

For style *delta*\ , a plus or minus change in the tilt factor is
specified.  The value can be in lattice or box distance units.  See
the discussion of the units keyword below.

All of these styles change the xy, xz, yz tilt factors.  In LAMMPS,
tilt factors (xy,xz,yz) for triclinic boxes are required to be no more
than half the distance of the parallel box length.  For example, if
xlo = 2 and xhi = 12, then the x box length is 10 and the xy tilt
factor must be between -5 and 5.  Similarly, both xz and yz must be
between -(xhi-xlo)/2 and +(yhi-ylo)/2.  Note that this is not a
limitation, since if the maximum tilt factor is 5 (as in this
example), then configurations with tilt = ..., -15, -5, 5, 15, 25,
... are all equivalent.  Any tilt factor specified by this command
must be within these limits.

----------

The *boundary* keyword takes arguments that have exactly the same
meaning as they do for the :doc:`boundary <boundary>` command.  In each
dimension, a single letter assigns the same style to both the lower
and upper face of the box.  Two letters assigns the first style to the
lower face and the second style to the upper face.

The style *p* means the box is periodic; the other styles mean
non-periodic. For style *f*\ , the position of the face is fixed.  For
style *s*\ , the position of the face is set so as to encompass the
atoms in that dimension (shrink-wrapping), no matter how far they
move.  For style *m*\ , shrink-wrapping occurs, but is bounded by the
current box edge in that dimension, so that the box will become no
smaller.  See the :doc:`boundary <boundary>` command for more
explanation of these style options.

Note that the "boundary" command itself can only be used before the
simulation box is defined via a :doc:`read_data <read_data>` or
:doc:`create_box <create_box>` or :doc:`read_restart <read_restart>`
command.  This command allows the boundary conditions to be changed
later in your input script.  Also note that the
:doc:`read_restart <read_restart>` will change boundary conditions to
match what is stored in the restart file.  So if you wish to change
them, you should use the change_box command after the read_restart
command.

----------

The *ortho* and *triclinic* keywords convert the simulation box to be
orthogonal or triclinic (non-orthogonal).

The simulation box is defined as either orthogonal or triclinic when
it is created via the :doc:`create_box <create_box>`,
:doc:`read_data <read_data>`, or :doc:`read_restart <read_restart>`
commands.

These keywords allow you to toggle the existing simulation box from
orthogonal to triclinic and vice versa.  For example, an initial
equilibration simulation can be run in an orthogonal box, the box can
be toggled to triclinic, and then a :doc:`non-equilibrium MD (NEMD) simulation <Howto_nemd>` can be run with deformation via the :doc:`fix deform <fix_deform>` command.

If the simulation box is currently triclinic and has non-zero tilt in
xy, yz, or xz, then it cannot be converted to an orthogonal box.

----------

The *set* keyword saves the current box size/shape.  This can be
useful if you wish to use the *remap* keyword more than once or if you
wish it to be applied to an intermediate box size/shape in a sequence
of keyword operations.  Note that the box size/shape is saved before
any of the keywords are processed, i.e. the box size/shape at the time
the create_box command is encountered in the input script.

The *remap* keyword remaps atom coordinates from the last saved box
size/shape to the current box state.  For example, if you stretch the
box in the x dimension or tilt it in the xy plane via the *x* and *xy*
keywords, then the *remap* command will dilate or tilt the atoms to
conform to the new box size/shape, as if the atoms moved with the box
as it deformed.

Note that this operation is performed without regard to periodic
boundaries.  Also, any shrink-wrapping of non-periodic boundaries (see
the :doc:`boundary <boundary>` command) occurs after all keywords,
including this one, have been processed.

Only atoms in the specified group are remapped.

----------

The *units* keyword determines the meaning of the distance units used
to define various arguments.  A *box* value selects standard distance
units as defined by the :doc:`units <units>` command, e.g. Angstroms for
units = real or metal.  A *lattice* value means the distance units are
in lattice spacings.  The :doc:`lattice <lattice>` command must have
been previously used to define the lattice spacing.

----------

Restrictions
""""""""""""

If you use the *ortho* or *triclinic* keywords, then at the point in
the input script when this command is issued, no :doc:`dumps <dump>` can
be active, nor can a :doc:`fix deform <fix_deform>` be active.  This is
because these commands test whether the simulation box is orthogonal
when they are first issued.  Note that these commands can be used in
your script before a change_box command is issued, so long as an
:doc:`undump <undump>` or :doc:`unfix <unfix>` command is also used to
turn them off.

Related commands
""""""""""""""""

:doc:`fix deform <fix_deform>`, :doc:`boundary <boundary>`

Default
"""""""

The option default is units = lattice.
