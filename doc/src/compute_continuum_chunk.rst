.. index:: compute continuum/chunk

compute continuum/chunk command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID continuum/chunk chunkID cutoff width value1 value2 ... keyword args ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* continuum/chunk = style name of this compute command
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
* cutoff = cutoff for the truncated Gaussian kernel
* width = standard deviation of the corresponding untruncated Gaussian kernel
* one or more input values can be listed
* value = *density*, *volume/fraction*, *momentum/a*, *velocity/a*, *momentum/grad/ab*, *velocity/grad/ab*, *strain/rate/ab*, *stress/ab*, *stress/ke/ab*, *stress/contacts/ab*, *boundary/force/a*, *fabric/ab*, *temperature*

  .. parsed-literal::

       *density* = density field
       *volume/fraction* = volume fraction field
       *momentum/a* = a-component of the momentum field
       *velocity/a* = a-component of the velocity field
       *momentum/grad/ab* = ab-component of the momentum gradient field
       *velocity/grad/ab* = ab-component of the velocity gradient field
       *strain/rate/ab* = ab-component of the strain-rate field
       *stress/ab* = ab-component of the total stress field
       *stress/ke/ab* = ab-component of the kinetic stress field
       *stress/contacts/ab* = ab-component of the contact stress field
       *boundary/force/a* = a-component of the boundary force density
       *fabric/ab* = ab-component of the fabric tensor field
       *temperature* = granular temperature field

* zero or more keyword/arg pairs may be appended
* keyword = *boundary/atom* or *boundary/fix*

  .. parsed-literal::

       *boundary/atom* arg = *groupbound*
         groupbound = group-ID for atoms that make up a boundary
       *boundary/fix* arg = none
         enables boundary corrections from fix wall/gran

Examples
""""""""

.. code-block:: LAMMPS

   compute cc1 all chunk/atom bin/2d x 0.0 1.0 y 0.0 1.0 units box
   compute prop1 all property/chunk cc1 coord1 coord2 count
   compute cont1 all continuum/chunk cc1 1.0 0.3 density velocity/* stress/*
   fix 1 all ave/time 10 10 100 c_prop1[*] c_cont1[*] file continuum.dat mode vector

   compute cc2 all chunk/atom bin/1d z lower 0.5 units reduced
   compute cont2 flow continuum/chunk cc2 2.0 1.0 volume/fraction temperature boundary/force/* boundary/fix

Description
"""""""""""

.. versionadded:: TBD

Define a computation that calculates coarse-grained continuum fields for
chunks of atoms on demand using the construction in
:ref:`(Goldhirsch) <_compute_continuum_chunk_goldhirsch>`.  The fields
are evaluated at the chunk centers with a truncated Gaussian kernel.
This compute does no time averaging and writes no files itself.  To time
average or write the per-chunk data, use :doc:`fix ave/time
<fix_ave_time>` exactly as with other chunk-based compute styles.

In LAMMPS, chunks are collections of atoms defined by a :doc:`compute
chunk/atom <compute_chunk_atom>` command, which assigns each atom to a
single chunk (or no chunk).  The ID of that command is specified as
chunkID.  For example, a chunk may be a molecule or a spatial bin.  See
:doc:`compute chunk/atom <compute_chunk_atom>` and :doc:`Howto chunk
<Howto_chunk>` for details.

This compute is only compatible with the binning styles of
:doc:`compute chunk/atom <compute_chunk_atom>` (*bin/1d*, *bin/2d*, or
*bin/3d*).  If binning is not performed along one of the box dimensions,
all outputs are normalized by the box length in that dimension.  For
example, if a 3d system is binned only along *z*, the reported fields are
normalized by *Lx* and *Ly*.

The available fields include scalar, vector, and tensor quantities.  A
vector field such as the velocity may be requested as an individual
component such as *velocity/x* or with a wildcard such as
*velocity/\**.  This expands to all components (two in 2d, three in
3d).  Tensor fields such as the stress use component names like
*stress/xy*.  Wildcards such as *stress/\** expand to all components
(four in 2d, nine in 3d).  Tensors are not assumed to be symmetric.

Only atoms in the specified group contribute to the calculations.  In
pair sums, both atoms must be in the group.  The referenced
:doc:`compute chunk/atom <compute_chunk_atom>` command defines its own
group and optional region.  Atoms with chunk ID = 0 are not assigned to
any chunk and do not contribute.

----------

The *density* field is

.. math::

   \sum_i m_i W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

where the summation is across all atoms :math:`i` in the chunk. :math:`m_i`
is the atom mass, :math:`\vec{r}_\mathrm{chunk}` is the chunk center,
:math:`\vec{r}_i` is the atom position, and :math:`W` is the kernel.

The *volume/fraction* field is

.. math::

   \sum_i V_i W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

where :math:`V_i` is the finite sized particle volume in 3d or area in 2d.

The *momentum* field is

.. math::

   \sum_i m_i v_{i,a} W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

where :math:`v_{i,a}` is the :math:`a` component of the atom velocity.
The *momentum/grad* field is then obtained with centered finite
differences between neighboring chunks.  Gradient values are zero for
chunks adjacent to a nonperiodic boundary in the corresponding
direction.

The *velocity* field is the ratio of the *momentum* and *density*
fields.  The *velocity/grad* field is then obtained with centered finite
differences.  Thus, if a box dimension is represented by a single chunk,
gradients along that dimension are zero. Gradient values are also zero on
bins that are adjacent to a nonperiodic boundary.

The *boundary/force* field is the interaction force density of
boundaries as defined in :ref:`(Weinhart) <_compute_continuum_chunk_weinhart>`:

.. math::

   \sum_i \sum_k f_{ik,a} W(\vec{r}_\mathrm{chunk} - \vec{r}_{\mathrm{contact},ik})

where the sum over :math:`k` is over boundary elements and
:math:`\vec{r}_{\mathrm{contact},ik}` is the contact point between atom
:math:`i` and boundary element :math:`k`.  At least one of the
*boundary/atom* or *boundary/fix* keywords must be used to request this
quantity.

The *stress/ke* field is the kinetic contribution to the stress:

.. math::

   -\sum_i m_i (v_{i,a} - v_{\mathrm{chunk},a}) (v_{i,b} - v_{\mathrm{chunk},b})
   W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

where :math:`v_{i,a}` is the :math:`a`-th component of the velocity of atom :math:`i` and
:math:`v_{\mathrm{chunk},a}` is the :math:`a`-th component of the average velocity
of the chunk defined by the *velocity* option above.

The *stress/contacts* field is the contact contribution to the stress:

.. math::

   -\sum_{i,j} f_{ij,a} r_{ij,b} \int_0^1 ds\, W(\vec{r}_\mathrm{chunk} -
   \vec{r}_i + s \vec{r}_{ij})

where :math:`f_{ij,a}` is the force on atom :math:`i` from atom
:math:`j` and :math:`\vec{r}_{ij}` is the displacement between the two
atoms.

The *stress* field is the sum of the kinetic and contact contributions.

The *fabric* field is

.. math::

   \sum_{i,j} V_i r_{ij,a} r_{ij,b} \int_0^1 ds\, W(\vec{r}_\mathrm{chunk} -
   \vec{r}_i + s \vec{r}_{ij})

where :math:`V_i` is the volume of the atom in 3D and area in 2D.

The *strain/rate* field is

.. math::

   \frac{1}{2} \left( \grad_{ab} v + \grad_{ba} v \right)

where :math:`\grad_{ab} v` is the :math:`ab` component of the
*velocity/grad* field.

The *temperature* field is a local granular temperature:

.. math::

   \frac{1}{2} \sum_i m_i (v_i - v_\mathrm{chunk})^2
   W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

The optional *boundary/atom* and *boundary/fix* keywords turn on the
boundary corrections for *stress* and *stress/contacts* described in
:ref:`(Weinhart) <_compute_continuum_chunk_weinhart>`.  The
*boundary/atom* keyword designates a group of atoms as a boundary.
Those atoms are removed from the atom and pair sums above, except for the
*boundary/force* contribution.  The *boundary/fix* keyword applies the
analogous correction for boundaries created by
:doc:`fix wall/gran <fix_wall_gran>`.

Output info
""""""""""

This compute calculates a global array where the number of rows is the
number of chunks :math:`N_\text{chunk}` defined by the referenced
:doc:`compute chunk/atom <compute_chunk_atom>` command.  The number of
columns is the number of requested values after wildcard expansion.  The
columns appear in the same order as specified in the command.  Internal
intermediate quantities needed to evaluate derived fields are not exposed
as output columns.

These values can be accessed by any command that uses global arrays from
a compute as input.  See :doc:`Howto output <Howto_output>` for an
overview of output options.  The array values are intensive.  Units
depend on the requested fields.

Restrictions
"""""""""""

This compute is part of the GRANULAR package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more information.

Only *bin/1d*, *bin/2d*, and *bin/3d* styles of
:doc:`compute chunk/atom <compute_chunk_atom>` are supported.

The *volume/fraction*, *stress*, *stress/contacts*, *boundary/force*,
and *fabric* values require particles with a radius attribute.

Pair-dependent quantities require a pair style that supports
``pair->single()``.  The *boundary/fix* keyword requires at least one
:doc:`fix wall/gran <fix_wall_gran>` instance with the *contacts*
keyword enabled.

Related commands
""""""""""""""

:doc:`compute property/chunk <compute_property_chunk>`,
:doc:`compute msd/chunk <compute_msd_chunk>`,
:doc:`fix ave/time <fix_ave_time>`,
:doc:`fix wall/gran <fix_wall_gran>`

Default
"""""""

No boundary corrections are applied unless *boundary/atom* or
*boundary/fix* is specified.

----------

.. _compute_continuum_chunk_goldhirsch:

**(Goldhirsch)** Goldhirsch, Granular Matter, 12, 3, 239-252 (2010).

.. _compute_continuum_chunk_weinhart:

**(Weinhart)** Weinhart, Thornton, Luding, Bokhove, Granular Matter, 14,
2, 289-294 (2012).
