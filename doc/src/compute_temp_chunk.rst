.. index:: compute temp/chunk

compute temp/chunk command
==========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID temp/chunk chunkID value1 value2 ... keyword value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/chunk = style name of this compute command
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
* zero or more values can be listed as value1,value2,etc
* value = *temp* or *kecom* or *internal*
  
  .. parsed-literal::
  
       temp = temperature of each chunk
       kecom = kinetic energy of each chunk based on velocity of center of mass
       internal = internal kinetic energy of each chunk

* zero or more keyword/value pairs may be appended
* keyword = *com* or *bias* or *adof* or *cdof*
  
  .. parsed-literal::
  
       *com* value = *yes* or *no*
         yes = subtract center-of-mass velocity from each chunk before calculating temperature
         no = do not subtract center-of-mass velocity
       *bias* value = bias-ID
         bias-ID = ID of a temperature compute that removes a velocity bias
       *adof* value = dof_per_atom
         dof_per_atom = define this many degrees-of-freedom per atom
       *cdof* value = dof_per_chunk
         dof_per_chunk = define this many degrees-of-freedom per chunk



Examples
""""""""


.. parsed-literal::

   compute 1 fluid temp/chunk molchunk
   compute 1 fluid temp/chunk molchunk temp internal
   compute 1 fluid temp/chunk molchunk bias tpartial adof 2.0

Description
"""""""""""

Define a computation that calculates the temperature of a group of
atoms that are also in chunks, after optionally subtracting out the
center-of-mass velocity of each chunk.  By specifying optional values,
it can also calculate the per-chunk temperature or energies of the
multiple chunks of atoms.

In LAMMPS, chunks are collections of atoms defined by a :doc:`compute chunk/atom <compute_chunk_atom>` command, which assigns each atom
to a single chunk (or no chunk).  The ID for this command is specified
as chunkID.  For example, a single chunk could be the atoms in a
molecule or atoms in a spatial bin.  See the :doc:`compute chunk/atom <compute_chunk_atom>` and :doc:`Howto chunk <Howto_chunk>`
doc pages for details of how chunks can be defined and examples of how
they can be used to measure properties of a system.

The temperature is calculated by the formula KE = DOF/2 k T, where KE =
total kinetic energy of all atoms assigned to chunks (sum of 1/2 m
v\^2), DOF = the total number of degrees of freedom for those atoms, k
= Boltzmann constant, and T = temperature.

The DOF is calculated as N\*adof + Nchunk\*cdof, where N = number of
atoms contributing to the KE, adof = degrees of freedom per atom, and
cdof = degrees of freedom per chunk.  By default adof = 2 or 3 =
dimensionality of system, as set via the :doc:`dimension <dimension>`
command, and cdof = 0.0.  This gives the usual formula for
temperature.

A kinetic energy tensor, stored as a 6-element vector, is also
calculated by this compute for use in the computation of a pressure
tensor.  The formula for the components of the tensor is the same as
the above formula, except that v\^2 is replaced by vx\*vy for the xy
component, etc.  The 6 components of the vector are ordered xx, yy,
zz, xy, xz, yz.

Note that the number of atoms contributing to the temperature is
calculated each time the temperature is evaluated since it is assumed
the atoms may be dynamically assigned to chunks.  Thus there is no
need to use the *dynamic* option of the
:doc:`compute_modify <compute_modify>` command for this compute style.

If any optional values are specified, then per-chunk quantities are
also calculated and stored in a global array, as described below.

The *temp* value calculates the temperature for each chunk by the
formula KE = DOF/2 k T, where KE = total kinetic energy of the chunk
of atoms (sum of 1/2 m v\^2), DOF = the total number of degrees of
freedom for all atoms in the chunk, k = Boltzmann constant, and T =
temperature.

The DOF in this case is calculated as N\*adof + cdof, where N = number
of atoms in the chunk, adof = degrees of freedom per atom, and cdof =
degrees of freedom per chunk.  By default adof = 2 or 3 =
dimensionality of system, as set via the :doc:`dimension <dimension>`
command, and cdof = 0.0.  This gives the usual formula for
temperature.

The *kecom* value calculates the kinetic energy of each chunk as if
all its atoms were moving with the velocity of the center-of-mass of
the chunk.

The *internal* value calculates the internal kinetic energy of each
chunk.  The interal KE is summed over the atoms in the chunk using an
internal "thermal" velocity for each atom, which is its velocity minus
the center-of-mass velocity of the chunk.


----------


Note that currently the global and per-chunk temperatures calculated
by this compute only include translational degrees of freedom for each
atom.  No rotational degrees of freedom are included for finite-size
particles.  Also no degrees of freedom are subtracted for any velocity
bias or constraints that are applied, such as :doc:`compute temp/partial <compute_temp_partial>`, or :doc:`fix shake <fix_shake>`
or :doc:`fix rigid <fix_rigid>`.  This is because those degrees of
freedom (e.g. a constrained bond) could apply to sets of atoms that
are both included and excluded from a specific chunk, and hence the
concept is somewhat ill-defined.  In some cases, you can use the
*adof* and *cdof* keywords to adjust the calculated degrees of freedom
appropriately, as explained below.

Note that the per-chunk temperature calculated by this compute and the
:doc:`fix ave/chunk temp <fix_ave_chunk>` command can be different.
This compute calculates the temperature for each chunk for a single
snapshot.  Fix ave/chunk can do that but can also time average those
values over many snapshots, or it can compute a temperature as if the
atoms in the chunk on different timesteps were collected together as
one set of atoms to calculate their temperature.  This compute allows
the center-of-mass velocity of each chunk to be subtracted before
calculating the temperature; fix ave/chunk does not.

.. note::

   Only atoms in the specified group contribute to the calculations
   performed by this compute.  The :doc:`compute chunk/atom <compute_chunk_atom>` command defines its own group;
   atoms will have a chunk ID = 0 if they are not in that group,
   signifying they are not assigned to a chunk, and will thus also not
   contribute to this calculation.  You can specify the "all" group for
   this command if you simply want to include atoms with non-zero chunk
   IDs.

The simplest way to output the per-chunk results of the compute
temp/chunk calculation to a file is to use the :doc:`fix ave/time <fix_ave_time>` command, for example:


.. parsed-literal::

   compute cc1 all chunk/atom molecule
   compute myChunk all temp/chunk cc1 temp
   fix 1 all ave/time 100 1 100 c_myChunk file tmp.out mode vector


----------


The keyword/value option pairs are used in the following ways.

The *com* keyword can be used with a value of *yes* to subtract the
velocity of the center-of-mass for each chunk from the velocity of the
atoms in that chunk, before calculating either the global or per-chunk
temperature.  This can be useful if the atoms are streaming or
otherwise moving collectively, and you wish to calculate only the
thermal temperature.

For the *bias* keyword, *bias-ID* refers to the ID of a temperature
compute that removes a "bias" velocity from each atom.  This also
allows calculation of the global or per-chunk temperature using only
the thermal temperature of atoms in each chunk after the translational
kinetic energy components have been altered in a prescribed way,
e.g. to remove a velocity profile.  It also applies to the calculation
of the other per-chunk values, such as *kecom* or *internal*\ , which
involve the center-of-mass velocity of each chunk, which is calculated
after the velocity bias is removed from each atom.  Note that the
temperature compute will apply its bias globally to the entire system,
not on a per-chunk basis.

The *adof* and *cdof* keywords define the values used in the degree of
freedom (DOF) formulas used for the global or per-chunk temperature,
as described above.  They can be used to calculate a more appropriate
temperature for some kinds of chunks.  Here are 3 examples:

If spatially binned chunks contain some number of water molecules and
:doc:`fix shake <fix_shake>` is used to make each molecule rigid, then
you could calculate a temperature with 6 degrees of freedom (DOF) (3
translational, 3 rotational) per molecule by setting *adof* to 2.0.

If :doc:`compute temp/partial <compute_temp_partial>` is used with the
*bias* keyword to only allow the x component of velocity to contribute
to the temperature, then *adof* = 1.0 would be appropriate.

If each chunk consists of a large molecule, with some number of its
bonds constrained by :doc:`fix shake <fix_shake>` or the entire molecule
by :doc:`fix rigid/small <fix_rigid>`, *adof* = 0.0 and *cdof* could be
set to the remaining degrees of freedom for the entire molecule
(entire chunk in this case), e.g. 6 for 3d, or 3 for 2d, for a rigid
molecule.


----------


**Output info:**

This compute calculates a global scalar (the temperature) and a global
vector of length 6 (KE tensor), which can be accessed by indices 1-6.
These values can be used by any command that uses global scalar or
vector values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

This compute also optionally calculates a global array, if one or more
of the optional values are specified.  The number of rows in the array
= the number of chunks *Nchunk* as calculated by the specified
:doc:`compute chunk/atom <compute_chunk_atom>` command.  The number of
columns is the number of specified values (1 or more).  These values
can be accessed by any command that uses global array values from a
compute as input.  Again, see the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The scalar value calculated by this compute is "intensive".  The
vector values are "extensive".  The array values are "intensive".

The scalar value will be in temperature :doc:`units <units>`.  The
vector values will be in energy :doc:`units <units>`.  The array values
will be in temperature :doc:`units <units>` for the *temp* value, and in
energy :doc:`units <units>` for the *kecom* and *internal* values.

Restrictions
""""""""""""


The *com* and *bias* keywords cannot be used together.

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`, :doc:`fix ave/chunk temp <fix_ave_chunk>`

Default
"""""""

The option defaults are com no, no bias, adof = dimensionality of the
system (2 or 3), and cdof = 0.0.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
