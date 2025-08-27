.. index:: fix continuum/chunk

fix continuum/chunk command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID continuum/chunk Nevery Nrepeat Nfreq chunkID cutoff width value1 value2 ... keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* ave/chunk = style name of this fix command
* Nevery = use input values every this many timesteps
* Nrepeat = # of times to use input values for calculating averages
* Nfreq = calculate averages every this many timesteps
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
* cutoff = cutoff for Gaussian kernel
* width = standard deviation of Gaussian kernel
* one or more input values can be listed
* value = *density*, *volume/fraction*, *v/a*, *boundary/force/a*, *stress/ab*, *stress/ke/ab*, *stress/contacts/ab*, *fabric/ab*

  .. parsed-literal::

       *density* = density field
       *volume/fraction* = volume fraction field
       *v/a* = a-component of the velocity field
       *boundary/force/a* = a-component of the boundary force density
       *stress/ab* = ab-component of the total stress field
       *stress/ke/ab* = ab-component of the kinetic stress field
       *stress/contacts/ab* = ab-component of the contacts stress field
       *fabric/ab* = ab-component of the fabric tensor field

* zero or more keyword/arg pairs may be appended
* keyword = *ave* or *boundary/atom* or *boundary/fix* or *file* or *append* or *overwrite* or *format* or *title1* or *title2* or *title3*

  .. parsed-literal::

       *ave* args = *one* or *running* or *window M*
         one = output new average value every Nfreq steps
         running = output cumulative average of all previous Nfreq steps
         window M = output average of M most recent Nfreq steps
       *boundary/atom* arg = *groupbound* = enable boundary corrections from atoms
         groupbound = group-ID for atoms that make up a boundary
       *boundary/fix* arg = none = enable boundary corrections from fixes
       *file* arg = filename
         filename = file to write results to
       *append* arg = filename
         filename = file to append results to
       *overwrite* arg = none = overwrite output file with only latest output
       *format* arg = string
         string = C-style format string
       *title1* arg = string
         string = text to print as 1st line of output file
       *title2* arg = string
         string = text to print as 2nd line of output file
       *title3* arg = string
         string = text to print as 3rd line of output file

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all continuum/chunk 10000 1 10000 binchunk c_myCentro 2.0 1.0 density stress/* title1 "My output values"
   fix 1 flow continuum/chunk 100 10 1000 molchunk 3.0 2.0 v/x v/y file vel.profile


Description
"""""""""""

Calculate coarse-grained fields of one or more values every few timesteps
using the construction in :ref:`(Goldhirsch) <_fix_continuum_chunk_goldhirsch>`
using chunks. Coarse graining is performed using a Gaussian kernel. The *width*
of the kernel is the typical standard deviation of the corresponding
non-truncated Gaussian distribution (an infinite cutoff). Fields are then
averaged over longer timescales. The resulting chunk fields can be used by
other :doc:`output commands <Howto_output>` such as :doc:`thermo_style custom
<thermo_style>`, and can also be written to a file.

The available fields include scalar, vector, and tensor values described below.
An individual component of a vector field such as the velocity can be requested
as *v/x*, *v/y*, or *v/z* or one can request all three components (two in 2D)
using a wildcard: *v/\**. Tensor fields, such as the stress, and require two
components to be specified. E.g., *stress/xx*, *stress/xy*, *stress/yx*, and
*stress/yy* in 2D. Note that the fix does not assume tensors are symmetric
and independently calculates each diagonal component. A wildcard can be used
to access all nine components (four in 2D) such as *stress/\**.

In LAMMPS, chunks are collections of atoms defined by a :doc:`compute
chunk/atom <compute_chunk_atom>` command, which assigns each atom to a
single chunk (or no chunk).  The ID for this command is specified as
chunkID.  For example, a single chunk could be the atoms in a molecule
or atoms in a spatial bin.  See the :doc:`compute chunk/atom
<compute_chunk_atom>` page and the :doc:`Howto chunk <Howto_chunk>`
page for details of how chunks can be defined and examples of how they
can be used to measure properties of a system.

This fix is ONLY compatible with the binning styles in compute chunk/atom
(*bin/1d*, *bin/2d*, or *bin/3d*). This fix also requires bins to be
wider than the specified *cutoff* in all dimensions. If binning is not
performed along one of the dimensions of the box, then all outputs are
normalized by the length of the simulation box in that dimension. For
instance, if a 3D system is only binned along the z dimension using
*bin/1d*, then all outputs are normalized by *Lx* and *Ly*.

Note that only atoms in the specified group contribute to the calculations.
In sums of pairs of atoms, both atoms must be in the group to contribute.
The :doc:`compute chunk/atom <compute_chunk_atom>` command defines its own
group as well as an optional region.  Atoms will have a chunk ID = 0,
meaning they belong to no chunk, if they are not in that group or region.
Thus you can specify the "all" group for this command if you simply want
to use the chunk definitions provided by chunkID.

----------

The :math:`N_\text{every}`, :math:`N_\text{repeat}`, and :math:`N_\text{freq}`
arguments specify on what time steps the input values will be accessed and
contribute to the average.  The final averaged quantities are generated on time
steps that are a multiples of :math:`N_\text{freq}`\ .  The average is over
:math:`N_\text{repeat}` quantities, computed in the preceding portion of the
simulation every :math:`N_\text{every}` time steps.  :math:`N_\text{freq}`
must be a multiple of :math:`N_\text{every}` and :math:`N_\text{every}` must be
non-zero even if :math:`N_\text{repeat} = 1`\ .  Also, the time steps
contributing to the average value cannot overlap (i.e.,
:math:`N_\text{repeat} \times N_\text{every}` cannot exceed :math:`N_\text{freq}`).

For example, if :math:`N_\text{every}=2`, :math:`N_\text{repeat}=6`, and
:math:`N_\text{freq}=100`, then values on
time steps 90, 92, 94, 96, 98, 100 will be used to compute the final average
on time step 100.  Similarly for time steps 190, 192, 194, 196, 198, 200 on
time step 200, etc.  If :math:`N_\text{repeat}=1` and
:math:`N_\text{freq} = 100`, then no time averaging is done; values are simply
generated on time steps 100, 200, etc.

.. note::

   To perform per-chunk averaging within a :math:`N_\text{freq}` time window,
   the number of chunks :math:`N_\text{chunk}` defined by the
   :doc:`compute chunk/atom <compute_chunk_atom>` command must remain
   constant.  If the *ave* keyword is set to *running* or *window* then
   :math:`N_\text{chunk}` must remain constant for the duration of the
   simulation.  This fix forces the chunk/atom compute specified by chunkID to
   hold :math:`N_\text{chunk}` constant for the appropriate time windows,
   by not allowing it to re-calculate :math:`N_\text{chunk}`, which can also
   affect how it assigns chunk IDs to atoms.  This is particularly important to
   understand if the chunks defined by the :doc:`compute chunk/atom
   <compute_chunk_atom>` command are spatial bins.  If its *units*
   keyword is set to *box* or *lattice*, then the number of bins
   :math:`N_\text{chunk}` and size of each bin will be fixed over the
   :math:`N_\text{freq}` time window, which can affect which atoms are
   discarded if the simulation box size changes.  If its *units* keyword is set
   to *reduced*, then the number of bins :math:`N_\text{chunk}` will still be
   fixed, but the size of each bin can vary at each time step if the
   simulation box size changes (e.g., for an NPT simulation).

----------

Additional optional keywords also affect the operation of this fix
and its outputs.

----------

The *ave* keyword determines how the per-chunk values produced every
:math:`N_\text{freq}` steps are averaged with values produced on previous steps
that were multiples of :math:`N_\text{freq}`, before they are accessed by
another output command or written to a file.

If the *ave* setting is *one*, which is the default, then the chunk
values produced on timesteps that are multiples of :math:`N_\text{freq}` are
independent of each other; they are output as-is without further averaging.

If the *ave* setting is *running*, then the chunk values produced on
timesteps that are multiples of :math:`N_\text{freq}` are summed and averaged
in a cumulative sense before being output.  Each output chunk value is thus
the average of the chunk value produced on that timestep with all
preceding values for the same chunk.  This running average begins when
the fix is defined; it can only be restarted by deleting the fix via
the :doc:`unfix <unfix>` command, or re-defining the fix by re-specifying it.

If the *ave* setting is *window*, then the chunk values produced on
timesteps that are multiples of :math:`N_\text{freq}` are summed and averaged
within a moving "window" of time, so that the last :math:`M` values for the
same chunk are used to produce the output.  For example, if :math:`M = 3` and
:math:`N_\text{freq} = 1000`, then the output on step 10000 will be the average
of the individual chunk values on time steps 8000, 9000, and 10000.  Outputs on
early steps will average over less than :math:`M` values if they are not
available.

----------

The *density* field is given by

.. math::

   \sum_i m_i W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

where the summation is across all atoms :math:`i` in the chunk. :math:`m_i`
is the atom mass, :math:`\vec{r}_\mathrm{chunk}` is the location of the center of the
chunk, :math:`\vec{r}_i` is the atom position, and :math:`W` is the Gaussian
kernel.

The *volume/fraction* field is given by

.. math::

   \sum_i V_i W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

where :math:`V_i` is the volume of the finite sized atom in 3D and the area in 2D.

The *momentum* field is given by

.. math::

   \sum_i m_i v_{i,a} W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

where :math:`v_{i,a}` is the :math:`a`-component of the atom velocity.

The *velocity* field is calculated as the ratio of the *momentum* and *density*
fields.

The *boundary/force* field is the interaction force density of boundaries as
defined in :ref:`(Weinhart)`. It is calculated as

.. math::

   \sum_i \sum_k f_{ik,a} W(\vec{r}_\mathrm{chunk} - \vec{r}_{\mathrm{contact},ik})

where the sum over :math:`k` is a summation across boundary elements and
:math:`\vec{r}_{\mathrm{contact},ik}` is the position of the contact between
atom :math:`i` and boundary element :math:`k`. A boundary element can either
be an atom or a fix as defined below in the discussion of the *boundary* options.
At least one of the two *boundary* options must be enabled to compute this quantity.

The *stress/ke* field is the kinetic contribution to the stress and is
given by

.. math::

   -\sum_i m_i v_{i,a} v_{i,b} W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

The *stress/contacts* field is the contact force contribution to the stress
and is given by

.. math::

   -\sum_{i,j} f_{ij,a} r_{ij,b} \int_0^1 ds W(\vec{r}_\mathrm{chunk} - \vec{r}_i + s \vec{r}_{ij})

where the summation is over all pairs of atoms in the chunk that are in contact (nonzero force), :math:`f_{ij,a}` is the force on atom :math:`i`
from atom :math:`j`, and :math:`\vec{r}_{ij}` is the displacement between
the two atoms.

The *stress* field is the summation of the kinetic and contact contributions.

The *fabric* field is defined as

.. math::

   \sum_{i,j} V_i r_{ij,a} r_{ij,b} \int_0^1 ds W(\vec{r}_\mathrm{chunk} - \vec{r}_i + s \vec{r}_{ij})

where :math:`V_i` is the volume of the atom in 3D and area in 2D.

TBD: The *momentum/grad* field is defined as

.. math::

   \sum_i V_i (p_{c,a} - m_i v_{i,a}) \grad_b W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

where :math:`p_{c,a}` is the :math:`a`-component of the *momentum* field and :math:`\grad_b W` is the :math:`b`-component of the of the gradient
of the kernel.

TBD: The *velocity/grad* field is defined as

.. math::

   \sum_i V_i (v_{c,a} - v_{i,a}) \grad_b W(\vec{r}_\mathrm{chunk} - \vec{r}_i)

where :math:`v_{c,a}` is the :math:`a`-component of the *velocity* field.

TBD: The *strain/rate* is defined as

.. math::

   \frac{1}{2} \left(\grad_{ab} v + \grad_{ba} v)

where :math:`\grad_{ab} v` is the :math:`ab` component of the *velocity/grad*
field.

----------

The optional *boundary/atom* and *boundary/fix* keywords turn on corrections
to the *stress* and *stress/contacts* fields using the construction from
:ref:`(Weinhart)`. At least one of these two keywords must be used to compute
the *boundary/force* vector.

The *boundary/atom* allows a group-ID to be specified. All atoms in this
group-ID will no longer contribute to any of the above summations, except
for *boundary/force*. In addition, these atoms will add an alternate contribution
to the *stress* and *stress/contacts* fields:

.. math::

   -\sum_{i,k} f_{ik,a} a_{ik,b} \int_0^1 ds W(\vec{r}_\mathrm{chunk} - \vec{r}_i + s \vec{a}_{ik})

where :math:`\vec{a}_{ik}` is the vector between the position of non-boundary
atom :math:`i` and the location of its contact with boundary atom :math:`k`.

The *boundary/fix* keyword applies analogous corrections from boundaries created
with instances of :doc:`fix wall/gran<fix_wall_gran>` or
:doc:`fix wall/gran/region<fix_wall_gran_region>`. An error will be triggered
if no such fixes are detected.

----------

The *file* or *append* keywords allow a filename to be specified.  If
*file* is used, then the filename is overwritten if it already exists.
If *append* is used, then the filename is appended to if it already
exists, or created if it does not exist.  Every :math:`N_\text{freq}`
timesteps, a section of chunk info will be written to a text file in the
following format.  A line with the timestep and number of chunks is
written.  Then one line per chunk is written, containing the chunk ID
:math:`(1-N_\text{chunk}),` an optional original ID value, optional
coordinate values for chunks that represent spatial bins, the number of
atoms in the chunk, and one or more calculated values.  More explanation
of the optional values is given below.  The number of values in each
line corresponds to the number of values specified in the fix ave/chunk
command.  The number of atoms and the value(s) are summed or average
quantities, as explained above.

The *overwrite* keyword will continuously overwrite the output file
with the latest output, so that it only contains one timestep worth of
output.  This option can only be used with the *ave running* setting.

The *format* keyword sets the numeric format of each value when it is
printed to a file via the *file* keyword.  Note that all values are
floating point quantities.  The default format is " %g".  You can specify
a higher precision if desired (e.g., " %20.16g").

The *title1* and *title2* and *title3* keywords allow specification of
the strings that will be printed as the first three lines of the output
file, assuming the *file* keyword was used.  LAMMPS uses default
values for each of these, so they do not need to be specified.

By default, these header lines are as follows:

.. parsed-literal::

   # Chunk-averaged data for fix ID and group name
   # Timestep Number-of-chunks
   # Chunk (OrigID) (Coord1) (Coord2) (Coord3) Ncount value1 value2 ...

In the first line, ID and name are replaced with the fix-ID and group
name.  The second line describes the two values that are printed at
the first of each section of output.  In the third line the values are
replaced with the appropriate value names (e.g., *v/x* or *stress/xy*).

The words in parenthesis only appear with corresponding columns if the
chunk style specified for the :doc:`compute chunk/atom
<compute_chunk_atom>` command supports them.  The OrigID column is
only used if the *compress* keyword was set to *yes* for the
:doc:`compute chunk/atom <compute_chunk_atom>` command.  This means
that the original chunk IDs (e.g., molecule IDs) will have been
compressed to remove chunk IDs with no atoms assigned to them.  Thus a
compressed chunk ID of 3 may correspond to an original chunk ID or
molecule ID of 415.  The OrigID column will list 415 for the third chunk.

The CoordN columns depends on the *binning* style.  For *bin/1d*,
*bin/2d*, and *bin/3d* styles the column values are the center point
of the bin in the corresponding dimension.  Just Coord1 is used for
*bin/1d*, Coord2 is added for *bin/2d*, Coord3 is added for *bin/3d*\ .

Note that if the value of the *units* keyword used in the
:doc:`compute chunk/atom command <compute_chunk_atom>` is *box* or
*lattice*, the coordinate values will be in distance :doc:`units <units>`.
If the value of the *units* keyword is *reduced*, the
coordinate values will be in unitless reduced units (0--1).

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.

This fix computes a global array of values which can be accessed by
various :doc:`output commands <Howto_output>`.  The values can only be
accessed on timesteps that are multiples of :math:`N_\text{freq}`, since
that is when averaging is performed.  The global array has # of rows =
the number of chunks :math:`N_\text{chunk}`, as calculated by the
specified :doc:`compute chunk/atom <compute_chunk_atom>` command.  The #
of columns is :math:`M+1+N_\text{values}`, where :math:`M \in
\{1,\dotsc,4\}`, depending on whether the optional columns for OrigID
and CoordN are used, as explained above. Note that a wildcard increases
the number of values output as columns. Following the optional
columns, the next column contains the count of atoms in the chunk, and
the remaining columns are the Nvalue quantities.  When the array is
accessed with a row :math:`I` that exceeds the current number of chunks,
than a 0.0 is returned by the fix instead of an error, since the number
of chunks can vary as a simulation runs depending on how that value is
computed by the compute chunk/atom command.

The array values calculated by this fix are treated as "intensive",
since they are typically already normalized by the count of atoms in
each chunk.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the GRANULAR package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute <compute>`, :doc:`fix ave/atom <fix_ave_atom>`,
:doc:`fix ave/histo <fix_ave_histo>`, :doc:`fix ave/time <fix_ave_time>`,
:doc:`variable <variable>`, :doc:`fix ave/correlate <fix_ave_correlate>`,
:doc:`fix ave/grid <fix_ave_grid>`

Default
"""""""

The option defaults are ave = one, no file
output, and title 1,2,3 = strings as described above.

----------

.. _fix_continuum_chunk_goldhirsch:

**(Goldhirsch)** Goldhirsch, Granular Matter, 12, 3, 239-252 (2010).

**(Weinhart)** Weinhart, Thornton, Luding, Bokhove, Granular Matter, 14, 2, 289-294 (2012).
