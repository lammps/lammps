.. index:: fix ave/chunk

fix ave/chunk command
=====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID ave/chunk Nevery Nrepeat Nfreq chunkID value1 value2 ... keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* ave/chunk = style name of this fix command
* Nevery = use input values every this many timesteps
* Nrepeat = # of times to use input values for calculating averages
* Nfreq = calculate averages every this many timesteps
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
* one or more input values can be listed
* value = vx, vy, vz, fx, fy, fz, density/mass, density/number, temp, c\_ID, c\_ID[I], f\_ID, f\_ID[I], v\_name
  
  .. parsed-literal::
  
       vx,vy,vz,fx,fy,fz = atom attribute (velocity, force component)
       density/number, density/mass = number or mass density
       temp = temperature
       c_ID = per-atom vector calculated by a compute with ID
       c_ID[I] = Ith column of per-atom array calculated by a compute with ID, I can include wildcard (see below)
       f_ID = per-atom vector calculated by a fix with ID
       f_ID[I] = Ith column of per-atom array calculated by a fix with ID, I can include wildcard (see below)
       v_name = per-atom vector calculated by an atom-style variable with name

* zero or more keyword/arg pairs may be appended
* keyword = *norm* or *ave* or *bias* or *adof* or *cdof* or *file* or *overwrite* or *title1* or *title2* or *title3*
  
  .. parsed-literal::
  
       *norm* arg = *all* or *sample* or *none* = how output on *Nfreq* steps is normalized
         all = output is sum of atoms across all *Nrepeat* samples, divided by atom count
         sample = output is sum of *Nrepeat* sample averages, divided by *Nrepeat*
         none = output is sum of *Nrepeat* sample sums, divided by *Nrepeat*
       *ave* args = *one* or *running* or *window M*
         one = output new average value every Nfreq steps
         running = output cumulative average of all previous Nfreq steps
         window M = output average of M most recent Nfreq steps
       *bias* arg = bias-ID
         bias-ID = ID of a temperature compute that removes a velocity bias for temperature calculation
       *adof* value = dof_per_atom
         dof_per_atom = define this many degrees-of-freedom per atom for temperature calculation
       *cdof* value = dof_per_chunk
         dof_per_chunk = define this many degrees-of-freedom per chunk for temperature calculation
       *file* arg = filename
         filename = file to write results to
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


.. parsed-literal::

   fix 1 all ave/chunk 10000 1 10000 binchunk c_myCentro title1 "My output values"
   fix 1 flow ave/chunk 100 10 1000 molchunk vx vz norm sample file vel.profile
   fix 1 flow ave/chunk 100 5 1000 binchunk density/mass ave running
   fix 1 flow ave/chunk 100 5 1000 binchunk density/mass ave running

**NOTE:**

If you are trying to replace a deprecated fix ave/spatial command
with the newer, more flexible fix ave/chunk and :doc:`compute chunk/atom <compute_chunk_atom>` commands, you simply need to split
the fix ave/spatial arguments across the two new commands.  For
example, this command:


.. parsed-literal::

   fix 1 flow ave/spatial 100 10 1000 y 0.0 1.0 vx vz norm sample file vel.profile

could be replaced by:


.. parsed-literal::

   compute cc1 flow chunk/atom bin/1d y 0.0 1.0
   fix 1 flow ave/chunk 100 10 1000 cc1 vx vz norm sample file vel.profile

Description
"""""""""""

Use one or more per-atom vectors as inputs every few timesteps, sum
the values over the atoms in each chunk at each timestep, then average
the per-chunk values over longer timescales.  The resulting chunk
averages can be used by other :doc:`output commands <Howto_output>` such
as :doc:`thermo\_style custom <thermo_style>`, and can also be written to
a file.

In LAMMPS, chunks are collections of atoms defined by a :doc:`compute chunk/atom <compute_chunk_atom>` command, which assigns each atom
to a single chunk (or no chunk).  The ID for this command is specified
as chunkID.  For example, a single chunk could be the atoms in a
molecule or atoms in a spatial bin.  See the :doc:`compute chunk/atom <compute_chunk_atom>` doc page and the :doc:`Howto chunk <Howto_chunk>` doc page for details of how chunks can be
defined and examples of how they can be used to measure properties of
a system.

Note that only atoms in the specified group contribute to the summing
and averaging calculations.  The :doc:`compute chunk/atom <compute_chunk_atom>` command defines its own group as
well as an optional region.  Atoms will have a chunk ID = 0, meaning
they belong to no chunk, if they are not in that group or region.
Thus you can specify the "all" group for this command if you simply
want to use the chunk definitions provided by chunkID.

Each specified per-atom value can be an atom attribute (position,
velocity, force component), a mass or number density, or the result of
a :doc:`compute <compute>` or :doc:`fix <fix>` or the evaluation of an
atom-style :doc:`variable <variable>`.  In the latter cases, the
compute, fix, or variable must produce a per-atom quantity, not a
global quantity.  Note that the :doc:`compute property/atom <compute_property_atom>` command provides access to
any attribute defined and stored by atoms.  If you wish to
time-average global quantities from a compute, fix, or variable, then
see the :doc:`fix ave/time <fix_ave_time>` command.

The per-atom values of each input vector are summed and averaged
independently of the per-atom values in other input vectors.

:doc:`Computes <compute>` that produce per-atom quantities are those
which have the word *atom* in their style name.  See the doc pages for
individual :doc:`fixes <fix>` to determine which ones produce per-atom
quantities.  :doc:`Variables <variable>` of style *atom* are the only
ones that can be used with this fix since all other styles of variable
produce global quantities.

Note that for values from a compute or fix, the bracketed index I can
be specified using a wildcard asterisk with the index to effectively
specify multiple values.  This takes the form "\*" or "\*n" or "n\*" or
"m\*n".  If N = the size of the vector (for *mode* = scalar) or the
number of columns in the array (for *mode* = vector), then an asterisk
with no numeric values means all indices from 1 to N.  A leading
asterisk means all indices from 1 to n (inclusive).  A trailing
asterisk means all indices from n to N (inclusive).  A middle asterisk
means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual columns of the array
had been listed one by one.  E.g. these 2 fix ave/chunk commands are
equivalent, since the :doc:`compute property/atom <compute_property_atom>` command creates, in this
case, a per-atom array with 3 columns:


.. parsed-literal::

   compute myAng all property/atom angmomx angmomy angmomz
   fix 1 all ave/chunk 100 1 100 cc1 c_myAng[\*] file tmp.angmom
   fix 2 all ave/chunk 100 1 100 cc1 c_myAng[1] c_myAng[2] c_myAng[3] file tmp.angmom

.. note::

   This fix works by creating an array of size *Nchunk* by Nvalues
   on each processor.  *Nchunk* is the number of chunks which is defined
   by the :doc:`compute chunk/atom <compute_chunk_atom>` command.
   Nvalues is the number of input values specified.  Each processor loops
   over its atoms, tallying its values to the appropriate chunk.  Then
   the entire array is summed across all processors.  This means that
   using a large number of chunks will incur an overhead in memory and
   computational cost (summing across processors), so be careful to
   define a reasonable number of chunks.


----------


The *Nevery*\ , *Nrepeat*\ , and *Nfreq* arguments specify on what
timesteps the input values will be accessed and contribute to the
average.  The final averaged quantities are generated on timesteps
that are a multiples of *Nfreq*\ .  The average is over *Nrepeat*
quantities, computed in the preceding portion of the simulation every
*Nevery* timesteps.  *Nfreq* must be a multiple of *Nevery* and
*Nevery* must be non-zero even if *Nrepeat* is 1.  Also, the timesteps
contributing to the average value cannot overlap, i.e. Nrepeat\*Nevery
can not exceed Nfreq.

For example, if Nevery=2, Nrepeat=6, and Nfreq=100, then values on
timesteps 90,92,94,96,98,100 will be used to compute the final average
on timestep 100.  Similarly for timesteps 190,192,194,196,198,200 on
timestep 200, etc.  If Nrepeat=1 and Nfreq = 100, then no time
averaging is done; values are simply generated on timesteps
100,200,etc.

Each input value can also be averaged over the atoms in each chunk.
The way the averaging is done across the *Nrepeat* timesteps to
produce output on the *Nfreq* timesteps, and across multiple *Nfreq*
outputs, is determined by the *norm* and *ave* keyword settings, as
discussed below.

.. note::

   To perform per-chunk averaging within a *Nfreq* time window, the
   number of chunks *Nchunk* defined by the :doc:`compute chunk/atom <compute_chunk_atom>` command must remain constant.  If
   the *ave* keyword is set to *running* or *window* then *Nchunk* must
   remain constant for the duration of the simulation.  This fix forces
   the chunk/atom compute specified by chunkID to hold *Nchunk* constant
   for the appropriate time windows, by not allowing it to re-calculate
   *Nchunk*\ , which can also affect how it assigns chunk IDs to atoms.
   This is particularly important to understand if the chunks defined by
   the :doc:`compute chunk/atom <compute_chunk_atom>` command are spatial
   bins.  If its *units* keyword is set to *box* or *lattice*\ , then the
   number of bins *Nchunk* and size of each bin will be fixed over the
   *Nfreq* time window, which can affect which atoms are discarded if the
   simulation box size changes.  If its *units* keyword is set to
   *reduced*\ , then the number of bins *Nchunk* will still be fixed, but
   the size of each bin can vary at each timestep if the simulation box
   size changes, e.g. for an NPT simulation.


----------


The atom attribute values (vx,vy,vz,fx,fy,fz) are self-explanatory.
As noted above, any other atom attributes can be used as input values
to this fix by using the :doc:`compute property/atom <compute_property_atom>` command and then specifying
an input value from that compute.

The *density/number* value means the number density is computed for
each chunk, i.e. number/volume.  The *density/mass* value means the
mass density is computed for each chunk, i.e. total-mass/volume.  The
output values are in units of 1/volume or density (mass/volume).  See
the :doc:`units <units>` command doc page for the definition of density
for each choice of units, e.g. gram/cm\^3.  If the chunks defined by
the :doc:`compute chunk/atom <compute_chunk_atom>` command are spatial
bins, the volume is the bin volume.  Otherwise it is the volume of the
entire simulation box.

The *temp* value means the temperature is computed for each chunk, by
the formula KE = DOF/2 k T, where KE = total kinetic energy of the
chunk of atoms (sum of 1/2 m v\^2), DOF = the total number of degrees
of freedom for all atoms in the chunk, k = Boltzmann constant, and T =
temperature.

The DOF is calculated as N\*adof + cdof, where N = number of atoms in
the chunk, adof = degrees of freedom per atom, and cdof = degrees of
freedom per chunk.  By default adof = 2 or 3 = dimensionality of
system, as set via the :doc:`dimension <dimension>` command, and cdof =
0.0.  This gives the usual formula for temperature.

Note that currently this temperature only includes translational
degrees of freedom for each atom.  No rotational degrees of freedom
are included for finite-size particles.  Also no degrees of freedom
are subtracted for any velocity bias or constraints that are applied,
such as :doc:`compute temp/partial <compute_temp_partial>`, or :doc:`fix shake <fix_shake>` or :doc:`fix rigid <fix_rigid>`.  This is because
those degrees of freedom (e.g. a constrained bond) could apply to sets
of atoms that are both included and excluded from a specific chunk,
and hence the concept is somewhat ill-defined.  In some cases, you can
use the *adof* and *cdof* keywords to adjust the calculated degrees of
freedom appropriately, as explained below.

Also note that a bias can be subtracted from atom velocities before
they are used in the above formula for KE, by using the *bias*
keyword.  This allows, for example, a thermal temperature to be
computed after removal of a flow velocity profile.

Note that the per-chunk temperature calculated by this fix and the
:doc:`compute temp/chunk <compute_temp_chunk>` command can be different.
The compute calculates the temperature for each chunk for a single
snapshot.  This fix can do that but can also time average those values
over many snapshots, or it can compute a temperature as if the atoms
in the chunk on different timesteps were collected together as one set
of atoms to calculate their temperature.  The compute allows the
center-of-mass velocity of each chunk to be subtracted before
calculating the temperature; this fix does not.

If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  If no bracketed integer is
appended, the per-atom vector calculated by the compute is used.  If a
bracketed integer is appended, the Ith column of the per-atom array
calculated by the compute is used.  Users can also write code for
their own compute styles and :doc:`add them to LAMMPS <Modify>`.
See the discussion above for how I can be specified with a wildcard
asterisk to effectively specify multiple values.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  If no bracketed integer is
appended, the per-atom vector calculated by the fix is used.  If a
bracketed integer is appended, the Ith column of the per-atom array
calculated by the fix is used.  Note that some fixes only produce
their values on certain timesteps, which must be compatible with
*Nevery*\ , else an error results.  Users can also write code for their
own fix styles and :doc:`add them to LAMMPS <Modify>`.  See the
discussion above for how I can be specified with a wildcard asterisk
to effectively specify multiple values.

If a value begins with "v\_", a variable name must follow which has
been previously defined in the input script.  Variables of style
*atom* can reference thermodynamic keywords and various per-atom
attributes, or invoke other computes, fixes, or variables when they
are evaluated, so this is a very general means of generating per-atom
quantities to average within chunks.


----------


Additional optional keywords also affect the operation of this fix
and its outputs.

The *norm* keyword affects how averaging is done for the per-chunk
values that are output every *Nfreq* timesteps.

It the *norm* setting is *all*\ , which is the default, a chunk value is
summed over all atoms in all *Nrepeat* samples, as is the count of
atoms in the chunk.  The averaged output value for the chunk on the
*Nfreq* timesteps is Total-sum / Total-count.  In other words it is an
average over atoms across the entire *Nfreq* timescale.  For the
*density/number* and *density/mass* values, the volume (bin volume or
system volume) used in the final normalization will be the volume at
the final *Nfreq* timestep.

If the *norm* setting is *sample*\ , the chunk value is summed over
atoms for each sample, as is the count, and an "average sample value"
is computed for each sample, i.e. Sample-sum / Sample-count.  The
output value for the chunk on the *Nfreq* timesteps is the average of
the *Nrepeat* "average sample values", i.e. the sum of *Nrepeat*
"average sample values" divided by *Nrepeat*\ .  In other words it is an
average of an average.  For the *density/number* and *density/mass*
values, the volume (bin volume or system volume) used in the
per-sample normalization will be the current volume at each sampling
step.

If the *norm* setting is *none*\ , a similar computation as for the
*sample* setting is done, except the individual "average sample
values" are "summed sample values".  A summed sample value is simply
the chunk value summed over atoms in the sample, without dividing by
the number of atoms in the sample.  The output value for the chunk on
the *Nfreq* timesteps is the average of the *Nrepeat* "summed sample
values", i.e. the sum of *Nrepeat* "summed sample values" divided by
*Nrepeat*\ .  For the *density/number* and *density/mass* values, the
volume (bin volume or system volume) used in the per-sample sum
normalization will be the current volume at each sampling step.

The *ave* keyword determines how the per-chunk values produced every
*Nfreq* steps are averaged with values produced on previous steps that
were multiples of *Nfreq*\ , before they are accessed by another output
command or written to a file.

If the *ave* setting is *one*\ , which is the default, then the chunk
values produced on timesteps that are multiples of *Nfreq* are
independent of each other; they are output as-is without further
averaging.

If the *ave* setting is *running*\ , then the chunk values produced on
timesteps that are multiples of *Nfreq* are summed and averaged in a
cumulative sense before being output.  Each output chunk value is thus
the average of the chunk value produced on that timestep with all
preceding values for the same chunk.  This running average begins when
the fix is defined; it can only be restarted by deleting the fix via
the :doc:`unfix <unfix>` command, or re-defining the fix by
re-specifying it.

If the *ave* setting is *window*\ , then the chunk values produced on
timesteps that are multiples of *Nfreq* are summed and averaged within
a moving "window" of time, so that the last M values for the same
chunk are used to produce the output.  E.g. if M = 3 and Nfreq = 1000,
then the output on step 10000 will be the average of the individual
chunk values on steps 8000,9000,10000.  Outputs on early steps will
average over less than M values if they are not available.

The *bias* keyword specifies the ID of a temperature compute that
removes a "bias" velocity from each atom, specified as *bias-ID*\ .  It
is only used when the *temp* value is calculated, to compute the
thermal temperature of each chunk after the translational kinetic
energy components have been altered in a prescribed way, e.g.  to
remove a flow velocity profile.  See the doc pages for individual
computes that calculate a temperature to see which ones implement a
bias.

The *adof* and *cdof* keywords define the values used in the degree of
freedom (DOF) formula described above for temperature calculation
for each chunk.  They are only used when the *temp* value is
calculated.  They can be used to calculate a more appropriate
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

The *file* keyword allows a filename to be specified.  Every *Nfreq*
timesteps, a section of chunk info will be written to a text file in
the following format.  A line with the timestep and number of chunks
is written.  Then one line per chunk is written, containing the chunk
ID (1-Nchunk), an optional original ID value, optional coordinate
values for chunks that represent spatial bins, the number of atoms in
the chunk, and one or more calculated values.  More explanation of the
optional values is given below.  The number of values in each line
corresponds to the number of values specified in the fix ave/chunk
command.  The number of atoms and the value(s) are summed or average
quantities, as explained above.

The *overwrite* keyword will continuously overwrite the output file
with the latest output, so that it only contains one timestep worth of
output.  This option can only be used with the *ave running* setting.

The *format* keyword sets the numeric format of each value when it is
printed to a file via the *file* keyword.  Note that all values are
floating point quantities.  The default format is %g.  You can specify
a higher precision if desired, e.g. %20.16g.

The *title1* and *title2* and *title3* keywords allow specification of
the strings that will be printed as the first 3 lines of the output
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
replaced with the appropriate value names, e.g. fx or c\_myCompute\ **2**\ .

The words in parenthesis only appear with corresponding columns if the
chunk style specified for the :doc:`compute chunk/atom <compute_chunk_atom>` command supports them.  The OrigID
column is only used if the *compress* keyword was set to *yes* for the
:doc:`compute chunk/atom <compute_chunk_atom>` command.  This means that
the original chunk IDs (e.g. molecule IDs) will have been compressed
to remove chunk IDs with no atoms assigned to them.  Thus a compressed
chunk ID of 3 may correspond to an original chunk ID or molecule ID of
415.  The OrigID column will list 415 for the 3rd chunk.

The CoordN columns only appear if a *binning* style was used in the
:doc:`compute chunk/atom <compute_chunk_atom>` command.  For *bin/1d*\ ,
*bin/2d*\ , and *bin/3d* styles the column values are the center point
of the bin in the corresponding dimension.  Just Coord1 is used for
*bin/1d*\ , Coord2 is added for *bin/2d*\ , Coord3 is added for *bin/3d*\ .
For *bin/sphere*\ , just Coord1 is used, and it is the radial
coordinate.  For *bin/cylinder*\ , Coord1 and Coord2 are used.  Coord1
is the radial coordinate (away from the cylinder axis), and coord2 is
the coordinate along the cylinder axis.

Note that if the value of the *units* keyword used in the :doc:`compute chunk/atom command <compute_chunk_atom>` is *box* or *lattice*\ , the
coordinate values will be in distance :doc:`units <units>`.  If the
value of the *units* keyword is *reduced*\ , the coordinate values will
be in unitless reduced units (0-1).  This is not true for the Coord1 value
of style *bin/sphere* or *bin/cylinder* which both represent radial
dimensions.  Those values are always in distance :doc:`units <units>`.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix\_modify <fix_modify>` options
are relevant to this fix.

This fix computes a global array of values which can be accessed by
various :doc:`output commands <Howto_output>`.  The values can only be
accessed on timesteps that are multiples of *Nfreq* since that is when
averaging is performed.  The global array has # of rows = the number
of chunks *Nchunk* as calculated by the specified :doc:`compute chunk/atom <compute_chunk_atom>` command.  The # of columns =
M+1+Nvalues, where M = 1 to 4, depending on whether the optional
columns for OrigID and CoordN are used, as explained above.  Following
the optional columns, the next column contains the count of atoms in
the chunk, and the remaining columns are the Nvalue quantities.  When
the array is accessed with a row I that exceeds the current number of
chunks, than a 0.0 is returned by the fix instead of an error, since
the number of chunks can vary as a simulation runs depending on how
that value is computed by the compute chunk/atom command.

The array values calculated by this fix are treated as "intensive",
since they are typically already normalized by the count of atoms in
each chunk.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute <compute>`, :doc:`fix ave/atom <fix_ave_atom>`, :doc:`fix ave/histo <fix_ave_histo>`, :doc:`fix ave/time <fix_ave_time>`,
:doc:`variable <variable>`, :doc:`fix ave/correlate <fix_ave_correlate>`

Default
"""""""

The option defaults are norm = all, ave = one, bias = none, no file output, and
title 1,2,3 = strings as described above.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
