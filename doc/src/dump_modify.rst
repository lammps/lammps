.. index:: dump_modify

dump_modify command
===================

:doc:`dump_modify <dump_image>` command for image/movie options
===============================================================

Syntax
""""""

.. code-block:: LAMMPS

   dump_modify dump-ID keyword values ...

* dump-ID = ID of dump to modify
* one or more keyword/value pairs may be appended

* these keywords apply to various dump styles
* keyword = *append* or *at* or *balance* or *buffer* or *delay* or *element* or *every* or *every/time* or *fileper* or *first* or *flush* or *format* or *header* or *image* or *label* or *maxfiles* or *nfile* or *pad* or *pbc* or *precision* or *region* or *refresh* or *scale* or *sfactor* or *skip* or *sort* or *tfactor* or *thermo* or *thresh* or *time* or *units* or *unwrap*

  .. parsed-literal::

       *append* arg = *yes* or *no*
       *at* arg = N
         N = index of frame written upon first dump
       *balance* arg = *yes* or *no*
       *buffer* arg = *yes* or *no*
       *colname* values =  ID string, or *default*
         string = new column header name
         ID = integer from 1 to N, or integer from -1 to -N, where N = # of quantities being output
              *or* a custom dump keyword or reference to compute, fix, property or variable.
       *delay* arg = Dstep
         Dstep = delay output until this timestep
       *element* args = E1 E2 ... EN, where N = # of atom types
         E1,...,EN = element name (e.g., C or Fe or Ga)
       *every* arg = N
         N = dump on timesteps which are a multiple of N
         N can be a variable (see below)
       *every/time* arg = Delta
         Delta = dump once every Delta interval of simulation time (time units)
         Delta can be a variable (see below)
       *fileper* arg = Np
         Np = write one file for every this many processors
       *first* arg = *yes* or *no*
       *flush* arg = *yes* or *no*
       *format* args = *line* string, *int* string, *float* string, ID string, or *none*
         string = C-style format string
         ID = integer from 1 to N, or integer from -1 to -N, where N = # of quantities being output
              *or* a custom dump keyword or reference to compute, fix, property or variable.
       *header* arg = *yes* or *no*
         *yes* to write the header
         *no* to not write the header
       *image* arg = *yes* or *no*
       *label* arg = string
         string = character string (e.g., BONDS) to use in header of dump local file
       *maxfiles* arg = Fmax
         Fmax = keep only the most recent *Fmax* snapshots (one snapshot per file)
       *nfile* arg = Nf
         Nf = write this many files, one from each of Nf processors
       *pad* arg = Nchar = # of characters to convert timestep to
       *pbc* arg = *yes* or *no* = remap atoms via periodic boundary conditions
       *precision* arg = power-of-10 value from 10 to 1000000
       *region* arg = region-ID or "none"
       *refresh* arg = c_ID = compute ID that supports a refresh operation
       *scale* arg = *yes* or *no*
       *sfactor* arg = coordinate scaling factor (> 0.0)
       *skip* arg = v_name
         v_name = variable with name which evaluates to non-zero (skip) or 0
       *sort* arg = *off* or *id* or N or -N
          off = no sorting of per-atom lines within a snapshot
          id = sort per-atom lines by atom ID
          N = sort per-atom lines in ascending order by the Nth column
          -N = sort per-atom lines in descending order by the Nth column
       *tfactor* arg = time scaling factor (> 0.0)
       *thermo* arg = *yes* or *no*
       *time* arg = *yes* or *no*
       *thresh* args = attribute operator value
         attribute = same attributes (x,fy,etotal,sxx,etc) used by dump custom style
         operator = "<" or "<=" or ">" or ">=" or "==" or "!=" or "\|\^"
         value = numeric value to compare to, or LAST
         these 3 args can be replaced by the word "none" to turn off thresholding
       *units* arg = *yes* or *no*
       *unwrap* arg = *yes* or *no*

* these keywords apply only to the *image* and *movie* :doc:`styles <dump_image>`
* keyword = *acolor* or *adiam* or *amap* or *backcolor* or *bcolor* or *bdiam* or *boxcolor* or *color* or *bitrate* or *framerate*

  .. parsed-literal::

       see the :doc:`dump image <dump_image>` doc page for details

* these keywords apply only to the */gz* and */zstd* dump styles
* keyword = *compression_level*

  .. parsed-literal::

       *compression_level* args = level
         level = integer specifying the compression level that should be used (see below for supported levels)

* these keywords apply only to the */zstd* dump styles
* keyword = *checksum*

  .. parsed-literal::

       *checksum* args = *yes* or *no* (add checksum at end of zst file)

Examples
""""""""

.. code-block:: LAMMPS

   dump_modify 1 format line "%d %d %20.15g %g %g" scale yes
   dump_modify 1 format float %20.15g scale yes
   dump_modify myDump image yes scale no flush yes
   dump_modify 1 region mySphere thresh x < 0.0 thresh fx >= 3.2
   dump_modify xtcdump precision 10000 sfactor 0.1
   dump_modify 1 every 1000 nfile 20
   dump_modify 1 every v_myVar

Description
"""""""""""

Modify the parameters of a previously defined dump command.  Not all
parameters are relevant to all dump styles.

As explained on the :doc:`dump <dump>` doc page, the *atom/mpiio*,
*custom/mpiio*, and *xyz/mpiio* dump styles are identical in command
syntax and in the format of the dump files they create, to the
corresponding styles without "mpiio", except the single dump file they
produce is written in parallel via the MPI-IO library.  Thus if a
dump_modify option below is valid for the *atom* style, it is also
valid for the *atom/mpiio* style, and similarly for the other styles
which allow for use of MPI-IO.

----------

Unless otherwise noted, the following keywords apply to all the
various dump styles, including the :doc:`dump image <dump_image>` and
:doc:`dump movie <dump_image>` styles.

----------

The *append* keyword applies to all dump styles except *cfg* and *xtc*
and *dcd*\ .  It also applies only to text output files, not to binary
or gzipped or image/movie files.  If specified as *yes*, then dump
snapshots are appended to the end of an existing dump file.  If
specified as *no*, then a new dump file will be created which will
overwrite an existing file with the same name.

----------

The *at* keyword only applies to the *netcdf* dump style.  It can only
be used if the *append yes* keyword is also used.  The *N* argument is
the index of which frame to append to.  A negative value can be
specified for *N*, which means a frame counted from the end of the
file.  The *at* keyword can only be used if the dump_modify command is
before the first command that causes dump snapshots to be output
(e.g., a :doc:`run <run>` or :doc:`minimize <minimize>` command).  Once the
dump file has been opened, this keyword has no further effect.

----------

The *buffer* keyword applies only to dump styles *atom*, *cfg*,
*custom*, *local*, and *xyz*\ .  It also applies only to text output
files, not to binary or gzipped files.  If specified as *yes*, which
is the default, then each processor writes its output into an internal
text buffer, which is then sent to the processor(s) which perform file
writes, and written by those processors(s) as one large chunk of text.
If specified as *no*, each processor sends its per-atom data in binary
format to the processor(s) which perform file wirtes, and those
processor(s) format and write it line by line into the output file.

The buffering mode is typically faster since each processor does the
relatively expensive task of formatting the output for its own atoms.
However it requires about twice the memory (per processor) for the
extra buffering.

----------

.. versionadded:: 4May2022

The *colname* keyword can be used to change the default header keyword
for dump styles: *atom*, *custom*, and *cfg* and their compressed, ADIOS,
and MPIIO variants.  The setting for *ID string* replaces the default
text with the provided string.  *ID* can be a positive integer when it
represents the column number counting from the left, a negative integer
when it represents the column number from the right (i.e. -1 is the last
column/keyword), or a custom dump keyword (or compute, fix, property, or
variable reference) and then it replaces the string for that specific
keyword. For *atom* dump styles only the keywords "id", "type", "x",
"y", "z", "ix", "iy", "iz" can be accessed via string regardless of
whether scaled or unwrapped coordinates were enabled or disabled, and
it always assumes 8 columns for indexing regardless of whether image
flags are enabled or not.  For dump style *cfg* only changes to the
"auxiliary" keywords (6th or later keyword) will become visible.

The *colname* keyword can be used multiple times. If multiple *colname*
settings refer to the same keyword, the last setting has precedence.  A
setting of *default* clears all previous settings, reverting all values
to their default names. Using the *scale* or *image* keyword will also
reset all header keywords to their default values.

----------

The *delay* keyword applies to all dump styles.  No snapshots will be
output until the specified *Dstep* timestep or later.  Specifying
*Dstep* < 0 is the same as turning off the delay setting.  This is a
way to turn off unwanted output early in a simulation, for example,
during an equilibration phase.

----------

The *element* keyword applies only to the dump *cfg*, *xyz*, and
*image* styles.  It associates element names (e.g., H, C, Fe) with
LAMMPS atom types.  See the list of element names at the bottom of
this page.

In the case of dump *cfg*, this allows the `AtomEye <atomeye_>`_
visualization package to read the dump file and render atoms with the
appropriate size and color.

In the case of dump *image*, the output images will follow the same
`AtomEye <atomeye_>`_ convention.  An element name is specified for each
atom type (1 to Ntype) in the simulation.  The same element name can
be given to multiple atom types.

In the case of *xyz* format dumps, there are no restrictions to what
label can be used as an element name.  Any white-space separated text
will be accepted.

.. _atomeye: http://li.mit.edu/Archive/Graphics/A/

----------

The *every* keyword can be used with any dump style except the *dcd*
and *xtc* styles.  It specifies that the output of dump snapshots will
now be performed on timesteps which are a multiple of a new :math:`N`
value, This overrides the dump frequency originally specified by the
:doc:`dump <dump>` command.

The *every* keyword can be specified in one of two ways.  It can be a
numeric value in which case it must be > 0.  Or it can be an
:doc:`equal-style variable <variable>`, which should be specified as
v_name, where name is the variable name.

In this case, the variable is evaluated at the beginning of a run to
determine the next timestep at which a dump snapshot will be written
out.  On that timestep the variable will be evaluated again to
determine the next timestep, etc.  Thus the variable should return
timestep values.  See the stagger() and logfreq() and stride() math
functions for :doc:`equal-style variables <variable>`, as examples of
useful functions to use in this context.  Other similar math functions
could easily be added as options for :doc:`equal-style variables
<variable>`.  Also see the next() function, which allows use of a
file-style variable which reads successive values from a file, each
time the variable is evaluated.  Used with the *every* keyword, if the
file contains a list of ascending timesteps, you can output snapshots
whenever you wish.

Note that when using the variable option with the *every* keyword, you
need to use the *first* option if you want an initial snapshot written
to the dump file.  The *every* keyword cannot be used with the dump
*dcd* style.

For example, the following commands will
write snapshots at timesteps 0,10,20,30,100,200,300,1000,2000,etc:

.. code-block:: LAMMPS

   variable        s equal logfreq(10,3,10)
   dump            1 all atom 100 tmp.dump
   dump_modify     1 every v_s first yes

The following commands would write snapshots at the timesteps listed
in file tmp.times:

.. code-block:: LAMMPS

   variable        f file tmp.times
   variable        s equal next(f)
   dump            1 all atom 100 tmp.dump
   dump_modify     1 every v_s

.. note::

   When using a file-style variable with the *every* keyword, the
   file of timesteps must list a first timestep that is beyond the
   current timestep (e.g., it cannot be 0).  And it must list one or more
   timesteps beyond the length of the run you perform.  This is because
   the dump command will generate an error if the next timestep it reads
   from the file is not a value greater than the current timestep.  Thus
   if you wanted output on steps 0,15,100 of a 100-timestep run, the file
   should contain the values 15,100,101 and you should also use the
   dump_modify first command.  Any final value > 100 could be used in
   place of 101.

----------

The *every/time* keyword can be used with any dump style except the
*dcd* and *xtc* styles.  It changes the frequency of dump snapshots
from being based on the current timestep to being determined by
elapsed simulation time, i.e. in time units of the :doc:`units
<units>` command, and specifies *Delta* for the interval between
snapshots.  This can be useful when the timestep size varies during a
simulation run, e.g. by use of the :doc:`fix dt/reset <fix_dt_reset>`
command.  The default is to perform output on timesteps which a
multiples of specified timestep value :math:`N`; see the *every*
keyword.

The *every/time* keyword can be used with any dump style except the
*dcd* and *xtc* styles.  It does two things.  It specifies that the
interval between dump snapshots will be set in simulation time
(i.e. in time units of the :doc:`units <units>` command).  This can be
useful when the timestep size varies during a simulation run (e.g., by
use of the :doc:`fix dt/reset <fix_dt_reset>` command).  The default is
to specify the interval in timesteps; see the *every* keyword.  The
*every/time* command also sets the interval value.

.. note::

   If you wish dump styles *atom*, *custom*, *local*, or *xyz* to
   include the simulation time as a field in the header portion of
   each snapshot, you also need to use the dump_modify *time* keyword
   with a setting of *yes*.  See its documentation below.

Note that since snapshots are output on simulation steps, each
snapshot will be written on the first timestep whose associated
simulation time is >= the exact snapshot time value.

As with the *every* option, the *Delta* value can be specified in one
of two ways.  It can be a numeric value in which case it must be >
0.0.  Or it can be an :doc:`equal-style variable <variable>`, which
should be specified as v_name, where name is the variable name.

In this case, the variable is evaluated at the beginning of a run to
determine the next simulation time at which a dump snapshot will be
written out.  On that timestep the variable will be evaluated again to
determine the next simulation time, etc.  Thus the variable should
return values in time units.  Note the current timestep or simulation
time can be used in an :doc:`equal-style variables <variable>` since
they are both thermodynamic keywords.  Also see the next() function,
which allows use of a file-style variable which reads successive
values from a file, each time the variable is evaluated.  Used with
the *every/time* keyword, if the file contains a list of ascending
simulation times, you can output snapshots whenever you wish.

Note that when using the variable option with the *every/time*
keyword, you need to use the *first* option if you want an initial
snapshot written to the dump file.  The *every/time* keyword cannot be
used with the dump *dcd* style.

For example, the following commands will write snapshots at successive
simulation times which grow by a factor of 1.5 with each interval.
The dt value used in the variable is to avoid a zero result when the
initial simulation time is 0.0.

.. code-block:: LAMMPS

   variable        increase equal 1.5*(time+dt)
   dump            1 all atom 100 tmp.dump
   dump_modify     1 every/time v_increase first yes

The following commands would write snapshots at the times listed in
file tmp.times:

.. code-block:: LAMMPS

   variable        f file tmp.times
   variable        s equal next(f)
   dump            1 all atom 100 tmp.dump
   dump_modify     1 every/time v_s

.. note::

   When using a file-style variable with the *every/time* keyword, the
   file of timesteps must list a first time that is beyond the time
   associated with the current timestep (e.g., it cannot be 0.0).  And
   it must list one or more times beyond the length of the run you
   perform.  This is because the dump command will generate an error
   if the next time it reads from the file is not a value greater than
   the current time.  Thus if you wanted output at times 0,15,100 of a
   run of length 100 in simulation time, the file should contain the
   values 15,100,101 and you should also use the dump_modify first
   command.  Any final value > 100 could be used in place of 101.

----------

The *first* keyword determines whether a dump snapshot is written on
the very first timestep after the dump command is invoked.  This will
always occur if the current timestep is a multiple of $N$, the
frequency specified in the :doc:`dump <dump>` command or
:doc:`dump_modify every <dump_modify>` command, including timestep 0.
It will also always occur if the current simulation time is a multiple
of *Delta*, the time interval specified in the :doc:`dump_modify
every/time <dump_modify>` command.

But if this is not the case, a dump snapshot will only be written if
the setting of this keyword is *yes*\ .  If it is *no*, which is the
default, then it will not be written.

Note that if the argument to the :doc:`dump_modify every
<dump_modify>` :doc:`dump_modify every/time <dump_modify>` commands is
a variable and not a numeric value, then specifying *first yes* is the
only way to write a dump snapshot on the first timestep after the dump
command is invoked.

----------

The *flush* keyword determines whether a flush operation is invoked
after a dump snapshot is written to the dump file.  A flush insures
the output in that file is current (no buffering by the OS), even if
LAMMPS halts before the simulation completes.  Flushes cannot be
performed with dump style *xtc*\ .

----------

The *format* keyword can be used to change the default numeric format output
by the text-based dump styles: *atom*, *local*, *custom*, *cfg*, and
*xyz* styles, and their MPIIO variants. Only the *line* or *none*
options can be used with the *atom* and *xyz* styles.

All the specified format strings are C-style formats, such as used by
the C/C++ printf() command.  The *line* keyword takes a single
argument which is the format string for an entire line of output for
each atom (do not include a trailing "\n"), with :math:`N` fields, which you
must enclose in quotes if there is more than one field.  The *int* and
*float* keywords take a single format argument and are applied to all
integer or floating-point quantities output.  The setting for *M string*
also takes a single format argument which is used for the :math:`M`\ th
value output in each line (e.g., the fifth column is output in high
precision by "format 5 %20.15g").

.. note::

   When using the *line* keyword for the *cfg* style, the first two
   fields (atom ID and type) are not actually written into the CFG file,
   however you must include formats for them in the format string.

The *format* keyword can be used multiple times.  The precedence is
that for each value in a line of output, the *M* format (if specified)
is used, else the *int* or *float* setting (if specified) is used,
else the *line* setting (if specified) for that value is used, else
the default setting is used.  A setting of *none* clears all previous
settings, reverting all values to their default format.

.. note::

   Atom and molecule IDs are stored internally as 4-byte or 8-byte
   signed integers, depending on how LAMMPS was compiled.  When
   specifying the *format int* option you can use a "%d"-style format
   identifier in the format string and LAMMPS will convert this to the
   corresponding 8-byte form if it is needed when outputting those
   values.  However, when specifying the *line* option or *format M
   string* option for those values, you should specify a format string
   appropriate for an 8-byte signed integer (e.g., one with "%ld") if
   LAMMPS was compiled with the -DLAMMPS_BIGBIG option for 8-byte IDs.

.. note::

   Any value written to a text-based dump file that is a per-atom
   quantity calculated by a :doc:`compute <compute>` or :doc:`fix <fix>` is
   stored internally as a floating-point value.  If the value is actually
   an integer and you wish it to appear in the text dump file as a
   (large) integer, then you need to use an appropriate format.  For
   example, these commands:

.. code-block:: LAMMPS

   compute     1 all property/local batom1 batom2
   dump        1 all local 100 tmp.bonds index c_1[1] c_1[2]
   dump_modify 1 format line "%d %0.0f %0.0f"

will output the two atom IDs for atoms in each bond as integers.  If
the dump_modify command were omitted, they would appear as
floating-point values, assuming they were large integers (more than six
digits).  The "index" keyword should use the "%d" format since it is
not generated by a compute or fix, and is stored internally as an
integer.

----------

The *fileper* keyword is documented below with the *nfile* keyword.

----------

The *header* keyword toggles whether the dump file will include a
header.  Excluding a header will reduce the size of the dump file for
data produced by :doc:`pair tracker <pair_tracker>` or
:doc:`bpm bond styles <Howto_bpm>` which may not require the
information typically written to the header.

----------

The *image* keyword applies only to the dump *atom* style.  If the
image value is *yes*, three flags are appended to each atom's coords which
are the absolute box image of the atom in each dimension.  For
example, an :math:`x` image flag of :math:`-2` with a normalized coord of 0.5
means the atom is in the center of the box, but has passed through the box
boundary twice and is really two box lengths to the left of its
current coordinate.  Note that for dump style *custom* these various
values can be printed in the dump file by using the appropriate atom
attributes in the dump command itself.
Using this keyword will reset all custom header names set with
*dump_modify colname* to their respective default values.

----------

The *label* keyword applies only to the dump *local* style.
When it writes local information, such as bond or angle topology
to a dump file, it will use the specified *label* to format the header.
By default this includes two lines:

.. parsed-literal::

   ITEM: NUMBER OF ENTRIES
   ITEM: ENTRIES ...

The word "ENTRIES" will be replaced with the string specified
(e.g., BONDS or ANGLES).

----------

The *maxfiles* keyword can only be used when a '\*' wildcard is
included in the dump file name (i.e., when writing a new file(s) for
each snapshot).  The specified *Fmax* is how many snapshots will be
kept.  Once this number is reached, the file(s) containing the oldest
snapshot is deleted before a new dump file is written.  If the
specified :math:`\text{Fmax} \le 0`, then all files are retained.

This can be useful for debugging, especially if you do not know on what
timestep something bad will happen (e.g., when LAMMPS will exit with an
error).  You can dump every time step and limit the number of dump
files produced, even if you run for thousands of steps.

----------

The *nfile* or *fileper* keywords can be used in conjunction with the
"%" wildcard character in the specified dump file name, for all dump
styles except the *dcd*, *image*, *movie*, *xtc*, and *xyz* styles
(for which "%" is not allowed).  As explained on the :doc:`dump <dump>`
command doc page, the "%" character causes the dump file to be written
in pieces, one piece for each of :math:`P` processors.  By default, :math:`P`
is the number of processors the simulation is running on.  The *nfile* or
*fileper* keyword can be used to set :math:`P` to a smaller value, which can
be more efficient when running on a large number of processors.

The *nfile* keyword sets :math:`P` to the specified :math:`N_f` value.
For example, if :math:`N_f = 4`, and the simulation is running on 100
processors, four files will be written by processors 0, 25, 50, and 75.
Each will collect information from itself and the next 24 processors and write
it to a dump file.

For the *fileper* keyword, the specified value of :math:`N_p` means write one
file for every :math:`N_p` processors.  For example, if :math:`N_p = 4`,
every fourth processor (0, 4, 8, 12, etc.) will collect information from itself
and the next three processors and write it to a dump file.

----------

The *pad* keyword only applies when the dump filename is specified
with a wildcard "\*" character which becomes the timestep.  If *pad* is
0, which is the default, the timestep is converted into a string of
unpadded length (e.g., 100 or 12000 or 2000000).  When *pad* is
specified with *Nchar* :math:`>` 0, the string is padded with leading zeroes
so they are all the same length = *Nchar*\ .  For example, pad 7 would
yield 0000100, 0012000, 2000000.  This can be useful so that
post-processing programs can easily read the files in ascending
timestep order.

----------

The *pbc* keyword applies to all the dump styles.  As explained on the
:doc:`dump <dump>` doc page, atom coordinates in a dump file may be
slightly outside the simulation box.  This is because periodic
boundary conditions are enforced only on timesteps when neighbor lists
are rebuilt, which will not typically coincide with the timesteps dump
snapshots are written.  If the setting of this keyword is set to
*yes*, then all atoms will be remapped to the periodic box before the
snapshot is written, then restored to their original position.  If it
is set to *no* they will not be.  The *no* setting is the default
because it requires no extra computation.

----------

The *precision* keyword only applies to the dump *xtc* style.  A
specified value of :math:`N` means that coordinates are stored to :math:`1/N`
nanometer accuracy (e.g., for :math:`N = 1000`, the coordinates are written to
:math:`1/1000` nanometer accuracy).

----------

The *refresh* keyword only applies to the dump *custom*, *cfg*,
*image*, and *movie* styles.  It allows an "incremental" dump file to
be written, by refreshing a compute that is used as a threshold for
determining which atoms are included in a dump snapshot.  The
specified *c_ID* gives the ID of the compute.  It is prefixed by "c\_"
to indicate a compute, which is the only current option.  At some
point, other options may be added (e.g., fixes or variables).

.. note::

   This keyword can only be specified once for a dump.  Refreshes
   of multiple computes cannot yet be performed.

The definition and motivation of an incremental dump file is as
follows.  Instead of outputting all atoms at each snapshot (with some
associated values), you may only wish to output the subset of atoms
with a value that has changed in some way compared to the value the
last time that atom was output.  In some scenarios this can result in
a dramatically smaller dump file.  If desired, by post-processing the
sequence of snapshots, the values for all atoms at all timesteps can
be inferred.

A concrete example is a simulation of atom diffusion in a solid,
represented as atoms on a lattice.  Diffusive hops are rare.  Imagine
that when a hop occurs an atom moves more than a distance *Dhop*\ .  For
any snapshot we only want to output atoms that have hopped since the
last snapshot.  This can be accomplished with something the following
commands:

.. code-block:: LAMMPS

   variable        Dhop equal 0.6
   variable        check atom "c_dsp[4] > v_Dhop"
   compute         dsp all displace/atom refresh check
   dump            1 all custom 20 tmp.dump id type x y z
   dump_modify     1 append yes thresh c_dsp[4] > ${Dhop} refresh c_dsp

The :doc:`compute displace/atom <compute_displace_atom>` command
calculates the displacement of each atom from its reference position.
The "4" index is the scalar displacement; 1, 2, and 3 are the :math:`xyz`
components of the displacement.  The :doc:`dump_modify thresh <dump_modify>`
command will cause only atoms that have displaced more than
:math:`0.6~\mathrm{\mathring A}` to be output on a given snapshot (assuming
metal units).  However, note that when an atom is output, we also need to
update the reference position for that atom to its new coordinates.  So that it
will not be output in every snapshot thereafter.  That reference position is
stored by :doc:`compute displace/atom <compute_displace_atom>`.  So the
dump_modify *refresh* option triggers a call to compute displace/atom at the
end of every dump to perform that update.  The *refresh check* option
shown as part of the :doc:`compute displace/atom <compute_displace_atom>`
command enables the compute to respond to the call from the dump command, and
update the appropriate reference positions.  This is done be defining an
:doc:`atom-style variable <variable>`, *check* in this example, which
calculates a Boolean value (0 or 1) for each atom, based on the same
criterion used by dump_modify thresh.

See the :doc:`compute displace/atom <compute_displace_atom>` command for
more details, including an example of how to produce output that
includes an initial snapshot with the reference position of all atoms.

Note that only computes with a *refresh* option will work with
dump_modify refresh.  See individual compute doc pages for details.
Currently, only compute displace/atom supports this option.  Others
may be added at some point.  If you use a compute that does not support
refresh operations, LAMMPS will not complain; dump_modify refresh will
simply do nothing.

----------

The *region* keyword only applies to the dump *custom*, *cfg*,
*image*, and *movie* styles.  If specified, only atoms in the region
will be written to the dump file or included in the image/movie.  Only
one region can be applied as a filter (the last one specified).  See
the :doc:`region <region>` command for more details.  Note that a region
can be defined as the "inside" or "outside" of a geometric shape, and
it can be the "union" or "intersection" of a series of simpler
regions.

----------

The *scale* keyword applies only to the dump *atom* style.  A scale
value of *yes* means atom coords are written in normalized units from
0.0 to 1.0 in each box dimension.  If the simulation box is triclinic
(tilted), then all atom coords will still be between 0.0 and 1.0.  A
value of *no* means they are written in absolute distance units
(e.g., :math:`\mathrm{\mathring A}` or :math:`\sigma`).
Using this keyword will reset all custom header names set with
*dump_modify colname* to their respective default values.

----------

The *sfactor* and *tfactor* keywords only apply to the dump *xtc*
style.  They allow customization of the unit conversion factors used
when writing to XTC files.  By default, they are initialized for
whatever :doc:`units <units>` style is being used, to write out
coordinates in nanometers and time in picoseconds.  For example, for *real*
units, LAMMPS defines *sfactor* = 0.1 and *tfactor* = 0.001, since the
:math:`\mathrm{\mathring A}` and fs used by *real* units are 0.1 nm and
0.001 ps, respectively.  If you are using a units system with distance and time
units far from nm and ps, you may wish to write XTC files with
different units, since the compression algorithm used in XTC files is
most effective when the typical magnitude of position data is between
10.0 and 0.1.

----------

.. versionadded:: 15Sep2022

The *skip* keyword can be used with all dump styles.  It allows a dump
snapshot to be skipped (not written to the dump file), if a condition
is met.  The condition is computed by an :doc:`equal-style variable
<variable>`, which should be specified as v_name, where name is the
variable name.  If the variable evaluation returns a non-zero value,
then the dump snapshot is skipped.  If it returns zero, the dump
proceeds as usual.  Note that :doc:`equal-style variable <variable>`
can contain Boolean operators which effectively evaluate as a true
(non-zero) or false (zero) result.

The *skip* keyword can be useful for debugging purposes, e.g. to dump
only on a particular timestep.  Or to limit output to conditions of
interest, e.g. only when the force on some atom exceeds a threshold
value.

----------

The *sort* keyword determines whether lines of per-atom output in a
snapshot are sorted or not.  A sort value of *off* means they will
typically be written in indeterminate order, either in serial or
parallel.  This is the case even in serial if the :doc:`atom_modify sort
<atom_modify>` option is turned on, which it is by default, to improve
performance.  A sort value of *id* means sort the output by atom ID.  A
sort value of :math:`N` or :math:`-N` means sort the output by the value
in the :math:`N`\ th column of per-atom info in either ascending or
descending order.

The dump *local* style cannot be sorted by atom ID, since there are
typically multiple lines of output per atom.  Some dump styles, such
as *dcd* and *xtc*, require sorting by atom ID to format the output
file correctly.  If multiple processors are writing the dump file, via
the "%" wildcard in the dump filename and the *nfile* or *fileper*
keywords are set to non-default values (i.e., the number of dump file
pieces is not equal to the number of procs), then sorting cannot be
performed.

In a parallel run, the per-processor dump file pieces can have
significant imbalance in number of lines of per-atom info. The *balance*
keyword determines whether the number of lines in each processor
snapshot are balanced to be nearly the same. A balance value of *no*
means no balancing will be done, while *yes* means balancing will be
performed. This balancing preserves dump sorting order. For a serial
run, this option is ignored since the output is already balanced.

.. note::

   Unless it is required by the dump style, sorting dump file
   output requires extra overhead in terms of CPU and communication cost,
   as well as memory, versus unsorted output.

----------

The *thermo* keyword only applies the dump styles *netcdf* and *yaml*.
It triggers writing of :doc:`thermo <thermo>` information to the dump file
alongside per-atom data.  The values included in the dump file are
identical to the values specified by :doc:`thermo_style <thermo_style>`.

----------

The *thresh* keyword only applies to the dump *custom*, *cfg*,
*image*, and *movie* styles.  Multiple thresholds can be specified.
Specifying *none* turns off all threshold criteria.  If thresholds are
specified, only atoms whose attributes meet all the threshold criteria
are written to the dump file or included in the image.  The possible
attributes that can be tested for are the same as those that can be
specified in the :doc:`dump custom <dump>` command, with the exception
of the *element* attribute, since it is not a numeric value.  Note
that a different attributes can be used than those output by the
:doc:`dump custom <dump>` command.  For example, you can output the
coordinates and stress of atoms whose energy is above some threshold.

If an atom-style variable is used as the attribute, then it can
produce continuous numeric values or effective Boolean 0/1 values,
which may be useful for the comparison operator.  Boolean values can
be generated by variable formulas that use comparison or Boolean math
operators or special functions like gmask() and rmask() and grmask().
See the :doc:`variable <variable>` command page for details.

The specified value must be a simple numeric value or the word LAST.
If LAST is used, it refers to the value of the attribute the last time
the dump command was invoked to produce a snapshot.  This is a way to
only dump atoms whose attribute has changed (or not changed).
Three examples follow.

.. code-block:: LAMMPS

   dump_modify ... thresh ix != LAST

This will dump atoms which have crossed the periodic :math:`x` boundary of the
simulation box since the last dump.  (Note that atoms that crossed
once and then crossed back between the two dump timesteps would not be
included.)

.. code-block:: LAMMPS

   region foo sphere 10 20 10 15
   variable inregion atom rmask(foo)
   dump_modify ... thresh v_inregion |^ LAST

This will dump atoms which crossed the boundary of the spherical
region since the last dump.

.. code-block:: LAMMPS

   variable charge atom "(q > 0.5) || (q < -0.5)"
   dump_modify ... thresh v_charge |^ LAST

This will dump atoms whose charge has changed from an absolute value
less than :math:`\frac12` to greater than :math:`\frac12` (or vice versa) since the last dump (e.g., due to reactions and subsequent charge equilibration in a
reactive force field).

The choice of operators listed above are the usual comparison
operators.  The XOR operation (exclusive or) is also included as "\|\^".
In this context, XOR means that if either the attribute or value is
0.0 and the other is non-zero, then the result is "true" and the
threshold criterion is met.  Otherwise it is not met.

----------

The *time* keyword only applies to the dump *atom*, *custom*, *local*,
and *xyz* styles (and their COMPRESS package versions *atom/gz*,
*custom/gz* and *local/gz*\ ).  For the first three styles, if set to
*yes*, each frame will will contain two extra lines before the "ITEM:
TIMESTEP" entry:

.. parsed-literal::

   ITEM: TIME
   \<elapsed time\>

For the *xyz* style, the simulation time is included on the same line
as the timestep value.

This will output the current elapsed simulation time in current
time units equivalent to the :doc:`thermo keyword <thermo_style>` *time*\ .
This is to simplify post-processing of trajectories using a variable time
step (e.g., when using :doc:`fix dt/reset <fix_dt_reset>`).
The default setting is *no*\ .

----------

The *units* keyword only applies to the dump *atom*, *custom*, and
*local* styles (and their COMPRESS package versions *atom/gz*,
*custom/gz* and *local/gz*\ ). If set to *yes*, each individual dump
file will contain two extra lines at the very beginning with:

.. parsed-literal::

   ITEM: UNITS
   \<units style\>

This will output the current selected :doc:`units <units>` style
to the dump file and thus allows visualization and post-processing
tools to determine the choice of units of the data in the dump file.
The default setting is *no*\ .

----------

The *unwrap* keyword only applies to the dump *dcd* and *xtc* styles.
If set to *yes*, coordinates will be written "unwrapped" by the image
flags for each atom.  Unwrapped means that if the atom has passed through
a periodic boundary one or more times, the value is printed for what
the coordinate would be if it had not been wrapped back into the
periodic box.  Note that these coordinates may thus be far outside the
box size stored with the snapshot.

----------

The COMPRESS package offers both GZ and Zstd compression variants of
styles atom, custom, local, cfg, and xyz. When using these styles the
compression level can be controlled by the :code:`compression_level`
keyword. File names with these styles have to end in either
:code:`.gz` or :code:`.zst`.

GZ supports compression levels from :math:`-1` (default), 0 (no compression),
and 1 to 9, 9 being the best compression. The COMPRESS :code:`/gz` styles use 9
as default compression level.

Zstd offers a wider range of compression levels, including negative
levels that sacrifice compression for performance. 0 is the default,
positive levels are 1 to 22, with 22 being the most expensive
compression. Zstd promises higher compression/decompression speeds for
similar compression ratios. For more details see
`http://facebook.github.io/zstd/`.

In addition, Zstd compressed files can include a checksum of the
entire contents. The Zstd enabled dump styles enable this feature by
default and it can be disabled with the :code:`checksum` keyword.

----------

Restrictions
""""""""""""

Not all *dump_modify* options can be applied to all dump styles.
Details are in the discussions of the individual options.

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump image <dump_image>`, :doc:`undump <undump>`

Default
"""""""

The option defaults are

* append = no
* balance = no
* buffer = yes for dump styles *atom*, *custom*, *loca*, and *xyz*
* element = "C" for every atom type
* every = whatever it was set to via the :doc:`dump <dump>` command
* fileper = # of processors
* first = no
* flush = yes
* format = %d and %g for each integer or floating point value
* image = no
* label = ENTRIES
* maxfiles = -1
* nfile = 1
* pad = 0
* pbc = no
* precision = 1000
* region = none
* scale = yes
* sort = off for dump styles *atom*, *custom*, *cfg*, and *local*
* sort = id for dump styles *dcd*, *xtc*, and *xyz*
* thresh = none
* units = no
* unwrap = no

* compression_level = 9 (gz variants)
* compression_level = 0 (zstd variants)
* checksum = yes (zstd variants)

