.. index:: dump\_modify

dump\_modify command
====================

Syntax
""""""


.. parsed-literal::

   dump_modify dump-ID keyword values ...

* dump-ID = ID of dump to modify
* one or more keyword/value pairs may be appended
* these keywords apply to various dump styles
* keyword = *append* or *at* or *buffer* or *delay* or *element* or *every* or *fileper* or *first* or *flush* or *format* or *image* or *label* or *maxfiles* or *nfile* or *pad* or *pbc* or *precision* or *region* or *refresh* or *scale* or *sfactor* or *sort* or *tfactor* or *thermo* or *thresh* or *time* or *units* or *unwrap*
  
  .. parsed-literal::
  
       *append* arg = *yes* or *no*
       *at* arg = N
         N = index of frame written upon first dump
       *buffer* arg = *yes* or *no*
       *delay* arg = Dstep
         Dstep = delay output until this timestep
       *element* args = E1 E2 ... EN, where N = # of atom types
         E1,...,EN = element name, e.g. C or Fe or Ga
       *every* arg = N
         N = dump every this many timesteps
         N can be a variable (see below)
       *fileper* arg = Np
         Np = write one file for every this many processors
       *first* arg = *yes* or *no*
       *flush* arg = *yes* or *no*
       *format* args = *line* string, *int* string, *float* string, M string, or *none*
         string = C-style format string
         M = integer from 1 to N, where N = # of per-atom quantities being output
       *image* arg = *yes* or *no*
       *label* arg = string
         string = character string (e.g. BONDS) to use in header of dump local file
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
  
       *acolor* args = type color
         type = atom type or range of types (see below)
         color = name of color or color1/color2/...
       *adiam* args = type diam
         type = atom type or range of types (see below)
         diam = diameter of atoms of that type (distance units)
       *amap* args = lo hi style delta N entry1 entry2 ... entryN
         lo = number or *min* = lower bound of range of color map
         hi = number or *max* = upper bound of range of color map
         style = 2 letters = "c" or "d" or "s" plus "a" or "f"
           "c" for continuous
           "d" for discrete
           "s" for sequential
           "a" for absolute
           "f" for fractional
         delta = binsize (only used for style "s", otherwise ignored)
           binsize = range is divided into bins of this width
         N = # of subsequent entries
         entry = value color (for continuous style)
           value = number or *min* or *max* = single value within range
           color = name of color used for that value
         entry = lo hi color (for discrete style)
           lo/hi = number or *min* or *max* = lower/upper bound of subset of range
           color = name of color used for that subset of values
         entry = color (for sequential style)
           color = name of color used for a bin of values
       *backcolor* arg = color
         color = name of color for background
       *bcolor* args = type color
         type = bond type or range of types (see below)
         color = name of color or color1/color2/...
       *bdiam* args = type diam
         type = bond type or range of types (see below)
         diam = diameter of bonds of that type (distance units)
       *boxcolor* arg = color
         color = name of color for simulation box lines and processor sub-domain lines
       *color* args = name R G B
         name = name of color
         R,G,B = red/green/blue numeric values from 0.0 to 1.0
       *bitrate* arg = rate
         rate = target bitrate for movie in kbps
       *framerate* arg = fps
         fps = frames per second for movie



Examples
""""""""


.. parsed-literal::

   dump_modify 1 format line "%d %d %20.15g %g %g" scale yes
   dump_modify 1 format float %20.15g scale yes
   dump_modify myDump image yes scale no flush yes
   dump_modify 1 region mySphere thresh x < 0.0 thresh epair >= 3.2
   dump_modify xtcdump precision 10000 sfactor 0.1
   dump_modify 1 every 1000 nfile 20
   dump_modify 1 every v_myVar
   dump_modify 1 amap min max cf 0.0 3 min green 0.5 yellow max blue boxcolor red

Description
"""""""""""

Modify the parameters of a previously defined dump command.  Not all
parameters are relevant to all dump styles.

As explained on the :doc:`dump <dump>` doc page, the *atom/mpiio*\ ,
*custom/mpiio*\ , and *xyz/mpiio* dump styles are identical in command
syntax and in the format of the dump files they create, to the
corresponding styles without "mpiio", except the single dump file they
produce is written in parallel via the MPI-IO library.  Thus if a
dump\_modify option below is valid for the *atom* style, it is also
valid for the *atom/mpiio* style, and similarly for the other styles
which allow for use of MPI-IO.


----------


These keywords apply to various dump styles, including the :doc:`dump image <dump_image>` and :doc:`dump movie <dump_image>` styles.  The
description gives details.


----------


The *append* keyword applies to all dump styles except *cfg* and *xtc*
and *dcd*\ .  It also applies only to text output files, not to binary
or gzipped or image/movie files.  If specified as *yes*\ , then dump
snapshots are appended to the end of an existing dump file.  If
specified as *no*\ , then a new dump file will be created which will
overwrite an existing file with the same name.


----------


The *at* keyword only applies to the *netcdf* dump style.  It can only
be used if the *append yes* keyword is also used.  The *N* argument is
the index of which frame to append to.  A negative value can be
specified for *N*\ , which means a frame counted from the end of the
file.  The *at* keyword can only be used if the dump\_modify command is
before the first command that causes dump snapshots to be output,
e.g. a :doc:`run <run>` or :doc:`minimize <minimize>` command.  Once the
dump file has been opened, this keyword has no further effect.


----------


The *buffer* keyword applies only to dump styles *atom*\ , *cfg*\ ,
*custom*\ , *local*\ , and *xyz*\ .  It also applies only to text output
files, not to binary or gzipped files.  If specified as *yes*\ , which
is the default, then each processor writes its output into an internal
text buffer, which is then sent to the processor(s) which perform file
writes, and written by those processors(s) as one large chunk of text.
If specified as *no*\ , each processor sends its per-atom data in binary
format to the processor(s) which perform file wirtes, and those
processor(s) format and write it line by line into the output file.

The buffering mode is typically faster since each processor does the
relatively expensive task of formatting the output for its own atoms.
However it requires about twice the memory (per processor) for the
extra buffering.


----------


The *delay* keyword applies to all dump styles.  No snapshots will be
output until the specified *Dstep* timestep or later.  Specifying
*Dstep* < 0 is the same as turning off the delay setting.  This is a
way to turn off unwanted output early in a simulation, for example,
during an equilibration phase.


----------


The *element* keyword applies only to the dump *cfg*\ , *xyz*\ , and
*image* styles.  It associates element names (e.g. H, C, Fe) with
LAMMPS atom types.  See the list of element names at the bottom of
this page.

In the case of dump *cfg*\ , this allows the `AtomEye <atomeye_>`_
visualization package to read the dump file and render atoms with the
appropriate size and color.

In the case of dump *image*\ , the output images will follow the same
`AtomEye <atomeye_>`_ convention.  An element name is specified for each
atom type (1 to Ntype) in the simulation.  The same element name can
be given to multiple atom types.

In the case of *xyz* format dumps, there are no restrictions to what
label can be used as an element name.  Any white-space separated text
will be accepted.

.. _atomeye: http://mt.seas.upenn.edu/Archive/Graphics/A




----------


The *every* keyword changes the dump frequency originally specified by
the :doc:`dump <dump>` command to a new value.  The every keyword can be
specified in one of two ways.  It can be a numeric value in which case
it must be > 0.  Or it can be an :doc:`equal-style variable <variable>`,
which should be specified as v\_name, where name is the variable name.

In this case, the variable is evaluated at the beginning of a run to
determine the next timestep at which a dump snapshot will be written
out.  On that timestep the variable will be evaluated again to
determine the next timestep, etc.  Thus the variable should return
timestep values.  See the stagger() and logfreq() and stride() math
functions for :doc:`equal-style variables <variable>`, as examples of
useful functions to use in this context.  Other similar math functions
could easily be added as options for :doc:`equal-style variables <variable>`.  Also see the next() function, which allows
use of a file-style variable which reads successive values from a
file, each time the variable is evaluated.  Used with the *every*
keyword, if the file contains a list of ascending timesteps, you can
output snapshots whenever you wish.

Note that when using the variable option with the *every* keyword, you
need to use the *first* option if you want an initial snapshot written
to the dump file.  The *every* keyword cannot be used with the dump
*dcd* style.

For example, the following commands will
write snapshots at timesteps 0,10,20,30,100,200,300,1000,2000,etc:


.. parsed-literal::

   variable        s equal logfreq(10,3,10)
   dump            1 all atom 100 tmp.dump
   dump_modify     1 every v_s first yes

The following commands would write snapshots at the timesteps listed
in file tmp.times:


.. parsed-literal::

   variable        f file tmp.times
   variable        s equal next(f)
   dump            1 all atom 100 tmp.dump
   dump_modify     1 every v_s

.. note::

   When using a file-style variable with the *every* keyword, the
   file of timesteps must list a first timestep that is beyond the
   current timestep (e.g. it cannot be 0).  And it must list one or more
   timesteps beyond the length of the run you perform.  This is because
   the dump command will generate an error if the next timestep it reads
   from the file is not a value greater than the current timestep.  Thus
   if you wanted output on steps 0,15,100 of a 100-timestep run, the file
   should contain the values 15,100,101 and you should also use the
   dump\_modify first command.  Any final value > 100 could be used in
   place of 101.


----------


The *first* keyword determines whether a dump snapshot is written on
the very first timestep after the dump command is invoked.  This will
always occur if the current timestep is a multiple of N, the frequency
specified in the :doc:`dump <dump>` command, including timestep 0.  But
if this is not the case, a dump snapshot will only be written if the
setting of this keyword is *yes*\ .  If it is *no*\ , which is the
default, then it will not be written.


----------


The *flush* keyword determines whether a flush operation is invoked
after a dump snapshot is written to the dump file.  A flush insures
the output in that file is current (no buffering by the OS), even if
LAMMPS halts before the simulation completes.  Flushes cannot be
performed with dump style *xtc*\ .


----------


The *format* keyword can be used to change the default numeric format
output by the text-based dump styles: *atom*\ , *custom*\ , *cfg*\ , and
*xyz* styles, and their MPIIO variants.  Only the *line* or *none*
options can be used with the *atom* and *xyz* styles.

All the specified format strings are C-style formats, e.g. as used by
the C/C++ printf() command.  The *line* keyword takes a single
argument which is the format string for an entire line of output for
each atom (do not include a trailing "\n"), with N fields, which you
must enclose in quotes if it is more than one field.  The *int* and
*float* keywords take a single format argument and are applied to all
integer or floating-point quantities output.  The setting for *M
string* also takes a single format argument which is used for the Mth
value output in each line, e.g. the 5th column is output in high
precision for "format 5 %20.15g".

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
   appropriate for an 8-byte signed integer, e.g. one with "%ld", if
   LAMMPS was compiled with the -DLAMMPS\_BIGBIG option for 8-byte IDs.

.. note::

   Any value written to a text-based dump file that is a per-atom
   quantity calculated by a :doc:`compute <compute>` or :doc:`fix <fix>` is
   stored internally as a floating-point value.  If the value is actually
   an integer and you wish it to appear in the text dump file as a
   (large) integer, then you need to use an appropriate format.  For
   example, these commands:


.. parsed-literal::

   compute     1 all property/local batom1 batom2
   dump        1 all local 100 tmp.bonds index c_1[1] c_1[2]
   dump_modify 1 format "%d %0.0f %0.0f"

will output the two atom IDs for atoms in each bond as integers.  If
the dump\_modify command were omitted, they would appear as
floating-point values, assuming they were large integers (more than 6
digits).  The "index" keyword should use the "%d" format since it is
not generated by a compute or fix, and is stored internally as an
integer.


----------


The *fileper* keyword is documented below with the *nfile* keyword.


----------


The *image* keyword applies only to the dump *atom* style.  If the
image value is *yes*\ , 3 flags are appended to each atom's coords which
are the absolute box image of the atom in each dimension.  For
example, an x image flag of -2 with a normalized coord of 0.5 means
the atom is in the center of the box, but has passed through the box
boundary 2 times and is really 2 box lengths to the left of its
current coordinate.  Note that for dump style *custom* these various
values can be printed in the dump file by using the appropriate atom
attributes in the dump command itself.


----------


The *label* keyword applies only to the dump *local* style.  When
it writes local information, such as bond or angle topology
to a dump file, it will use the specified *label* to format
the header.  By default this includes 2 lines:


.. parsed-literal::

   ITEM: NUMBER OF ENTRIES
   ITEM: ENTRIES ...

The word "ENTRIES" will be replaced with the string specified,
e.g. BONDS or ANGLES.


----------


The *maxfiles* keyword can only be used when a '\*' wildcard is
included in the dump file name, i.e. when writing a new file(s) for
each snapshot.  The specified *Fmax* is how many snapshots will be
kept.  Once this number is reached, the file(s) containing the oldest
snapshot is deleted before a new dump file is written.  If the
specified *Fmax* <= 0, then all files are retained.

This can be useful for debugging, especially if you don't know on what
timestep something bad will happen, e.g. when LAMMPS will exit with an
error.  You can dump every timestep, and limit the number of dump
files produced, even if you run for 1000s of steps.


----------


The *nfile* or *fileper* keywords can be used in conjunction with the
"%" wildcard character in the specified dump file name, for all dump
styles except the *dcd*\ , *image*\ , *movie*\ , *xtc*\ , and *xyz* styles
(for which "%" is not allowed).  As explained on the :doc:`dump <dump>`
command doc page, the "%" character causes the dump file to be written
in pieces, one piece for each of P processors.  By default P = the
number of processors the simulation is running on.  The *nfile* or
*fileper* keyword can be used to set P to a smaller value, which can
be more efficient when running on a large number of processors.

The *nfile* keyword sets P to the specified Nf value.  For example, if
Nf = 4, and the simulation is running on 100 processors, 4 files will
be written, by processors 0,25,50,75.  Each will collect information
from itself and the next 24 processors and write it to a dump file.

For the *fileper* keyword, the specified value of Np means write one
file for every Np processors.  For example, if Np = 4, every 4th
processor (0,4,8,12,etc) will collect information from itself and the
next 3 processors and write it to a dump file.


----------


The *pad* keyword only applies when the dump filename is specified
with a wildcard "\*" character which becomes the timestep.  If *pad* is
0, which is the default, the timestep is converted into a string of
unpadded length, e.g. 100 or 12000 or 2000000.  When *pad* is
specified with *Nchar* > 0, the string is padded with leading zeroes
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
*yes*\ , then all atoms will be remapped to the periodic box before the
snapshot is written, then restored to their original position.  If it
is set to *no* they will not be.  The *no* setting is the default
because it requires no extra computation.


----------


The *precision* keyword only applies to the dump *xtc* style.  A
specified value of N means that coordinates are stored to 1/N
nanometer accuracy, e.g. for N = 1000, the coordinates are written to
1/1000 nanometer accuracy.


----------


The *refresh* keyword only applies to the dump *custom*\ , *cfg*\ ,
*image*\ , and *movie* styles.  It allows an "incremental" dump file to
be written, by refreshing a compute that is used as a threshold for
determining which atoms are included in a dump snapshot.  The
specified *c\_ID* gives the ID of the compute.  It is prefixed by "c\_"
to indicate a compute, which is the only current option.  At some
point, other options may be added, e.g. fixes or variables.

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


.. parsed-literal::

   variable        Dhop equal 0.6
   variable        check atom "c_dsp[4] > v_Dhop"
   compute         dsp all displace/atom refresh check
   dump            1 all custom 20 tmp.dump id type x y z
   dump_modify     1 append yes thresh c_dsp[4] > ${Dhop} refresh c_dsp

The :doc:`compute displace/atom <compute_displace_atom>` command
calculates the displacement of each atom from its reference position.
The "4" index is the scalar displacement; 1,2,3 are the xyz components
of the displacement.  The :doc:`dump_modify thresh <dump_modify>`
command will cause only atoms that have displaced more than 0.6
Angstroms to be output on a given snapshot (assuming metal units).
However, note that when an atom is output, we also need to update the
reference position for that atom to its new coordinates.  So that it
will not be output in every snapshot thereafter.  That reference
position is stored by :doc:`compute displace/atom <compute_displace_atom>`.  So the dump\_modify
*refresh* option triggers a call to compute displace/atom at the end
of every dump to perform that update.  The *refresh check* option
shown as part of the :doc:`compute displace/atom <compute_displace_atom>` command enables the compute
to respond to the call from the dump command, and update the
appropriate reference positions.  This is done be defining an
:doc:`atom-style variable <variable>`, *check* in this example, which
calculates a Boolean value (0 or 1) for each atom, based on the same
criterion used by dump\_modify thresh.

See the :doc:`compute displace/atom <compute_displace_atom>` command for
more details, including an example of how to produce output that
includes an initial snapshot with the reference position of all atoms.

Note that only computes with a *refresh* option will work with
dump\_modify refresh.  See individual compute doc pages for details.
Currently, only compute displace/atom supports this option.  Others
may be added at some point.  If you use a compute that doesn't support
refresh operations, LAMMPS will not complain; dump\_modify refresh will
simply do nothing.


----------


The *region* keyword only applies to the dump *custom*\ , *cfg*\ ,
*image*\ , and *movie* styles.  If specified, only atoms in the region
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
(e.g. Angstroms or sigma).


----------


The *sfactor* and *tfactor* keywords only apply to the dump *xtc*
style.  They allow customization of the unit conversion factors used
when writing to XTC files.  By default they are initialized for
whatever :doc:`units <units>` style is being used, to write out
coordinates in nanometers and time in picoseconds.  I.e. for *real*
units, LAMMPS defines *sfactor* = 0.1 and *tfactor* = 0.001, since the
Angstroms and fmsec used by *real* units are 0.1 nm and 0.001 psec
respectively.  If you are using a units system with distance and time
units far from nm and psec, you may wish to write XTC files with
different units, since the compression algorithm used in XTC files is
most effective when the typical magnitude of position data is between
10.0 and 0.1.


----------


The *sort* keyword determines whether lines of per-atom output in a
snapshot are sorted or not.  A sort value of *off* means they will
typically be written in indeterminate order, either in serial or
parallel.  This is the case even in serial if the :doc:`atom_modify sort <atom_modify>` option is turned on, which it is by default, to
improve performance.  A sort value of *id* means sort the output by
atom ID.  A sort value of N or -N means sort the output by the value
in the Nth column of per-atom info in either ascending or descending
order.

The dump *local* style cannot be sorted by atom ID, since there are
typically multiple lines of output per atom.  Some dump styles, such
as *dcd* and *xtc*\ , require sorting by atom ID to format the output
file correctly.  If multiple processors are writing the dump file, via
the "%" wildcard in the dump filename, then sorting cannot be
performed.

.. note::

   Unless it is required by the dump style, sorting dump file
   output requires extra overhead in terms of CPU and communication cost,
   as well as memory, versus unsorted output.


----------


The *thermo* keyword only applies the dump *netcdf* style.  It
triggers writing of :doc:`thermo <thermo>` information to the dump file
alongside per-atom data.  The values included in the dump file are
identical to the values specified by :doc:`thermo_style <thermo_style>`.


----------


The *thresh* keyword only applies to the dump *custom*\ , *cfg*\ ,
*image*\ , and *movie* styles.  Multiple thresholds can be specified.
Specifying *none* turns off all threshold criteria.  If thresholds are
specified, only atoms whose attributes meet all the threshold criteria
are written to the dump file or included in the image.  The possible
attributes that can be tested for are the same as those that can be
specified in the :doc:`dump custom <dump>` command, with the exception
of the *element* attribute, since it is not a numeric value.  Note
that a different attributes can be used than those output by the :doc:`dump custom <dump>` command.  E.g. you can output the coordinates and
stress of atoms whose energy is above some threshold.

If an atom-style variable is used as the attribute, then it can
produce continuous numeric values or effective Boolean 0/1 values
which may be useful for the comparison operator.  Boolean values can
be generated by variable formulas that use comparison or Boolean math
operators or special functions like gmask() and rmask() and grmask().
See the :doc:`variable <variable>` command doc page for details.

The specified value must be a simple numeric value or the word LAST.
If LAST is used, it refers to the value of the attribute the last time
the dump command was invoked to produce a snapshot.  This is a way to
only dump atoms whose attribute has changed (or not changed).
Three examples follow.


.. parsed-literal::

   dump_modify ... thresh ix != LAST

This will dump atoms which have crossed the periodic x boundary of the
simulation box since the last dump.  (Note that atoms that crossed
once and then crossed back between the two dump timesteps would not be
included.)


.. parsed-literal::

   region foo sphere 10 20 10 15
   variable inregion atom rmask(foo)
   dump_modify ... thresh v_inregion \|\^ LAST

This will dump atoms which crossed the boundary of the spherical
region since the last dump.


.. parsed-literal::

   variable charge atom "(q > 0.5) \|\| (q < -0.5)"
   dump_modify ... thresh v_charge \|\^ LAST

This will dump atoms whose charge has changed from an absolute value
less than 1/2 to greater than 1/2 (or vice versa) since the last dump.
E.g. due to reactions and subsequent charge equilibration in a
reactive force field.

The choice of operators listed above are the usual comparison
operators.  The XOR operation (exclusive or) is also included as "\|\^".
In this context, XOR means that if either the attribute or value is
0.0 and the other is non-zero, then the result is "true" and the
threshold criterion is met.  Otherwise it is not met.


----------


The *time* keyword only applies to the dump *atom*\ , *custom*\ , and
*local* styles (and their COMPRESS package versions *atom/gz*\ ,
*custom/gz* and *local/gz*\ ). If set to *yes*\ , each frame will will
contain two extra lines before the "ITEM: TIMESTEP" entry:


.. parsed-literal::

   ITEM: TIME
   \<elapsed time\>

This will output the current elapsed simulation time in current
time units equivalent to the :doc:`thermo keyword <thermo_style>` *time*\ .
This is to simplify post-processing of trajectories using a variable time
step, e.g. when using :doc:`fix dt/reset <fix_dt_reset>`.
The default setting is *no*\ .


----------


The *units* keyword only applies to the dump *atom*\ , *custom*\ , and
*local* styles (and their COMPRESS package versions *atom/gz*\ ,
*custom/gz* and *local/gz*\ ). If set to *yes*\ , each individual dump
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
If set to *yes*\ , coordinates will be written "unwrapped" by the image
flags for each atom.  Unwrapped means that if the atom has passed through
a periodic boundary one or more times, the value is printed for what
the coordinate would be if it had not been wrapped back into the
periodic box.  Note that these coordinates may thus be far outside the
box size stored with the snapshot.


----------


These keywords apply only to the :doc:`dump image <dump_image>` and
:doc:`dump movie <dump_image>` styles.  Any keyword that affects an
image, also affects a movie, since the movie is simply a collection of
images.  Some of the keywords only affect the :doc:`dump movie <dump_image>` style.  The descriptions give details.


----------


The *acolor* keyword can be used with the :doc:`dump image <dump_image>`
command, when its atom color setting is *type*\ , to set the color that
atoms of each type will be drawn in the image.

The specified *type* should be an integer from 1 to Ntypes = the
number of atom types.  A wildcard asterisk can be used in place of or
in conjunction with the *type* argument to specify a range of atom
types.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the
number of atom types, then an asterisk with no numeric values means
all types from 1 to N.  A leading asterisk means all types from 1 to n
(inclusive).  A trailing asterisk means all types from n to N
(inclusive).  A middle asterisk means all types from m to n
(inclusive).

The specified *color* can be a single color which is any of the 140
pre-defined colors (see below) or a color name defined by the
dump\_modify color option.  Or it can be two or more colors separated
by a "/" character, e.g. red/green/blue.  In the former case, that
color is assigned to all the specified atom types.  In the latter
case, the list of colors are assigned in a round-robin fashion to each
of the specified atom types.


----------


The *adiam* keyword can be used with the :doc:`dump image <dump_image>`
command, when its atom diameter setting is *type*\ , to set the size
that atoms of each type will be drawn in the image.  The specified
*type* should be an integer from 1 to Ntypes.  As with the *acolor*
keyword, a wildcard asterisk can be used as part of the *type*
argument to specify a range of atom types.  The specified *diam* is
the size in whatever distance :doc:`units <units>` the input script is
using, e.g. Angstroms.


----------


The *amap* keyword can be used with the :doc:`dump image <dump_image>`
command, with its *atom* keyword, when its atom setting is an
atom-attribute, to setup a color map.  The color map is used to assign
a specific RGB (red/green/blue) color value to an individual atom when
it is drawn, based on the atom's attribute, which is a numeric value,
e.g. its x-component of velocity if the atom-attribute "vx" was
specified.

The basic idea of a color map is that the atom-attribute will be
within a range of values, and that range is associated with a series
of colors (e.g. red, blue, green).  An atom's specific value (vx =
-3.2) can then mapped to the series of colors (e.g. halfway between
red and blue), and a specific color is determined via an interpolation
procedure.

There are many possible options for the color map, enabled by the
*amap* keyword.  Here are the details.

The *lo* and *hi* settings determine the range of values allowed for
the atom attribute.  If numeric values are used for *lo* and/or *hi*\ ,
then values that are lower/higher than that value are set to the
value.  I.e. the range is static.  If *lo* is specified as *min* or
*hi* as *max* then the range is dynamic, and the lower and/or
upper bound will be calculated each time an image is drawn, based
on the set of atoms being visualized.

The *style* setting is two letters, such as "ca".  The first letter is
either "c" for continuous, "d" for discrete, or "s" for sequential.
The second letter is either "a" for absolute, or "f" for fractional.

A continuous color map is one in which the color changes continuously
from value to value within the range.  A discrete color map is one in
which discrete colors are assigned to sub-ranges of values within the
range.  A sequential color map is one in which discrete colors are
assigned to a sequence of sub-ranges of values covering the entire
range.

An absolute color map is one in which the values to which colors are
assigned are specified explicitly as values within the range.  A
fractional color map is one in which the values to which colors are
assigned are specified as a fractional portion of the range.  For
example if the range is from -10.0 to 10.0, and the color red is to be
assigned to atoms with a value of 5.0, then for an absolute color map
the number 5.0 would be used.  But for a fractional map, the number
0.75 would be used since 5.0 is 3/4 of the way from -10.0 to 10.0.

The *delta* setting must be specified for all styles, but is only used
for the sequential style; otherwise the value is ignored.  It
specifies the bin size to use within the range for assigning
consecutive colors to.  For example, if the range is from -10.0 to
10.0 and a *delta* of 1.0 is used, then 20 colors will be assigned to
the range.  The first will be from -10.0 <= color1 < -9.0, then 2nd
from -9.0 <= color2 < -8.0, etc.

The *N* setting is how many entries follow.  The format of the entries
depends on whether the color map style is continuous, discrete or
sequential.  In all cases the *color* setting can be any of the 140
pre-defined colors (see below) or a color name defined by the
dump\_modify color option.

For continuous color maps, each entry has a *value* and a *color*\ .
The *value* is either a number within the range of values or *min* or
*max*\ .  The *value* of the first entry must be *min* and the *value*
of the last entry must be *max*\ .  Any entries in between must have
increasing values.  Note that numeric values can be specified either
as absolute numbers or as fractions (0.0 to 1.0) of the range,
depending on the "a" or "f" in the style setting for the color map.

Here is how the entries are used to determine the color of an
individual atom, given the value X of its atom attribute.  X will fall
between 2 of the entry values.  The color of the atom is linearly
interpolated (in each of the RGB values) between the 2 colors
associated with those entries.  For example, if X = -5.0 and the 2
surrounding entries are "red" at -10.0 and "blue" at 0.0, then the
atom's color will be halfway between "red" and "blue", which happens
to be "purple".

For discrete color maps, each entry has a *lo* and *hi* value and a
*color*\ .  The *lo* and *hi* settings are either numbers within the
range of values or *lo* can be *min* or *hi* can be *max*\ .  The *lo*
and *hi* settings of the last entry must be *min* and *max*\ .  Other
entries can have any *lo* and *hi* values and the sub-ranges of
different values can overlap.  Note that numeric *lo* and *hi* values
can be specified either as absolute numbers or as fractions (0.0 to
1.0) of the range, depending on the "a" or "f" in the style setting
for the color map.

Here is how the entries are used to determine the color of an
individual atom, given the value X of its atom attribute.  The entries
are scanned from first to last.  The first time that *lo* <= X <=
*hi*\ , X is assigned the color associated with that entry.  You can
think of the last entry as assigning a default color (since it will
always be matched by X), and the earlier entries as colors that
override the default.  Also note that no interpolation of a color RGB
is done.  All atoms will be drawn with one of the colors in the list
of entries.

For sequential color maps, each entry has only a *color*\ .  Here is how
the entries are used to determine the color of an individual atom,
given the value X of its atom attribute.  The range is partitioned
into N bins of width *binsize*\ .  Thus X will fall in a specific bin
from 1 to N, say the Mth bin.  If it falls on a boundary between 2
bins, it is considered to be in the higher of the 2 bins.  Each bin is
assigned a color from the E entries.  If E < N, then the colors are
repeated.  For example if 2 entries with colors red and green are
specified, then the odd numbered bins will be red and the even bins
green.  The color of the atom is the color of its bin.  Note that the
sequential color map is really a shorthand way of defining a discrete
color map without having to specify where all the bin boundaries are.

Here is an example of using a sequential color map to color all the
atoms in individual molecules with a different color.  See the
examples/pour/in.pour.2d.molecule input script for an example of how
this is used.


.. parsed-literal::

   variable        colors string &
                   "red green blue yellow white &
                   purple pink orange lime gray"
   variable        mol atom mol%10
   dump            1 all image 250 image.\*.jpg v_mol type &
                   zoom 1.6 adiam 1.5
   dump_modify     1 pad 5 amap 0 10 sa 1 10 ${colors}

In this case, 10 colors are defined, and molecule IDs are
mapped to one of the colors, even if there are 1000s of molecules.


----------


The *backcolor* sets the background color of the images.  The color
name can be any of the 140 pre-defined colors (see below) or a color
name defined by the dump\_modify color option.


----------


The *bcolor* keyword can be used with the :doc:`dump image <dump_image>`
command, with its *bond* keyword, when its color setting is *type*\ , to
set the color that bonds of each type will be drawn in the image.

The specified *type* should be an integer from 1 to Nbondtypes = the
number of bond types.  A wildcard asterisk can be used in place of or
in conjunction with the *type* argument to specify a range of bond
types.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the
number of bond types, then an asterisk with no numeric values means
all types from 1 to N.  A leading asterisk means all types from 1 to n
(inclusive).  A trailing asterisk means all types from n to N
(inclusive).  A middle asterisk means all types from m to n
(inclusive).

The specified *color* can be a single color which is any of the 140
pre-defined colors (see below) or a color name defined by the
dump\_modify color option.  Or it can be two or more colors separated
by a "/" character, e.g. red/green/blue.  In the former case, that
color is assigned to all the specified bond types.  In the latter
case, the list of colors are assigned in a round-robin fashion to each
of the specified bond types.


----------


The *bdiam* keyword can be used with the :doc:`dump image <dump_image>`
command, with its *bond* keyword, when its diam setting is *type*\ , to
set the diameter that bonds of each type will be drawn in the image.
The specified *type* should be an integer from 1 to Nbondtypes.  As
with the *bcolor* keyword, a wildcard asterisk can be used as part of
the *type* argument to specify a range of bond types.  The specified
*diam* is the size in whatever distance :doc:`units <units>` you are
using, e.g. Angstroms.


----------


The *bitrate* keyword can be used with the :doc:`dump movie <dump_image>` command to define the size of the resulting
movie file and its quality via setting how many kbits per second are
to be used for the movie file. Higher bitrates require less
compression and will result in higher quality movies.  The quality is
also determined by the compression format and encoder.  The default
setting is 2000 kbit/s, which will result in average quality with
older compression formats.

.. note::

   Not all movie file formats supported by dump movie allow the
   bitrate to be set.  If not, the setting is silently ignored.


----------


The *boxcolor* keyword sets the color of the simulation box drawn
around the atoms in each image as well as the color of processor
sub-domain boundaries.  See the "dump image box" command for how to
specify that a box be drawn via the *box* keyword, and the sub-domain
boundaries via the *subbox* keyword.  The color name can be any of the
140 pre-defined colors (see below) or a color name defined by the
dump\_modify color option.


----------


The *color* keyword allows definition of a new color name, in addition
to the 140-predefined colors (see below), and associates 3
red/green/blue RGB values with that color name.  The color name can
then be used with any other dump\_modify keyword that takes a color
name as a value.  The RGB values should each be floating point values
between 0.0 and 1.0 inclusive.

When a color name is converted to RGB values, the user-defined color
names are searched first, then the 140 pre-defined color names.  This
means you can also use the *color* keyword to overwrite one of the
pre-defined color names with new RBG values.


----------


The *framerate* keyword can be used with the :doc:`dump movie <dump_image>` command to define the duration of the resulting
movie file.  Movie files written by the dump *movie* command have a
default frame rate of 24 frames per second and the images generated
will be converted at that rate.  Thus a sequence of 1000 dump images
will result in a movie of about 42 seconds.  To make a movie run
longer you can either generate images more frequently or lower the
frame rate.  To speed a movie up, you can do the inverse.  Using a
frame rate higher than 24 is not recommended, as it will result in
simply dropping the rendered images. It is more efficient to dump
images less frequently.


----------


Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump image <dump_image>`, :doc:`undump <undump>`

Default
"""""""

The option defaults are

* append = no
* buffer = yes for dump styles *atom*\ , *custom*\ , *loca*\ , and *xyz*
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
* sort = off for dump styles *atom*\ , *custom*\ , *cfg*\ , and *local*
* sort = id for dump styles *dcd*\ , *xtc*\ , and *xyz*
* thresh = none
* units = no
* unwrap = no

* acolor = \* red/green/blue/yellow/aqua/cyan
* adiam = \* 1.0
* amap = min max cf 0.0 2 min blue max red
* backcolor = black
* bcolor = \* red/green/blue/yellow/aqua/cyan
* bdiam = \* 0.5
* bitrate = 2000
* boxcolor = yellow
* color = 140 color names are pre-defined as listed below
* framerate = 24


----------


These are the standard 109 element names that LAMMPS pre-defines for
use with the :doc:`dump image <dump_image>` and dump\_modify commands.

* 1-10 = "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"
* 11-20 = "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca"
* 21-30 = "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"
* 31-40 = "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr"
* 41-50 = "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn"
* 51-60 = "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd"
* 61-70 = "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"
* 71-80 = "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"
* 81-90 = "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th"
* 91-100 = "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm"
* 101-109 = "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"


----------


These are the 140 colors that LAMMPS pre-defines for use with the
:doc:`dump image <dump_image>` and dump\_modify commands.  Additional
colors can be defined with the dump\_modify color command.  The 3
numbers listed for each name are the RGB (red/green/blue) values.
Divide each value by 255 to get the equivalent 0.0 to 1.0 value.

+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| aliceblue = 240, 248, 255     | antiquewhite = 250, 235, 215         | aqua = 0, 255, 255              | aquamarine = 127, 255, 212     | azure = 240, 255, 255          |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| beige = 245, 245, 220         | bisque = 255, 228, 196               | black = 0, 0, 0                 | blanchedalmond = 255, 255, 205 | blue = 0, 0, 255               |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| blueviolet = 138, 43, 226     | brown = 165, 42, 42                  | burlywood = 222, 184, 135       | cadetblue = 95, 158, 160       | chartreuse = 127, 255, 0       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| chocolate = 210, 105, 30      | coral = 255, 127, 80                 | cornflowerblue = 100, 149, 237  | cornsilk = 255, 248, 220       | crimson = 220, 20, 60          |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| cyan = 0, 255, 255            | darkblue = 0, 0, 139                 | darkcyan = 0, 139, 139          | darkgoldenrod = 184, 134, 11   | darkgray = 169, 169, 169       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| darkgreen = 0, 100, 0         | darkkhaki = 189, 183, 107            | darkmagenta = 139, 0, 139       | darkolivegreen = 85, 107, 47   | darkorange = 255, 140, 0       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| darkorchid = 153, 50, 204     | darkred = 139, 0, 0                  | darksalmon = 233, 150, 122      | darkseagreen = 143, 188, 143   | darkslateblue = 72, 61, 139    |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| darkslategray = 47, 79, 79    | darkturquoise = 0, 206, 209          | darkviolet = 148, 0, 211        | deeppink = 255, 20, 147        | deepskyblue = 0, 191, 255      |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| dimgray = 105, 105, 105       | dodgerblue = 30, 144, 255            | firebrick = 178, 34, 34         | floralwhite = 255, 250, 240    | forestgreen = 34, 139, 34      |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| fuchsia = 255, 0, 255         | gainsboro = 220, 220, 220            | ghostwhite = 248, 248, 255      | gold = 255, 215, 0             | goldenrod = 218, 165, 32       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| gray = 128, 128, 128          | green = 0, 128, 0                    | greenyellow = 173, 255, 47      | honeydew = 240, 255, 240       | hotpink = 255, 105, 180        |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| indianred = 205, 92, 92       | indigo = 75, 0, 130                  | ivory = 255, 240, 240           | khaki = 240, 230, 140          | lavender = 230, 230, 250       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| lavenderblush = 255, 240, 245 | lawngreen = 124, 252, 0              | lemonchiffon = 255, 250, 205    | lightblue = 173, 216, 230      | lightcoral = 240, 128, 128     |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| lightcyan = 224, 255, 255     | lightgoldenrodyellow = 250, 250, 210 | lightgreen = 144, 238, 144      | lightgrey = 211, 211, 211      | lightpink = 255, 182, 193      |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| lightsalmon = 255, 160, 122   | lightseagreen = 32, 178, 170         | lightskyblue = 135, 206, 250    | lightslategray = 119, 136, 153 | lightsteelblue = 176, 196, 222 |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| lightyellow = 255, 255, 224   | lime = 0, 255, 0                     | limegreen = 50, 205, 50         | linen = 250, 240, 230          | magenta = 255, 0, 255          |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| maroon = 128, 0, 0            | mediumaquamarine = 102, 205, 170     | mediumblue = 0, 0, 205          | mediumorchid = 186, 85, 211    | mediumpurple = 147, 112, 219   |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| mediumseagreen = 60, 179, 113 | mediumslateblue = 123, 104, 238      | mediumspringgreen = 0, 250, 154 | mediumturquoise = 72, 209, 204 | mediumvioletred = 199, 21, 133 |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| midnightblue = 25, 25, 112    | mintcream = 245, 255, 250            | mistyrose = 255, 228, 225       | moccasin = 255, 228, 181       | navajowhite = 255, 222, 173    |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| navy = 0, 0, 128              | oldlace = 253, 245, 230              | olive = 128, 128, 0             | olivedrab = 107, 142, 35       | orange = 255, 165, 0           |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| orangered = 255, 69, 0        | orchid = 218, 112, 214               | palegoldenrod = 238, 232, 170   | palegreen = 152, 251, 152      | paleturquoise = 175, 238, 238  |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| palevioletred = 219, 112, 147 | papayawhip = 255, 239, 213           | peachpuff = 255, 239, 213       | peru = 205, 133, 63            | pink = 255, 192, 203           |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| plum = 221, 160, 221          | powderblue = 176, 224, 230           | purple = 128, 0, 128            | red = 255, 0, 0                | rosybrown = 188, 143, 143      |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| royalblue = 65, 105, 225      | saddlebrown = 139, 69, 19            | salmon = 250, 128, 114          | sandybrown = 244, 164, 96      | seagreen = 46, 139, 87         |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| seashell = 255, 245, 238      | sienna = 160, 82, 45                 | silver = 192, 192, 192          | skyblue = 135, 206, 235        | slateblue = 106, 90, 205       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| slategray = 112, 128, 144     | snow = 255, 250, 250                 | springgreen = 0, 255, 127       | steelblue = 70, 130, 180       | tan = 210, 180, 140            |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| teal = 0, 128, 128            | thistle = 216, 191, 216              | tomato = 253, 99, 71            | turquoise = 64, 224, 208       | violet = 238, 130, 238         |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| wheat = 245, 222, 179         | white = 255, 255, 255                | whitesmoke = 245, 245, 245      | yellow = 255, 255, 0           | yellowgreen = 154, 205, 50     |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
