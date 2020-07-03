Command-line options
====================

At run time, LAMMPS recognizes several optional command-line switches
which may be used in any order.  Either the full word or a one-or-two
letter abbreviation can be used:

* :ref:`-e or -echo <echo>`
* :ref:`-h or -help <help>`
* :ref:`-i or -in <file>`
* :ref:`-k or -kokkos <run-kokkos>`
* :ref:`-l or -log <log>`
* :ref:`-m or -mpicolor <mpicolor>`
* :ref:`-nc or -nocite <nocite>`
* :ref:`-pk or -package <package>`
* :ref:`-p or -partition <partition>`
* :ref:`-pl or -plog <plog>`
* :ref:`-ps or -pscreen <pscreen>`
* :ref:`-ro or -reorder <reorder>`
* :ref:`-r2data or -restart2data <restart2data>`
* :ref:`-r2dump or -restart2dump <restart2dump>`
* :ref:`-sc or -screen <screen>`
* :ref:`-sf or -suffix <suffix>`
* :ref:`-v or -var <var>`

For example, the lmp_mpi executable might be launched as follows:

.. code-block:: bash

   $ mpirun -np 16 lmp_mpi -v f tmp.out -l my.log -sc none -i in.alloy
   $ mpirun -np 16 lmp_mpi -var f tmp.out -log my.log -screen none -in in.alloy

----------

.. _echo:

**-echo style**

Set the style of command echoing.  The style can be *none* or *screen*
or *log* or *both*\ .  Depending on the style, each command read from
the input script will be echoed to the screen and/or logfile.  This
can be useful to figure out which line of your script is causing an
input error.  The default value is *log*\ .  The echo style can also be
set by using the :doc:`echo <echo>` command in the input script itself.

----------

.. _help:

**-help**

Print a brief help summary and a list of options compiled into this
executable for each LAMMPS style (atom_style, fix, compute,
pair_style, bond_style, etc).  This can tell you if the command you
want to use was included via the appropriate package at compile time.
LAMMPS will print the info and immediately exit if this switch is
used.

----------

.. _file:

**-in file**

Specify a file to use as an input script.  This is an optional switch
when running LAMMPS in one-partition mode.  If it is not specified,
LAMMPS reads its script from standard input, typically from a script
via I/O redirection; e.g. lmp_linux < in.run.  I/O redirection should
also work in parallel, but if it does not (in the unlikely case that
an MPI implementation does not support it), then use the -in flag.
Note that this is a required switch when running LAMMPS in
multi-partition mode, since multiple processors cannot all read from
stdin.

----------

.. _run-kokkos:

**-kokkos on/off keyword/value ...**

Explicitly enable or disable KOKKOS support, as provided by the KOKKOS
package.  Even if LAMMPS is built with this package, as described
in :doc:`Speed kokkos <Speed_kokkos>`, this switch must be set to enable
running with KOKKOS-enabled styles the package provides.  If the
switch is not set (the default), LAMMPS will operate as if the KOKKOS
package were not installed; i.e. you can run standard LAMMPS or with
the GPU or USER-OMP packages, for testing or benchmarking purposes.

Additional optional keyword/value pairs can be specified which
determine how Kokkos will use the underlying hardware on your
platform.  These settings apply to each MPI task you launch via the
"mpirun" or "mpiexec" command.  You may choose to run one or more MPI
tasks per physical node.  Note that if you are running on a desktop
machine, you typically have one physical node.  On a cluster or
supercomputer there may be dozens or 1000s of physical nodes.

Either the full word or an abbreviation can be used for the keywords.
Note that the keywords do not use a leading minus sign.  I.e. the
keyword is "t", not "-t".  Also note that each of the keywords has a
default setting.  Examples of when to use these options and what
settings to use on different platforms is given on the :doc:`Speed kokkos <Speed_kokkos>` doc page.

* d or device
* g or gpus
* t or threads
* n or numa

.. parsed-literal::

   device Nd

This option is only relevant if you built LAMMPS with CUDA=yes, you
have more than one GPU per node, and if you are running with only one
MPI task per node.  The Nd setting is the ID of the GPU on the node to
run on.  By default Nd = 0.  If you have multiple GPUs per node, they
have consecutive IDs numbered as 0,1,2,etc.  This setting allows you
to launch multiple independent jobs on the node, each with a single
MPI task per node, and assign each job to run on a different GPU.

.. parsed-literal::

   gpus Ng Ns

This option is only relevant if you built LAMMPS with CUDA=yes, you
have more than one GPU per node, and you are running with multiple MPI
tasks per node (up to one per GPU).  The Ng setting is how many GPUs
you will use.  The Ns setting is optional.  If set, it is the ID of a
GPU to skip when assigning MPI tasks to GPUs.  This may be useful if
your desktop system reserves one GPU to drive the screen and the rest
are intended for computational work like running LAMMPS.  By default
Ng = 1 and Ns is not set.

Depending on which flavor of MPI you are running, LAMMPS will look for
one of these 4 environment variables

.. parsed-literal::

   SLURM_LOCALID (various MPI variants compiled with SLURM support)
   MPT_LRANK (HPE MPI)
   MV2_COMM_WORLD_LOCAL_RANK (Mvapich)
   OMPI_COMM_WORLD_LOCAL_RANK (OpenMPI)

which are initialized by the "srun", "mpirun" or "mpiexec" commands.
The environment variable setting for each MPI rank is used to assign a
unique GPU ID to the MPI task.

.. parsed-literal::

   threads Nt

This option assigns Nt number of threads to each MPI task for
performing work when Kokkos is executing in OpenMP or pthreads mode.
The default is Nt = 1, which essentially runs in MPI-only mode.  If
there are Np MPI tasks per physical node, you generally want Np\*Nt =
the number of physical cores per node, to use your available hardware
optimally.  This also sets the number of threads used by the host when
LAMMPS is compiled with CUDA=yes.

.. parsed-literal::

   numa Nm

This option is only relevant when using pthreads with hwloc support.
In this case Nm defines the number of NUMA regions (typically sockets)
on a node which will be utilized by a single MPI rank.  By default Nm
= 1.  If this option is used the total number of worker-threads per
MPI rank is threads\*numa.  Currently it is always almost better to
assign at least one MPI rank per NUMA region, and leave numa set to
its default value of 1. This is because letting a single process span
multiple NUMA regions induces a significant amount of cross NUMA data
traffic which is slow.

----------

.. _log:

**-log file**

Specify a log file for LAMMPS to write status information to.  In
one-partition mode, if the switch is not used, LAMMPS writes to the
file log.lammps.  If this switch is used, LAMMPS writes to the
specified file.  In multi-partition mode, if the switch is not used, a
log.lammps file is created with high-level status information.  Each
partition also writes to a log.lammps.N file where N is the partition
ID.  If the switch is specified in multi-partition mode, the high-level
logfile is named "file" and each partition also logs information to a
file.N.  For both one-partition and multi-partition mode, if the
specified file is "none", then no log files are created.  Using a
:doc:`log <log>` command in the input script will override this setting.
Option -plog will override the name of the partition log files file.N.

----------

.. _mpicolor:

**-mpicolor** color

If used, this must be the first command-line argument after the LAMMPS
executable name.  It is only used when LAMMPS is launched by an mpirun
command which also launches another executable(s) at the same time.
(The other executable could be LAMMPS as well.)  The color is an
integer value which should be different for each executable (another
application may set this value in a different way).  LAMMPS and the
other executable(s) perform an MPI_Comm_split() with their own colors
to shrink the MPI_COMM_WORLD communication to be the subset of
processors they are actually running on.

Currently, this is only used in LAMMPS to perform client/server
messaging with another application.  LAMMPS can act as either a client
or server (or both).  More details are given on the :doc:`Howto client/server <Howto_client_server>` doc page.

Specifically, this refers to the "mpi/one" mode of messaging provided
by the :doc:`message <message>` command and the CSlib library LAMMPS
links with from the lib/message directory.  See the
:doc:`message <message>` command for more details.

----------

.. _nocite:

**-nocite**

Disable writing the log.cite file which is normally written to list
references for specific cite-able features used during a LAMMPS run.
See the `citation page <https://lammps.sandia.gov/cite.html>`_ for more
details.

----------

.. _package:

**-package style args ....**

Invoke the :doc:`package <package>` command with style and args.  The
syntax is the same as if the command appeared at the top of the input
script.  For example "-package gpu 2" or "-pk gpu 2" is the same as
:doc:`package gpu 2 <package>` in the input script.  The possible styles
and args are documented on the :doc:`package <package>` doc page.  This
switch can be used multiple times, e.g. to set options for the
USER-INTEL and USER-OMP packages which can be used together.

Along with the "-suffix" command-line switch, this is a convenient
mechanism for invoking accelerator packages and their options without
having to edit an input script.

----------

.. _partition:

**-partition 8x2 4 5 ...**

Invoke LAMMPS in multi-partition mode.  When LAMMPS is run on P
processors and this switch is not used, LAMMPS runs in one partition,
i.e. all P processors run a single simulation.  If this switch is
used, the P processors are split into separate partitions and each
partition runs its own simulation.  The arguments to the switch
specify the number of processors in each partition.  Arguments of the
form MxN mean M partitions, each with N processors.  Arguments of the
form N mean a single partition with N processors.  The sum of
processors in all partitions must equal P.  Thus the command
"-partition 8x2 4 5" has 10 partitions and runs on a total of 25
processors.

Running with multiple partitions can be useful for running
:doc:`multi-replica simulations <Howto_replica>`, where each replica
runs on one or a few processors.  Note that with MPI installed on a
machine (e.g. your desktop), you can run on more (virtual) processors
than you have physical processors.

To run multiple independent simulations from one input script, using
multiple partitions, see the :doc:`Howto multiple <Howto_multiple>` doc
page.  World- and universe-style :doc:`variables <variable>` are useful
in this context.

----------

.. _plog:

**-plog file**

Specify the base name for the partition log files, so partition N
writes log information to file.N. If file is none, then no partition
log files are created.  This overrides the filename specified in the
-log command-line option.  This option is useful when working with
large numbers of partitions, allowing the partition log files to be
suppressed (-plog none) or placed in a sub-directory (-plog
replica_files/log.lammps) If this option is not used the log file for
partition N is log.lammps.N or whatever is specified by the -log
command-line option.

----------

.. _pscreen:

**-pscreen file**

Specify the base name for the partition screen file, so partition N
writes screen information to file.N. If file is none, then no
partition screen files are created.  This overrides the filename
specified in the -screen command-line option.  This option is useful
when working with large numbers of partitions, allowing the partition
screen files to be suppressed (-pscreen none) or placed in a
sub-directory (-pscreen replica_files/screen).  If this option is not
used the screen file for partition N is screen.N or whatever is
specified by the -screen command-line option.

----------

.. _reorder:

**-reorder**

This option has 2 forms:

.. parsed-literal::

   -reorder nth N
   -reorder custom filename

Reorder the processors in the MPI communicator used to instantiate
LAMMPS, in one of several ways.  The original MPI communicator ranks
all P processors from 0 to P-1.  The mapping of these ranks to
physical processors is done by MPI before LAMMPS begins.  It may be
useful in some cases to alter the rank order.  E.g. to insure that
cores within each node are ranked in a desired order.  Or when using
the :doc:`run_style verlet/split <run_style>` command with 2 partitions
to insure that a specific Kspace processor (in the second partition) is
matched up with a specific set of processors in the first partition.
See the :doc:`Speed tips <Speed_tips>` doc page for more details.

If the keyword *nth* is used with a setting *N*\ , then it means every
Nth processor will be moved to the end of the ranking.  This is useful
when using the :doc:`run_style verlet/split <run_style>` command with 2
partitions via the -partition command-line switch.  The first set of
processors will be in the first partition, the second set in the second
partition.  The -reorder command-line switch can alter this so that
the first N procs in the first partition and one proc in the second partition
will be ordered consecutively, e.g. as the cores on one physical node.
This can boost performance.  For example, if you use "-reorder nth 4"
and "-partition 9 3" and you are running on 12 processors, the
processors will be reordered from

.. parsed-literal::

   0 1 2 3 4 5 6 7 8 9 10 11

to

.. parsed-literal::

   0 1 2 4 5 6 8 9 10 3 7 11

so that the processors in each partition will be

.. parsed-literal::

   0 1 2 4 5 6 8 9 10
   3 7 11

See the "processors" command for how to insure processors from each
partition could then be grouped optimally for quad-core nodes.

If the keyword is *custom*\ , then a file that specifies a permutation
of the processor ranks is also specified.  The format of the reorder
file is as follows.  Any number of initial blank or comment lines
(starting with a "#" character) can be present.  These should be
followed by P lines of the form:

.. parsed-literal::

   I J

where P is the number of processors LAMMPS was launched with.  Note
that if running in multi-partition mode (see the -partition switch
above) P is the total number of processors in all partitions.  The I
and J values describe a permutation of the P processors.  Every I and
J should be values from 0 to P-1 inclusive.  In the set of P I values,
every proc ID should appear exactly once.  Ditto for the set of P J
values.  A single I,J pairing means that the physical processor with
rank I in the original MPI communicator will have rank J in the
reordered communicator.

Note that rank ordering can also be specified by many MPI
implementations, either by environment variables that specify how to
order physical processors, or by config files that specify what
physical processors to assign to each MPI rank.  The -reorder switch
simply gives you a portable way to do this without relying on MPI
itself.  See the :doc:`processors out <processors>` command for how
to output info on the final assignment of physical processors to
the LAMMPS simulation domain.

----------

.. _restart2data:

**-restart2data restartfile [remap] datafile keyword value ...**

Convert the restart file into a data file and immediately exit.  This
is the same operation as if the following 2-line input script were
run:

.. code-block:: LAMMPS

   read_restart restartfile [remap]
   write_data datafile keyword value ...

The specified restartfile and/or datafile name may contain the wild-card
character "\*".  The restartfile name may also contain the wild-card
character "%".  The meaning of these characters is explained on the
:doc:`read_restart <read_restart>` and :doc:`write_data <write_data>` doc
pages.  The use of "%" means that a parallel restart file can be read.
Note that a filename such as file.\* may need to be enclosed in quotes or
the "\*" character prefixed with a backslash ("\") to avoid shell
expansion of the "\*" character.

Following restartfile argument, the optional word "remap" may be used.
This has the same effect like adding it to a
:doc:`read_restart <read_restart>` command, and operates as explained on
its doc page.  This is useful if reading the restart file triggers an
error that atoms have been lost.  In that case, use of the remap flag
should allow the data file to still be produced.

The syntax following restartfile (or remap), namely

.. parsed-literal::

   datafile keyword value ...

is identical to the arguments of the :doc:`write_data <write_data>`
command.  See its doc page for details.  This includes its
optional keyword/value settings.

----------

.. _restart2dump:

**-restart2dump restartfile [remap] group-ID dumpstyle dumpfile arg1 arg2 ...**

Convert the restart file into a dump file and immediately exit.  This
is the same operation as if the following 2-line input script were
run:

.. code-block:: LAMMPS

   read_restart restartfile [remap]
   write_dump group-ID dumpstyle dumpfile arg1 arg2 ...

Note that the specified restartfile and dumpfile names may contain
wild-card characters ("\*","%") as explained on the
:doc:`read_restart <read_restart>` and :doc:`write_dump <write_dump>` doc
pages.  The use of "%" means that a parallel restart file and/or
parallel dump file can be read and/or written.  Note that a filename
such as file.\* may need to be enclosed in quotes or the "\*" character
prefixed with a backslash ("\") to avoid shell expansion of the "\*"
character.

Note that following the restartfile argument, the optional word "remap"
can be used.  This has the effect as adding it to the
:doc:`read_restart <read_restart>` command, as explained on its doc page.
This is useful if reading the restart file triggers an error that atoms
have been lost.  In that case, use of the remap flag should allow the
dump file to still be produced.

The syntax following restartfile (or remap), namely

.. code-block:: LAMMPS

   group-ID dumpstyle dumpfile arg1 arg2 ...

is identical to the arguments of the :doc:`write_dump <write_dump>`
command.  See its doc page for details.  This includes what per-atom
fields are written to the dump file and optional dump_modify settings,
including ones that affect how parallel dump files are written, e.g.
the *nfile* and *fileper* keywords.  See the
:doc:`dump_modify <dump_modify>` doc page for details.

----------

.. _screen:

**-screen file**

Specify a file for LAMMPS to write its screen information to.  In
one-partition mode, if the switch is not used, LAMMPS writes to the
screen.  If this switch is used, LAMMPS writes to the specified file
instead and you will see no screen output.  In multi-partition mode,
if the switch is not used, high-level status information is written to
the screen.  Each partition also writes to a screen.N file where N is
the partition ID.  If the switch is specified in multi-partition mode,
the high-level screen dump is named "file" and each partition also
writes screen information to a file.N.  For both one-partition and
multi-partition mode, if the specified file is "none", then no screen
output is performed. Option -pscreen will override the name of the
partition screen files file.N.

----------

.. _suffix:

**-suffix style args**

Use variants of various styles if they exist.  The specified style can
be *gpu*\ , *intel*\ , *kk*\ , *omp*\ , *opt*\ , or *hybrid*\ .  These
refer to optional packages that LAMMPS can be built with, as described
in :doc:`Accelerate performance <Speed>`.  The "gpu" style corresponds to the
GPU package, the "intel" style to the USER-INTEL package, the "kk"
style to the KOKKOS package, the "opt" style to the OPT package, and
the "omp" style to the USER-OMP package. The hybrid style is the only
style that accepts arguments. It allows for two packages to be
specified. The first package specified is the default and will be used
if it is available. If no style is available for the first package,
the style for the second package will be used if available. For
example, "-suffix hybrid intel omp" will use styles from the
USER-INTEL package if they are installed and available, but styles for
the USER-OMP package otherwise.

Along with the "-package" command-line switch, this is a convenient
mechanism for invoking accelerator packages and their options without
having to edit an input script.

As an example, all of the packages provide a :doc:`pair_style lj/cut <pair_lj>` variant, with style names lj/cut/gpu,
lj/cut/intel, lj/cut/kk, lj/cut/omp, and lj/cut/opt.  A variant style
can be specified explicitly in your input script, e.g. pair_style
lj/cut/gpu.  If the -suffix switch is used the specified suffix
(gpu,intel,kk,omp,opt) is automatically appended whenever your input
script command creates a new :doc:`atom <atom_style>`,
:doc:`pair <pair_style>`, :doc:`fix <fix>`, :doc:`compute <compute>`, or
:doc:`run <run_style>` style.  If the variant version does not exist,
the standard version is created.

For the GPU package, using this command-line switch also invokes the
default GPU settings, as if the command "package gpu 1" were used at
the top of your input script.  These settings can be changed by using
the "-package gpu" command-line switch or the :doc:`package gpu <package>` command in your script.

For the USER-INTEL package, using this command-line switch also
invokes the default USER-INTEL settings, as if the command "package
intel 1" were used at the top of your input script.  These settings
can be changed by using the "-package intel" command-line switch or
the :doc:`package intel <package>` command in your script. If the
USER-OMP package is also installed, the hybrid style with "intel omp"
arguments can be used to make the omp suffix a second choice, if a
requested style is not available in the USER-INTEL package.  It will
also invoke the default USER-OMP settings, as if the command "package
omp 0" were used at the top of your input script.  These settings can
be changed by using the "-package omp" command-line switch or the
:doc:`package omp <package>` command in your script.

For the KOKKOS package, using this command-line switch also invokes
the default KOKKOS settings, as if the command "package kokkos" were
used at the top of your input script.  These settings can be changed
by using the "-package kokkos" command-line switch or the :doc:`package kokkos <package>` command in your script.

For the OMP package, using this command-line switch also invokes the
default OMP settings, as if the command "package omp 0" were used at
the top of your input script.  These settings can be changed by using
the "-package omp" command-line switch or the :doc:`package omp <package>` command in your script.

The :doc:`suffix <suffix>` command can also be used within an input
script to set a suffix, or to turn off or back on any suffix setting
made via the command line.

----------

.. _var:

**-var name value1 value2 ...**

Specify a variable that will be defined for substitution purposes when
the input script is read.  This switch can be used multiple times to
define multiple variables.  "Name" is the variable name which can be a
single character (referenced as $x in the input script) or a full
string (referenced as ${abc}).  An :doc:`index-style variable <variable>` will be created and populated with the
subsequent values, e.g. a set of filenames.  Using this command-line
option is equivalent to putting the line "variable name index value1
value2 ..."  at the beginning of the input script.  Defining an index
variable as a command-line argument overrides any setting for the same
index variable in the input script, since index variables cannot be
re-defined.

See the :doc:`variable <variable>` command for more info on defining
index and other kinds of variables and the :doc:`Commands parse <Commands_parse>` page for more info on using variables in
input scripts.

.. note::

   Currently, the command-line parser looks for arguments that
   start with "-" to indicate new switches.  Thus you cannot specify
   multiple variable values if any of them start with a "-", e.g. a
   negative numeric value.  It is OK if the first value1 starts with a
   "-", since it is automatically skipped.
