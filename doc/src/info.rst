.. index:: info

info command
============

Syntax
""""""


.. parsed-literal::

   info args

* args = one or more of the following keywords: *out*\ , *all*\ , *system*\ , *memory*\ , *communication*\ , *computes*\ , *dumps*\ , *fixes*\ , *groups*\ , *regions*\ , *variables*\ , *coeffs*\ , *styles*\ , *time*\ , or *configuration*
* *out* values = *screen*\ , *log*\ , *append* filename, *overwrite* filename
* *styles* values = *all*\ , *angle*\ , *atom*\ , *bond*\ , *compute*\ , *command*\ , *dump*\ , *dihedral*\ , *fix*\ , *improper*\ , *integrate*\ , *kspace*\ , *minimize*\ , *pair*\ , *region*

Examples
""""""""


.. parsed-literal::

   info system
   info groups computes variables
   info all out log
   info all out append info.txt
   info styles all
   info styles atom

Description
"""""""""""

Print out information about the current internal state of the running
LAMMPS process. This can be helpful when debugging or validating
complex input scripts.  Several output categories are available and
one or more output category may be requested.

The *out* flag controls where the output is sent. It can only be sent
to one target. By default this is the screen, if it is active. The
*log* argument selects the log file instead. With the *append* and
*overwrite* option, followed by a filename, the output is written
to that file, which is either appended to or overwritten, respectively.

The *all* flag activates printing all categories listed below.

The *configuration* category prints some information about the
LAMMPS version as well as architecture and OS it is run on.

The *memory* category prints some information about the current
memory allocation of MPI rank 0 (this the amount of dynamically
allocated memory reported by LAMMPS classes). Where supported,
also some OS specific information about the size of the reserved
memory pool size (this is where malloc() and the new operator
request memory from) and the maximum resident set size is reported
(this is the maximum amount of physical memory occupied so far).

The *system* category prints a general system overview listing.  This
includes the unit style, atom style, number of atoms, bonds, angles,
dihedrals, and impropers and the number of the respective types, box
dimensions and properties, force computing styles and more.

The *communication* category prints a variety of information about
communication and parallelization: the MPI library version level,
the number of MPI ranks and OpenMP threads, the communication style
and layout, the processor grid dimensions, ghost atom communication
mode, cutoff, and related settings.

The *computes* category prints a list of all currently defined
computes, their IDs and styles and groups they operate on.

The *dumps* category prints a list of all currently active dumps,
their IDs, styles, filenames, groups, and dump frequencies.

The *fixes* category prints a list of all currently defined fixes,
their IDs and styles and groups they operate on.

The *groups* category prints a list of all currently defined groups.

The *regions* category prints a list of all currently defined regions,
their IDs and styles and whether "inside" or "outside" atoms are
selected.

The *variables* category prints a list of all currently defined
variables, their names, styles, definition and last computed value, if
available.

The *coeffs* category prints a list for each defined force style
(pair, bond, angle, dihedral, improper) indicating which of the
corresponding coefficients have been set. This can be very helpful
to debug error messages like "All pair coeffs are not set".

The *styles* category prints the list of styles available in the
current LAMMPS binary. It supports one of the following options
to control which category of styles is printed out:

* all
* angle
* atom
* bond
* compute
* command
* dump
* dihedral
* fix
* improper
* integrate
* kspace
* minimize
* pair
* region

The *time* category prints the accumulated CPU and wall time for the
process that writes output (usually MPI rank 0).

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`print <print>`

Default
"""""""

The *out* option has the default *screen*\ .

The *styles* option has the default *all*\ .
