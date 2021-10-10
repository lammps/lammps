.. index:: processors

processors command
==================

Syntax
""""""

.. parsed-literal::

   processors Px Py Pz keyword args ...

* Px,Py,Pz = # of processors in each dimension of 3d grid overlaying the simulation domain
* zero or more keyword/arg pairs may be appended
* keyword = *grid* or *map* or *part* or *file*

  .. parsed-literal::

       *grid* arg = gstyle params ...
         gstyle = *onelevel* or *twolevel* or *numa* or *custom*
           onelevel params = none
           twolevel params = Nc Cx Cy Cz
             Nc = number of cores per node
             Cx,Cy,Cz = # of cores in each dimension of 3d sub-grid assigned to each node
           numa params = none
           custom params = infile
             infile = file containing grid layout
       *map* arg = *cart* or *cart/reorder* or *xyz* or *xzy* or *yxz* or *yzx* or *zxy* or *zyx*
          cart = use MPI_Cart() methods to map processors to 3d grid with reorder = 0
          cart/reorder = use MPI_Cart() methods to map processors to 3d grid with reorder = 1
          xyz,xzy,yxz,yzx,zxy,zyx = map processors to 3d grid in IJK ordering
       *numa* arg = none
       *part* args = Psend Precv cstyle
         Psend = partition # (1 to Np) which will send its processor layout
         Precv = partition # (1 to Np) which will recv the processor layout
         cstyle = *multiple*
           *multiple* = Psend grid will be multiple of Precv grid in each dimension
       *file* arg = outfile
         outfile = name of file to write 3d grid of processors to

Examples
""""""""

.. code-block:: LAMMPS

   processors * * 5
   processors 2 4 4
   processors * * 8 map xyz
   processors * * * grid numa
   processors * * * grid twolevel 4 * * 1
   processors 4 8 16 grid custom myfile
   processors * * * part 1 2 multiple

Description
"""""""""""

Specify how processors are mapped as a regular 3d grid to the global
simulation box.  The mapping involves 2 steps.  First if there are P
processors it means choosing a factorization P = Px by Py by Pz so
that there are Px processors in the x dimension, and similarly for the
y and z dimensions.  Second, the P processors are mapped to the
regular 3d grid.  The arguments to this command control each of these
2 steps.

The Px, Py, Pz parameters affect the factorization.  Any of the 3
parameters can be specified with an asterisk "\*", which means LAMMPS
will choose the number of processors in that dimension of the grid.
It will do this based on the size and shape of the global simulation
box so as to minimize the surface-to-volume ratio of each processor's
sub-domain.

Choosing explicit values for Px or Py or Pz can be used to override
the default manner in which LAMMPS will create the regular 3d grid of
processors, if it is known to be sub-optimal for a particular problem.
E.g. a problem where the extent of atoms will change dramatically in a
particular dimension over the course of the simulation.

The product of Px, Py, Pz must equal P, the total # of processors
LAMMPS is running on.  For a :doc:`2d simulation <dimension>`, Pz must
equal 1.

Note that if you run on a prime number of processors P, then a grid
such as 1 x P x 1 will be required, which may incur extra
communication costs due to the high surface area of each processor's
sub-domain.

Also note that if multiple partitions are being used then P is the
number of processors in this partition; see the :doc:`-partition command-line switch <Run_options>` page for details.  Also note
that you can prefix the processors command with the
:doc:`partition <partition>` command to easily specify different
Px,Py,Pz values for different partitions.

You can use the :doc:`partition <partition>` command to specify
different processor grids for different partitions, e.g.

.. code-block:: LAMMPS

   partition yes 1 processors 4 4 4
   partition yes 2 processors 2 3 2

.. note::

   This command only affects the initial regular 3d grid created
   when the simulation box is first specified via a
   :doc:`create_box <create_box>` or :doc:`read_data <read_data>` or
   :doc:`read_restart <read_restart>` command.  Or if the simulation box is
   re-created via the :doc:`replicate <replicate>` command.  The same
   regular grid is initially created, regardless of which
   :doc:`comm_style <comm_style>` command is in effect.

If load-balancing is never invoked via the :doc:`balance <balance>` or
:doc:`fix balance <fix_balance>` commands, then the initial regular grid
will persist for all simulations.  If balancing is performed, some of
the methods invoked by those commands retain the logical topology of
the initial 3d grid, and the mapping of processors to the grid
specified by the processors command.  However the grid spacings in
different dimensions may change, so that processors own sub-domains of
different sizes.  If the :doc:`comm_style tiled <comm_style>` command is
used, methods invoked by the balancing commands may discard the 3d
grid of processors and tile the simulation domain with sub-domains of
different sizes and shapes which no longer have a logical 3d
connectivity.  If that occurs, all the information specified by the
processors command is ignored.

----------

The *grid* keyword affects the factorization of P into Px,Py,Pz and it
can also affect how the P processor IDs are mapped to the 3d grid of
processors.

The *onelevel* style creates a 3d grid that is compatible with the
Px,Py,Pz settings, and which minimizes the surface-to-volume ratio of
each processor's sub-domain, as described above.  The mapping of
processors to the grid is determined by the *map* keyword setting.

The *twolevel* style can be used on machines with multicore nodes to
minimize off-node communication.  It insures that contiguous
sub-sections of the 3d grid are assigned to all the cores of a node.
For example if *Nc* is 4, then 2x2x1 or 2x1x2 or 1x2x2 sub-sections of
the 3d grid will correspond to the cores of each node.  This affects
both the factorization and mapping steps.

The *Cx*, *Cy*, *Cz* settings are similar to the *Px*, *Py*, *Pz*
settings, only their product should equal *Nc*\ .  Any of the 3
parameters can be specified with an asterisk "\*", which means LAMMPS
will choose the number of cores in that dimension of the node's
sub-grid.  As with Px,Py,Pz, it will do this based on the size and
shape of the global simulation box so as to minimize the
surface-to-volume ratio of each processor's sub-domain.

.. note::

   For the *twolevel* style to work correctly, it assumes the MPI
   ranks of processors LAMMPS is running on are ordered by core and then
   by node.  E.g. if you are running on 2 quad-core nodes, for a total of
   8 processors, then it assumes processors 0,1,2,3 are on node 1, and
   processors 4,5,6,7 are on node 2.  This is the default rank ordering
   for most MPI implementations, but some MPIs provide options for this
   ordering, e.g. via environment variable settings.

The *numa* style operates similar to the *twolevel* keyword except
that it auto-detects which cores are running on which nodes.
Currently, it does this in only 2 levels, but it may be extended in
the future to account for socket topology and other non-uniform memory
access (NUMA) costs.  It also uses a different algorithm than the
*twolevel* keyword for doing the two-level factorization of the
simulation box into a 3d processor grid to minimize off-node
communication, and it does its own MPI-based mapping of nodes and
cores to the regular 3d grid.  Thus it may produce a different layout
of the processors than the *twolevel* options.

The *numa* style will give an error if the number of MPI processes is
not divisible by the number of cores used per node, or any of the Px
or Py of Pz values is greater than 1.

.. note::

   Unlike the *twolevel* style, the *numa* style does not require
   any particular ordering of MPI ranks i norder to work correctly.  This
   is because it auto-detects which processes are running on which nodes.

The *custom* style uses the file *infile* to define both the 3d
factorization and the mapping of processors to the grid.

The file should have the following format.  Any number of initial
blank or comment lines (starting with a "#" character) can be present.
The first non-blank, non-comment line should have
3 values:

.. parsed-literal::

   Px Py Py

These must be compatible with the total number of processors
and the Px, Py, Pz settings of the processors command.

This line should be immediately followed by
P = Px\*Py\*Pz lines of the form:

.. parsed-literal::

   ID I J K

where ID is a processor ID (from 0 to P-1) and I,J,K are the
processors location in the 3d grid.  I must be a number from 1 to Px
(inclusive) and similarly for J and K.  The P lines can be listed in
any order, but no processor ID should appear more than once.

----------

The *map* keyword affects how the P processor IDs (from 0 to P-1) are
mapped to the 3d grid of processors.  It is only used by the
*onelevel* and *twolevel* grid settings.

The *cart* style uses the family of MPI Cartesian functions to perform
the mapping, namely MPI_Cart_create(), MPI_Cart_get(),
MPI_Cart_shift(), and MPI_Cart_rank().  It invokes the
MPI_Cart_create() function with its reorder flag = 0, so that MPI is
not free to reorder the processors.

The *cart/reorder* style does the same thing as the *cart* style
except it sets the reorder flag to 1, so that MPI can reorder
processors if it desires.

The *xyz*, *xzy*, *yxz*, *yzx*, *zxy*, and *zyx* styles are all
similar.  If the style is IJK, then it maps the P processors to the
grid so that the processor ID in the I direction varies fastest, the
processor ID in the J direction varies next fastest, and the processor
ID in the K direction varies slowest.  For example, if you select
style *xyz* and you have a 2x2x2 grid of 8 processors, the assignments
of the 8 octants of the simulation domain will be:

.. parsed-literal::

   proc 0 = lo x, lo y, lo z octant
   proc 1 = hi x, lo y, lo z octant
   proc 2 = lo x, hi y, lo z octant
   proc 3 = hi x, hi y, lo z octant
   proc 4 = lo x, lo y, hi z octant
   proc 5 = hi x, lo y, hi z octant
   proc 6 = lo x, hi y, hi z octant
   proc 7 = hi x, hi y, hi z octant

Note that, in principle, an MPI implementation on a particular machine
should be aware of both the machine's network topology and the
specific subset of processors and nodes that were assigned to your
simulation.  Thus its MPI_Cart calls can optimize the assignment of
MPI processes to the 3d grid to minimize communication costs.  In
practice, however, few if any MPI implementations actually do this.
So it is likely that the *cart* and *cart/reorder* styles simply give
the same result as one of the IJK styles.

Also note, that for the *twolevel* grid style, the *map* setting is
used to first map the nodes to the 3d grid, then again to the cores
within each node.  For the latter step, the *cart* and *cart/reorder*
styles are not supported, so an *xyz* style is used in their place.

----------

The *part* keyword affects the factorization of P into Px,Py,Pz.

It can be useful when running in multi-partition mode, e.g. with the
:doc:`run_style verlet/split <run_style>` command.  It specifies a
dependency between a sending partition *Psend* and a receiving
partition *Precv* which is enforced when each is setting up their own
mapping of their processors to the simulation box.  Each of *Psend*
and *Precv* must be integers from 1 to Np, where Np is the number of
partitions you have defined via the :doc:`-partition command-line switch <Run_options>`.

A "dependency" means that the sending partition will create its
regular 3d grid as Px by Py by Pz and after it has done this, it will
send the Px,Py,Pz values to the receiving partition.  The receiving
partition will wait to receive these values before creating its own
regular 3d grid and will use the sender's Px,Py,Pz values as a
constraint.  The nature of the constraint is determined by the
*cstyle* argument.

For a *cstyle* of *multiple*, each dimension of the sender's processor
grid is required to be an integer multiple of the corresponding
dimension in the receiver's processor grid.  This is a requirement of
the :doc:`run_style verlet/split <run_style>` command.

For example, assume the sending partition creates a 4x6x10 grid = 240
processor grid.  If the receiving partition is running on 80
processors, it could create a 4x2x10 grid, but it will not create a
2x4x10 grid, since in the y-dimension, 6 is not an integer multiple of
4.

.. note::

   If you use the :doc:`partition <partition>` command to invoke
   different "processors" commands on different partitions, and you also
   use the *part* keyword, then you must insure that both the sending and
   receiving partitions invoke the "processors" command that connects the
   2 partitions via the *part* keyword.  LAMMPS cannot easily check for
   this, but your simulation will likely hang in its setup phase if this
   error has been made.

----------

The *file* keyword writes the mapping of the factorization of P
processors and their mapping to the 3d grid to the specified file
*outfile*\ .  This is useful to check that you assigned physical
processors in the manner you desired, which can be tricky to figure
out, especially when running on multiple partitions or on, a multicore
machine or when the processor ranks were reordered by use of the
:doc:`-reorder command-line switch <Run_options>` or due to use of
MPI-specific launch options such as a config file.

If you have multiple partitions you should insure that each one writes
to a different file, e.g. using a :doc:`world-style variable <variable>`
for the filename.  The file has a self-explanatory header, followed by
one-line per processor in this format:

world-ID universe-ID original-ID: I J K: name

The IDs are the processor's rank in this simulation (the world), the
universe (of multiple simulations), and the original MPI communicator
used to instantiate LAMMPS, respectively.  The world and universe IDs
will only be different if you are running on more than one partition;
see the :doc:`-partition command-line switch <Run_options>`.  The
universe and original IDs will only be different if you used the
:doc:`-reorder command-line switch <Run_options>` to reorder the
processors differently than their rank in the original communicator
LAMMPS was instantiated with.

I,J,K are the indices of the processor in the regular 3d grid, each
from 1 to Nd, where Nd is the number of processors in that dimension
of the grid.

The *name* is what is returned by a call to MPI_Get_processor_name()
and should represent an identifier relevant to the physical processors
in your machine.  Note that depending on the MPI implementation,
multiple cores can have the same *name*\ .

----------

Restrictions
""""""""""""

This command cannot be used after the simulation box is defined by a
:doc:`read_data <read_data>` or :doc:`create_box <create_box>` command.
It can be used before a restart file is read to change the 3d
processor grid from what is specified in the restart file.

The *grid numa* keyword only currently works with the *map cart*
option.

The *part* keyword (for the receiving partition) only works with the
*grid onelevel* or *grid twolevel* options.

Related commands
""""""""""""""""

:doc:`partition <partition>`, :doc:`-reorder command-line switch <Run_options>`

Default
"""""""

The option defaults are Px Py Pz = \* \* \*, grid = onelevel, and map =
cart.
