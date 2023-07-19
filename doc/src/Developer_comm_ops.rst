Communication patterns
----------------------

This page describes various inter-processor communication operations
provided by LAMMPS, mostly in the core *Comm* class.  These are operations
for common tasks implemented using MPI library calls.  They are used by
other classes to perform communication of different kinds.  These
operations are useful to know about when writing new code for LAMMPS
that needs to communicate data between processors.

Owned and ghost atoms
^^^^^^^^^^^^^^^^^^^^^

As described on the :doc:`parallel partitioning algorithms
<Developer_par_part>` page, LAMMPS spatially decomposes the simulation
domain, either in a *brick* or *tiled* manner.  Each processor (MPI
task) owns atoms within its subdomain and additionally stores ghost
atoms within a cutoff distance of its subdomain.

Forward and reverse communication
=================================

As described on the :doc:`parallel communication algorithms
<Developer_par_comm>` page, the most common communication operations are
first, *forward communication* which sends owned atom information from
each processor to nearby processors to store with their ghost atoms.
The need to do this communication arises when data from the owned atoms
is updated (e.g. their positions) and this updated information needs to
be **copied** to the corresponding ghost atoms.

And second, *reverse communication*, which sends ghost atom information
from each processor to the owning processor to **accumulate** (sum)
the values with the corresponding owned atoms.  The need for this
arises when data is computed and also stored with ghost atoms
(e.g. forces when using a "half" neighbor list) and thus those terms
need to be added to their corresponding atoms on the process where
they are "owned" atoms.  Please note, that with the :doc:`newton off
<newton>` setting this does not happen and the neighbor lists are
constructed so that these interactions are computed on both MPI
processes containing one of the atoms and only the data pertaining to
the local atom is stored.

The time-integration classes in LAMMPS invoke these operations each
timestep via the *forward_comm()* and *reverse_comm()* methods in the
*Comm* class.  Which per-atom data is communicated depends on the
currently used :doc:`atom style <atom_style>` and whether
:doc:`comm_modify vel <comm_modify>` setting is "no" (default) or
"yes".

Similarly, *Pair* style classes can invoke the *forward_comm(this)*
and *reverse_comm(this)* methods in the *Comm* class to perform the
same operations on per-atom data that is generated and stored within
the pair style class. Note that this function requires passing the
``this`` pointer as the first argument to enable the *Comm* class to
call the "pack" and "unpack" functions discussed below.  An example of
the use of these functions are many-body pair styles like the
embedded-atom method (EAM) which compute intermediate values in the
first part of the compute() function that need to be stored by both
owned and ghost atoms for the second part of the force computation.
The *Comm* class methods perform the MPI communication for buffers of
per-atom data.  They "call back" to the *Pair* class, so it can *pack*
or *unpack* the buffer with data the *Pair* class owns.  There are 4
such methods that the *Pair* class must define, assuming it uses both
forward and reverse communication:

* pack_forward_comm()
* unpack_forward_comm()
* pack_reverse_comm()
* unpack_reverse_comm()

The arguments to these methods include the buffer and a list of atoms
to pack or unpack.  The *Pair* class also must set the *comm_forward*
and *comm_reverse* variables, which store the number of values stored
in the communication buffers for each operation.  This means, if
desired, it can choose to store multiple per-atom values in the
buffer, and they will be communicated together to minimize
communication overhead.  The communication buffers are defined vectors
containing ``double`` values.  To correctly store integers that may be
64-bit (bigint, tagint, imageint) in the buffer, you need to use the
:ref:`ubuf union <communication_buffer_coding_with_ubuf>` construct.

The *Fix*, *Compute*, and *Dump* classes can also invoke the same kind
of forward and reverse communication operations using the same *Comm*
class methods.  Likewise, the same pack/unpack methods and
comm_forward/comm_reverse variables must be defined by the calling
*Fix*, *Compute*, or *Dump* class.

For *Fix* classes, there is an optional second argument to the
*forward_comm()* and *reverse_comm()* call which can be used when the
fix performs multiple modes of communication, with different numbers
of values per atom.  The fix should set the *comm_forward* and
*comm_reverse* variables to the maximum value, but can invoke the
communication for a particular mode with a smaller value.  For this
to work, the *pack_forward_comm()*, etc methods typically use a class
member variable to choose which values to pack/unpack into/from the
buffer.

Finally, for reverse communications in *Fix* classes there is also the
*reverse_comm_variable()* method that allows the communication to have
a different amount of data per-atom.  It invokes these corresponding
callback methods:

* pack_reverse_comm_size()
* unpack_reverse_comm_size()

which have extra arguments to specify the amount of data stored
in the buffer for each atom.

Higher level communication
^^^^^^^^^^^^^^^^^^^^^^^^^^

There are also several higher-level communication operations provided
in LAMMPS which work for either *brick* or *tiled* decompositions.
They may be useful for a new class to invoke if it requires more
sophisticated communication than the *forward* and *reverse* methods
provide.  The 3 communication operations described here are

* ring
* irregular
* rendezvous

You can invoke these *grep* command in the LAMMPS src directory, to
see a list of classes that invoke the 3 operations.

* ``grep "\->ring" *.cpp */*.cpp``
* ``grep "irregular\->" *.cpp``
* ``grep "\->rendezvous" *.cpp */*.cpp``

Ring operation
==============

The *ring* operation is invoked via the *ring()* method in the *Comm*
class.

Each processor first creates a buffer with a list of values, typically
associated with a subset of the atoms it owns.  Now think of the *P*
processors as connected to each other in a *ring*.  Each processor *M*
sends data to the next *M+1* processor.  It receives data from the
preceding *M-1* processor.  The ring is periodic so that the last
processor sends to the first processor, and the first processor
receives from the last processor.

Invoking the *ring()* method passes each processor's buffer in *P*
steps around the ring.  At each step a *callback* method, provided as
an argument to ring(), in the caller is invoked.  This allows each
processor to examine the data buffer provided by every other
processor.  It may extract values needed by its atoms from the
buffers, or it may alter placeholder values in the buffer.  In the
latter case, when the *ring* operation is complete, each processor can
examine its original buffer to extract modified values.

Note that the *ring* operation is similar to an MPI_Alltoall()
operation, where every processor effectively sends and receives data to
every other processor.  The difference is that the *ring* operation
does it one step at a time, so the total volume of data does not need
to be stored by every processor.  However, the *ring* operation is
also less efficient than MPI_Alltoall() because of the *P* stages
required.  So it is typically only suitable for small data buffers and
occasional operations that are not time-critical.

Irregular operation
===================

The *irregular* operation is provided by the *Irregular* class.  What
LAMMPS terms irregular communication is when each processor knows what
data it needs to send to what processor, but does not know what
processors are sending it data.  An example is when load-balancing is
performed and each processor needs to send some of its atoms to new
processors.

The *Irregular* class provides 5 high-level methods useful in this
context:

* create_data()
* exchange_data()
* create_atom()
* exchange_atom()
* migrate_atoms()

For the *create_data()* method, each processor specifies a list of *N*
datums to send, each to a specified processor.  Internally, the method
creates efficient data structures for performing the communication.
The *exchange_data()* method triggers the communication to be
performed.  Each processor provides the vector of *N* datums to send,
and the size of each datum.  All datums must be the same size.

The *create_atom()* and *exchange_atom()* methods are similar, except
that the size of each datum can be different.  Typically, this is used
to communicate atoms, each with a variable amount of per-atom data, to
other processors.

The *migrate_atoms()* method is a convenience wrapper on the
*create_atom()* and *exchange_atom()* methods to simplify
communication of all the per-atom data associated with an atom so that
the atom can effectively migrate to a new owning processor.  It is
similar to the *exchange()* method in the *Comm* class invoked when
atoms move to neighboring processors (in the regular or tiled
decomposition) during timestepping, except that it allows atoms to
have moved arbitrarily long distances and still be properly
communicated to a new owning processor.

Rendezvous operation
====================

Finally, the *rendezvous* operation is invoked via the *rendezvous()*
method in the *Comm* class.  Depending on how much communication is
needed and how many processors a LAMMPS simulation is running on, it
can be a much more efficient choice than the *ring()* method.  It uses
the *irregular* operation internally once or twice to do its
communication.  The rendezvous algorithm is described in detail in
:ref:`(Plimpton) <Plimpton>`, including some LAMMPS use cases.

For the *rendezvous()* method, each processor specifies a list of *N*
datums to send and which processor to send each of them to.
Internally, this communication is performed as an irregular operation.
The received datums are returned to the caller via invocation of
*callback* function, provided as an argument to *rendezvous()*.  The
caller can then process the received datums and (optionally) assemble
a new list of datums to communicate to a new list of specific
processors.  When the callback function exits, the *rendezvous()*
method performs a second irregular communication on the new list of
datums.

Examples in LAMMPS of use of the *rendezvous* operation are the
:doc:`fix rigid/small <fix_rigid>` and :doc:`fix shake
<fix_shake>` commands (for one-time identification of the rigid body
atom clusters) and the identification of special_bond 1-2, 1-3 and 1-4
neighbors within molecules.  See the :doc:`special_bonds <special_bonds>`
command for context.

----------

.. _Plimpton:

**(Plimpton)** Plimpton and Knight, JPDC, 147, 184-195 (2021).
