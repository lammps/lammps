.. index:: dump atom/adios
.. index:: dump custom/adios
.. index:: dump local/adios

dump atom/adios  command
=========================

dump custom/adios command
=========================

dump local/adios command
========================

Syntax
""""""

.. code-block:: LAMMPS

   dump ID group-ID atom/adios N file.bp
   dump ID group-ID custom/adios N file.bp args
   dump ID group-ID local/adios N file.bp args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be dumped
* N = dump every this many timesteps
* file.bp = name of file/stream to write to
* args for ``custom/adios`` = same options as in :doc:`dump custom <dump>` command
* args for ``local/adios`` = same options as in :doc:`dump local <dump>` command

Examples
""""""""

.. code-block:: LAMMPS

   dump adios1 all atom/adios   100 atoms.bp
   dump 4a     all custom/adios 100 dump_adios.bp id v_p x y z
   dump 2 subgroup custom/adios 100 dump_adios.bp mass type xs ys zs vx vy vz
   dump bonds  all local/adios  100 bonds.bp c_myBondCompute[1] c_myBondCompute[2]

   # write trajectory data in single-precision floating point
   dump adios1 all custom/adios 100 dump.bp id type x y z
   dump_modify adios1 precision fp32

   # four replicas all writing trajectory data to one shared .bp file
   dump adios1 all custom/adios 100 shared.bp id type x y z
   dump_modify adios1 shared yes

   # local data + shared replicas + fp32
   dump bonds all local/adios 100 bonds.bp c_myBondCompute[1] c_myBondCompute[2]
   dump_modify bonds shared yes precision fp32

Description
"""""""""""

Dump a snapshot of atom coordinates (``atom/adios``), per-atom quantities
(``custom/adios``), or per-local-entity quantities (``local/adios``) every
:math:`N` timesteps in the `ADIOS <adios_>`_-based "BP" file format.
ADIOS-BP files are binary, portable, and self-describing.

.. _adios: https://github.com/ornladios/ADIOS2

.. note::

   To be able to use ADIOS, a file ``adios2_config.xml`` with specific
   configuration settings is expected in the current working directory.
   If the file is not present, LAMMPS will try to create a minimal
   default file.  Please refer to the ADIOS documentation for details on
   how to adjust this file for optimal performance and desired features.

**Precision control (custom/adios, local/adios):**

By default all numeric data (trajectory table, box bounds) is written in
64-bit double precision (fp64).  Use :doc:`dump_modify <dump_modify>` with
the ``precision`` keyword to switch to 32-bit single precision (fp32):

.. code-block:: LAMMPS

   dump_modify ID precision fp32

Valid values are ``fp32`` (or ``float``) and ``fp64`` (or ``double``).
The chosen precision is recorded as the ``LAMMPS/precision`` attribute
(``"float32"`` or ``"float64"``) in the BP file metadata so that downstream
readers can detect it automatically.

Box bound variables (``boxxlo``, ``boxxhi``, etc.) are also written at the
selected precision.  Integer metadata (``ntimestep``, ``nprocs``,
``natoms``, etc.) is unaffected.

**Multi-replica shared output (custom/adios, local/adios):**

When running LAMMPS with multiple replicas (the ``-partition`` command-line
flag), each replica independently writes its own BP file by default.  The
``shared`` keyword causes all replicas to write collectively
to a *single* BP file:

.. code-block:: LAMMPS

   dump_modify ID shared yes

In shared mode the BP file contains one ADIOS step per global timestep.
Each replica writes to its own per-replica ADIOS variable so that data from
different replicas can be distinguished.  Variables are named with a
``_replicaN`` suffix (e.g. ``atoms_replica0``, ``atoms_replica1``,
``boxxlo_replica0``, ``boxxlo_replica1``, …).  The ``nreplicas`` attribute
records the total number of replicas.

The ``ntimestep`` and ``nprocs`` global scalars are written only once per
step (by replica 0) and carry no suffix.

This feature requires that all replicas hit the ``dump`` call at the same
timestep, which is automatically satisfied in standard partition runs.

**dump local/adios specifics:**

``dump local/adios`` writes local quantities computed by fixes or computes
(e.g. bond data, pair-data) in the same BP5 format as ``custom/adios``.
The local-data table is stored in a variable named ``local_data`` (or
``local_data_replicaN`` in shared mode).  Box bounds and step metadata are
included for convenience.  The ``LAMMPS/dump_style`` attribute is set to
``"local"``.

**Use from write_dump:**

It is possible to use these dump styles with the
:doc:`write_dump <write_dump>` command.  In this case, the sub-intervals
must not be set at all.  The write_dump command can be used to
create a new file at each individual dump.

.. code-block:: LAMMPS

   dump 4     all atom/adios 100 dump.bp
   write_dump all atom/adios singledump.bp

----------

Restrictions
""""""""""""

The number of atoms per snapshot **can** change with the adios style.
When using the ADIOS tool ``bpls`` to list the content of a .bp file,
``bpls`` will print ``*__*`` for the size of the output table indicating
that its size is changing every step.

The ``precision`` and ``shared`` keywords must be set (via ``dump_modify``)
**before the first run** that activates the dump.  Changing these options
after the ADIOS IO group has been initialised has no effect.

The ``shared yes`` option has no effect when running with a single replica
(``universe->nworlds == 1``).

The ``atom/adios``, ``custom/adios``, and ``local/adios`` dump styles are
part of the ADIOS package.  They are only enabled if LAMMPS was built with
that package.  See the :doc:`Build package <Build_package>` page for more
info.

----------

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`
