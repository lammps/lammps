Removed commands and packages
=============================

This page lists LAMMPS commands and packages that have been removed from
the distribution and provides suggestions for alternatives or
replacements.  LAMMPS has special dummy styles implemented, that will
stop LAMMPS and print a suitable error message in most cases, when a
style/command is used that has been removed or will replace the command
with the direct alternative (if available) and print a warning.

Fix ave/spatial and fix ave/spatial/sphere
------------------------------------------

.. deprecated:: 11Dec2015

The fixes ave/spatial and ave/spatial/sphere have been removed from LAMMPS
since they were superseded by the more general and extensible "chunk
infrastructure".  Here the system is partitioned in one of many possible
ways through the :doc:`compute chunk/atom <compute_chunk_atom>` command
and then averaging is done using :doc:`fix ave/chunk <fix_ave_chunk>`.
Please refer to the :doc:`chunk HOWTO <Howto_chunk>` section for an overview.

Box command
-----------

.. deprecated:: 22Dec2022

The *box* command has been removed and the LAMMPS code changed so it won't
be needed.  If present, LAMMPS will ignore the command and print a warning.

Reset_ids, reset_atom_ids, reset_mol_ids commands
-------------------------------------------------

.. deprecated:: 22Dec2022

The *reset_ids*, *reset_atom_ids*, and *reset_mol_ids* commands have
been folded into the :doc:`reset_atoms <reset_atoms>` command.  If
present, LAMMPS will replace the commands accordingly and print a
warning.

LATTE package
-------------

.. deprecated:: 15Jun2023

The LATTE package with the fix latte command was removed from LAMMPS.
This functionality has been superseded by :doc:`fix mdi/qm <fix_mdi_qm>`
and :doc:`fix mdi/qmmm <fix_mdi_qmmm>` from the :ref:`MDI package
<PKG-MDI>`.  These fixes are compatible with several quantum software
packages, including LATTE.  See the ``examples/QUANTUM`` dir and the
:doc:`MDI coupling HOWTO <Howto_mdi>` page.  MDI supports running LAMMPS
with LATTE as a plugin library (similar to the way fix latte worked), as
well as on a different set of MPI processors.

MEAM package
------------

The MEAM package in Fortran has been replaced by a C++ implementation.
The code in the :ref:`MEAM package <PKG-MEAM>` is a translation of the
Fortran code of MEAM into C++, which removes several restrictions
(e.g. there can be multiple instances in hybrid pair styles) and allows
for some optimizations leading to better performance.  The pair style
:doc:`meam <pair_meam>` has the exact same syntax.  For a transition
period the C++ version of MEAM was called USER-MEAMC so it could
coexist with the Fortran version.

Minimize style fire/old
-----------------------

.. deprecated:: 8Feb2023

Minimize style *fire/old* has been removed. Its functionality can be
reproduced with *fire* with specific options. Please see the
:doc:`min_modify command <min_modify>` documentation for details.

Pair style mesont/tpm, compute style mesont, atom style mesont
--------------------------------------------------------------

.. deprecated:: 8Feb2023

Pair style *mesont/tpm*, compute style *mesont*, and atom style
*mesont* have been removed from the :ref:`MESONT package <PKG-MESONT>`.
The same functionality is available through
:doc:`pair style mesocnt <pair_mesocnt>`,
:doc:`bond style mesocnt <bond_mesocnt>` and
:doc:`angle style mesocnt <angle_mesocnt>`.

MPIIO package
-------------

.. deprecated:: 21Nov2023

The MPIIO package has been removed from LAMMPS since it was unmaintained
for many years and thus not updated to incorporate required changes that
had been applied to the corresponding non-MPIIO commands. As a
consequence the MPIIO commands had become unreliable and sometimes
crashing LAMMPS or corrupting data.  Similar functionality is available
through the :ref:`ADIOS package <PKG-ADIOS>` and the :ref:`NETCDF
package <PKG-NETCDF>`.  Also, the :doc:`dump_modify nfile or dump_modify
fileper <dump_modify>` keywords may be used for an efficient way of
writing out dump files when running on large numbers of processors.
Similarly, the "nfile" and "fileper" keywords exist for restarts:
see :doc:`restart <restart>`, :doc:`read_restart <read_restart>`,
:doc:`write_restart <write_restart>`.


MSCG package
------------

.. deprecated:: 21Nov2023

The MSCG package has been removed from LAMMPS since it was unmaintained
for many years and instead superseded by the `OpenMSCG software
<https://software.rcc.uchicago.edu/mscg/>`_ of the Voth group at the
University of Chicago, which can be used independent from LAMMPS.

REAX package
------------

The REAX package has been removed since it was superseded by the
:ref:`REAXFF package <PKG-REAXFF>`.  The REAXFF package has been tested
to yield equivalent results to the REAX package, offers better
performance, supports OpenMP multi-threading via OPENMP, and GPU and
threading parallelization through KOKKOS.  The new pair styles are not
syntax compatible with the removed reax pair style, so input files will
have to be adapted.  The REAXFF package was originally called
USER-REAXC.

USER-REAXC package
------------------

.. deprecated:: 7Feb2024

The USER-REAXC package has been renamed to :ref:`REAXFF <PKG-REAXFF>`.
In the process also the pair style and related fixes were renamed to use
the "reaxff" string instead of "reax/c". For a while LAMMPS was maintaining
backward compatibility by providing aliases for the styles.  These have
been removed, so using "reaxff" is now *required*.

USER-CUDA package
-----------------

The USER-CUDA package had been removed, since it had been unmaintained
for a long time and had known bugs and problems.  Significant parts of
the design were transferred to the
:ref:`KOKKOS package <PKG-KOKKOS>`, which has similar
performance characteristics on NVIDIA GPUs. Both, the KOKKOS
and the :ref:`GPU package <PKG-GPU>` are maintained
and allow running LAMMPS with GPU acceleration.

i-PI tool
---------

.. versionchanged:: TBD

The i-PI tool has been removed from the LAMMPS distribution.  Instead,
instructions to install i-PI from PyPi via pip are provided.

restart2data tool
-----------------

The functionality of the restart2data tool has been folded into the
LAMMPS executable directly instead of having a separate tool.  A
combination of the commands :doc:`read_restart <read_restart>` and
:doc:`write_data <write_data>` can be used to the same effect.  For
added convenience this conversion can also be triggered by
:doc:`command line flags <Run_options>`
