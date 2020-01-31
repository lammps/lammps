Removed commands and packages
=============================

This page lists LAMMPS commands and packages that have been removed from
the distribution and provides suggestions for alternatives or replacements.
LAMMPS has special dummy styles implemented, that will stop LAMMPS and
print a suitable error message in most cases, when a style/command is used
that has been removed.

Fix ave/spatial and fix ave/spatial/sphere
------------------------------------------

The fixes ave/spatial and ave/spatial/sphere have been removed from LAMMPS
since they were superseded by the more general and extensible "chunk
infrastructure".  Here the system is partitioned in one of many possible
ways through the :doc:`compute chunk/atom <compute_chunk_atom>` command
and then averaging is done using :doc:`fix ave/chunk <fix_ave_chunk>`.
Please refer to the :doc:`chunk HOWTO <Howto_chunk>` section for an overview.

MEAM package
------------

The MEAM package has been removed since it was superseded by the
:ref:`USER-MEAMC package <PKG-USER-MEAMC>`. The code in
the USER-MEAMC package is a translation of the Fortran code of MEAM into C++,
which removes several restrictions (e.g. there can be multiple instances
in hybrid pair styles) and allows for some optimizations leading
to better performance.  The new pair style :doc:`meam/c <pair_meamc>` has
the exact same syntax as the old "meam" pair style and thus pair style
:doc:`meam <pair_meamc>` is an alias to the new style and backward
compatibility of old inputs is preserved.

REAX package
------------

The REAX package has been removed since it was superseded by the
:ref:`USER-REAXC package <PKG-USER-REAXC>`.  The USER-REAXC
package has been tested to yield equivalent results to the REAX package,
offers better performance, supports OpenMP multi-threading via USER-OMP,
and GPU and threading parallelization through KOKKOS.  The new pair styles
are not syntax compatible with the removed reax pair style, so input
files will have to be adapted.

USER-CUDA package
-----------------

The USER-CUDA package had been removed, since it had been unmaintained
for a long time and had known bugs and problems.  Significant parts of
the design were transferred to the
:ref:`KOKKOS package <PKG-KOKKOS>`, which has similar
performance characteristics on Nvidia GPUs. Both, the KOKKOS
and the :ref:`GPU package <PKG-GPU>` are maintained
and allow running LAMMPS with GPU acceleration.

restart2data tool
-----------------

The functionality of the restart2data tool has been folded into the
LAMMPS executable directly instead of having a separate tool.  A
combination of the commands :doc:`read_restart <read_restart>` and
:doc:`write_data <write_data>` can be used to the same effect.  For added
convenience this conversion can also be triggered by :doc:`command line flags <Run_options>`


