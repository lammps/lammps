.. index:: dump atom/adios
.. index:: dump custom/adios

dump atom/adios  command
=========================

dump custom/adios command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   dump ID group-ID atom/adios N file.bp
   dump ID group-ID custom/adios N file.bp args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be imaged
* adios = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file.bp = name of file/stream to write to
* args = same options as in :doc:`dump custom <dump>` command

Examples
""""""""

.. code-block:: LAMMPS

   dump adios1 all atom/adios   100 atoms.bp
   dump 4a     all custom/adios 100 dump_adios.bp id v_p x y z
   dump 2 subgroup custom/adios 100 dump_adios.bp mass type xs ys zs vx vy vz

Description
"""""""""""

Dump a snapshot of atom coordinates every :math:`N` timesteps in the `ADIOS
<adios_>`_-based "BP" file format, or using different I/O solutions in
ADIOS, to a stream that can be read on-line by another program.
ADIOS-BP files are binary, portable, and self-describing.

.. _adios: https://github.com/ornladios/ADIOS2

.. note::

   To be able to use ADIOS, a file ``adios2_config.xml`` with specific
   configuration settings is expected in the current working directory.
   If the file is not present, LAMMPS will try to create a minimal
   default file.  Please refer to the ADIOS documentation for details on
   how to adjust this file for optimal performance and desired features.

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
When using the ADIOS tool 'bpls' to list the content of a .bp file,
bpls will print *__* for the size of the output table indicating that
its size is changing every step.

The *atom/adios* and *custom/adios* dump styles are part of the ADIOS
package.  They are only enabled if LAMMPS was built with that package.
See the :doc:`Build package <Build_package>` page for more info.

----------

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`
