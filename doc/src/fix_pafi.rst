.. index:: fix pafi

fix pafi command
================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID pafi compute-ID Temp Tdamp seed overdamped <arg> com <arg>

* ID, group-ID are documented in :doc:`fix <fix>` command
* pafi = style name of this fix command
* compute-ID = ID of a :doc:`compute property/atom <compute_property_atom>` that holds data used by this fix
* Temp = desired temperature (temperature units)
* Tdamp = damping parameter (time units)
* seed = random number seed to use for white noise (positive integer)
* overdamped <arg> = enable (*arg* = 1) or disable (*arg* = 0) XXX
* com <arg> =  enable (*arg* = 1) or disable (*arg* = 0) XXX

Examples
""""""""

.. code-block:: LAMMPS

   compute pa all property/atom d_nx d_ny d_nz d_dnx d_dny d_dnz d_ddnx d_ddny d_ddnz
   run 0 post no
   fix hp all pafi pa 500.0 0.01 434 overdamped 0 com 1

Description
"""""""""""

TODO: add docs

**Restart, fix_modify, output, run start/stop, minimize info:**

Restrictions
""""""""""""

This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

**Default:** none
