.. index:: compute spin

compute spin command
====================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID spin

* ID, group-ID are documented in :doc:`compute <compute>` command
* spin = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute out_mag all spin

Description
"""""""""""

Define a computation that calculates magnetic quantities for a system
of atoms having spins.

This compute calculates the following 6 magnetic quantities:

* the three first quantities are the x,y and z coordinates of the total
  magnetization,
* the fourth quantity is the norm of the total magnetization,
* The fifth quantity is the magnetic energy (in eV),
* The sixth one is referred to as the spin temperature, according
  to the work of :ref:`(Nurdin) <Nurdin1>`.

The simplest way to output the results of the compute spin calculation
is to define some of the quantities as variables, and to use the thermo and
thermo_style commands, for example:

.. code-block:: LAMMPS

   compute out_mag         all spin

   variable mag_z          equal c_out_mag[3]
   variable mag_norm       equal c_out_mag[4]
   variable temp_mag       equal c_out_mag[6]

   thermo                  10
   thermo_style            custom step v_mag_z v_mag_norm v_temp_mag

This series of commands evaluates the total magnetization along z, the norm of
the total magnetization, and the magnetic temperature. Three variables are
assigned to those quantities. The thermo and thermo_style commands print them
every 10 timesteps.

Output info
"""""""""""

The array values are "intensive".  The array values will be in
metal units (:doc:`units <units>`).

Restrictions
""""""""""""

The *spin* compute is part of the SPIN package.  This compute is only
enabled if LAMMPS was built with this package.  See the :doc:`Build package <Build_package>` page for more info.  The atom_style
has to be "spin" for this compute to be valid.

**Related commands:**

none

Default
"""""""


none

----------

.. _Nurdin1:

**(Nurdin)** Nurdin and Schotte Phys Rev E, 61(4), 3579 (2000)
