.. index:: fix eos/cv

fix eos/cv command
==================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID eos/cv cv

* ID, group-ID are documented in :doc:`fix <fix>` command
* eos/cv = style name of this fix command
* cv = constant-volume heat capacity (energy/temperature units)

Examples
""""""""


.. parsed-literal::

   fix 1 all eos/cv 0.01

Description
"""""""""""

Fix *eos/cv* applies a mesoparticle equation of state to relate the
particle internal energy (u\_i) to the particle internal temperature
(dpdTheta\_i).  The *eos/cv* mesoparticle equation of state requires
the constant-volume heat capacity, and is defined as follows:

.. image:: Eqs/fix_eos-cv.jpg
   :align: center

where Cv is the constant-volume heat capacity, u\_cond is the internal
conductive energy, and u\_mech is the internal mechanical energy.  Note
that alternative definitions of the mesoparticle equation of state are
possible.


----------


Restrictions
""""""""""""


This command is part of the USER-DPD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This command also requires use of the :doc:`atom_style dpd <atom_style>`
command.

Related commands
""""""""""""""""

:doc:`fix shardlow <fix_shardlow>`, :doc:`pair dpd/fdt <pair_dpd_fdt>`

**Default:** none


----------


.. _Larentzos4:



**(Larentzos)** J.P. Larentzos, J.K. Brennan, J.D. Moore, and
W.D. Mattson, "LAMMPS Implementation of Constant Energy Dissipative
Particle Dynamics (DPD-E)", ARL-TR-6863, U.S. Army Research
Laboratory, Aberdeen Proving Ground, MD (2014).


