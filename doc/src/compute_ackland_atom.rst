.. index:: compute ackland/atom

compute ackland/atom command
============================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID ackland/atom keyword/value

* ID, group-ID are documented in :doc:`compute <compute>` command
* ackland/atom = style name of this compute command
* zero or more keyword/value pairs may be appended
* keyword = *legacy*

  .. parsed-literal::

       *legacy* yes/no = use (\ *yes*\ ) or do not use (\ *no*\ ) legacy ackland algorithm implementation



Examples
""""""""


.. parsed-literal::

   compute 1 all ackland/atom
   compute 1 all ackland/atom legacy yes

Description
"""""""""""

Defines a computation that calculates the local lattice structure
according to the formulation given in :ref:`(Ackland) <Ackland>`.
Historically, LAMMPS had two, slightly different implementations of
the algorithm from the paper. With the *legacy* keyword, it is
possible to switch between the pre-2015 (\ *legacy yes*\ ) and post-2015
implementation (\ *legacy no*\ ). The post-2015 variant is the default.

In contrast to the :doc:`centro-symmetry parameter <compute_centro_atom>` this method is stable against
temperature boost, because it is based not on the distance between
particles but the angles.  Therefore statistical fluctuations are
averaged out a little more.  A comparison with the Common Neighbor
Analysis metric is made in the paper.

The result is a number which is mapped to the following different
lattice structures:

* 0 = UNKNOWN
* 1 = BCC
* 2 = FCC
* 3 = HCP
* 4 = ICO

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e. each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently or to have multiple compute/dump commands, each of
which computes this quantity.-

**Output info:**

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

Restrictions
""""""""""""


This compute is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

The per-atom vector values will be unitless since they are the
integers defined above.

Related commands
""""""""""""""""

:doc:`compute centro/atom <compute_centro_atom>`

Default
"""""""
The keyword *legacy* defaults to *no*\ .


----------


.. _Ackland:



**(Ackland)** Ackland, Jones, Phys Rev B, 73, 054104 (2006).
