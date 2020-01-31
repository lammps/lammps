.. index:: compute edpd/temp/atom

compute edpd/temp/atom command
==============================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID edpd/temp/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* edpd/temp/atom = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all edpd/temp/atom

Description
"""""""""""

Define a computation that calculates the per-atom temperature
for each eDPD particle in a group.

The temperature is a local temperature derived from the internal energy
of each eDPD particle based on the local equilibrium hypothesis.
For more details please see :ref:`(Espanol1997) <Espanol1997>` and
:ref:`(Li2014) <Li2014a>`.

**Output info:**

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input. See the
:doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS
output options.

The per-atom vector values will be in temperature :doc:`units <units>`.

Restrictions
""""""""""""


This compute is part of the USER-MESO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_style edpd <pair_meso>`

**Default:** none


----------


.. _Espanol1997:



**(Espanol1997)** Espanol, Europhys Lett, 40(6): 631-636 (1997).  DOI:
10.1209/epl/i1997-00515-8

.. _Li2014a:



**(Li2014)** Li, Tang, Lei, Caswell, Karniadakis, J Comput Phys, 265:
113-127 (2014).  DOI: 10.1016/j.jcp.2014.02.003.


