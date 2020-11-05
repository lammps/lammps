.. index:: fix esfriction

fix esfriction command
======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID esfriction Z1 Z2 A1 A2 Ne_val E_pka ES_cutoff

* ID, group-ID are documented in :doc:`fix <fix>` command
* esfriction = style name of this fix command
* Z1, Z2 = atomic numbers of the PKA source and target atoms
* A1, A1 = atomic mass of the PKA source and target atoms
* Ne_val = no. of valance shell electrons
* E_pka = energy of PKA in metal units
* ES_cutoff = minimum kinetic energy for electronic stopping (metal units)



Examples
""""""""


.. parsed-literal::

   fix 2 mobile esfriction 26 26 55.847 55.847 2 5000 60.0


Description
"""""""""""

This fix offers an alternative implementation to apply electronic stopping (ES). The other implementation in LAMMPS is :doc:`fix electron/stopping <fix_electron_stopping>`. Our implementation computes and applies the frictional force due to ES to each atom (with energy above cutoff) on-the-fly. Which means, it does not require any pre-computed table as input.It applies effect of ES by damping the velocity of an atom by a viscous force:

.. math::

  \vec{F_{es}} = \beta \vec{v}

where :math:`\beta` is the drag coefficient having units of mass/time, and :math:`\vec{v}` is the velocity of the energetic atom. Here, the value of :math:`\beta` is calculated using Lindhard and Scharff's (LS) formula for energy dissipation in energetic ions as described in :ref:`our paper <esfriction_paper>`. 
The fix will also output the total energy lost due to ES (computed numerically) to screen and logfile. For validation, it also writes the expected energy lost due to ES as computed by the NRT model. 


**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` options are not supported.

This fix does not set any global scalar or vector.

The *start/stop* keywords of the :doc:`run <run>` command have no effect
on this fix.


Restrictions
""""""""""""
This fix is part of the USER-MISC package. It is only enabled if LAMMPS was 
built with that package. See the :doc:`Build package <Build_package>` doc 
page for more info. This fix only supports metal units.

Related commands
""""""""""""""""

:doc:`fix drag <fix_drag>`, :doc:`fix electron/stopping <fix_electron_stopping>`

**Default:** none

----------

.. _esfriction_paper:

**(ESFriction)** H. Hemani and A. Majalee and U. Bhardwaj and A. Arya and K. Nordlund and M. Warrier, Inclusion and validation of electronic stopping in the open source LAMMPS code, https://arxiv.org/abs/2005.11940.

