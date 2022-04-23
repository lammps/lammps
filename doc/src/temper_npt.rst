.. index:: temper/npt

temper/npt command
==================

Syntax
""""""

.. parsed-literal::

   temper/npt  N M temp fix-ID seed1 seed2 pressure index

* N = total # of timesteps to run
* M = attempt a tempering swap every this many steps
* temp = initial temperature for this ensemble
* fix-ID = ID of the fix that will control temperature and pressure during the run
* seed1 = random # seed used to decide on adjacent temperature to partner with
* seed2 = random # seed for Boltzmann factor in Metropolis swap
* pressure = setpoint pressure for the ensemble
* index = which temperature (0 to N-1) I am simulating (optional)

Examples
""""""""

.. code-block:: LAMMPS

   temper/npt 100000 100 $t nptfix 0 58728 1
   temper/npt 2500000 1000 300 nptfix  0 32285 $p
   temper/npt 5000000 2000 $t nptfix 0 12523 1 $w

Description
"""""""""""

Run a parallel tempering or replica exchange simulation using multiple
replicas (ensembles) of a system in the isothermal-isobaric (NPT)
ensemble.  The command temper/npt works like :doc:`temper <temper>` but
requires running replicas in the NPT ensemble instead of the canonical
(NVT) ensemble and allows for pressure to be set in the ensembles.
These multiple ensembles can run in parallel at different temperatures
or different pressures.  The acceptance criteria for temper/npt is
specific to the NPT ensemble and can be found in references
:ref:`(Okabe) <Okabe2>` and :ref:`(Mori) <Mori2>`.

Apart from the difference in acceptance criteria and the specification
of pressure, this command works much like the :doc:`temper <temper>`
command. See the documentation on :doc:`temper <temper>` for information
on how the parallel tempering is handled in general.

----------

Restrictions
""""""""""""

This command can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` page for more
info.

This command should be used with a fix that maintains the
isothermal-isobaric (NPT) ensemble.

Related commands
""""""""""""""""

:doc:`temper <temper>`, :doc:`variable <variable>`, :doc:`fix_npt <fix_nh>`

Default
"""""""

none

.. _Okabe2:

**(Okabe)** T. Okabe, M. Kawata, Y. Okamoto, M. Masuhiro, Chem. Phys. Lett., 335, 435-439 (2001).

.. _Mori2:

**(Mori)** Y. Mori, Y. Okamoto, J. Phys. Soc. Jpn., 7, 074003 (2010).
