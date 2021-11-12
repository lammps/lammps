.. index:: pair_style smatb
.. index:: pair_style smatb/single

pair_style smatb command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *smatb*
* args = none

.. parsed-literal::

     *smatb*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style smatb 2.5
   pair_coeff 1 1 2.88 10.35	4.178	0.210	1.818	4.07293506	4.9883063257983666

Description
"""""""""""

The *lj/cut* styles compute the standard 12/6 Lennard-Jones potential,
given by

.. math::

      E_i =
       \sum_{j,R_{ij}\leq R_c} A  e^{-p \left(\frac{R_{ij}}{R_{0}}-1\right)}
       -\sqrt{\sum_{j,R_{ij}\leq R_c}\Xi^2 e^{-2q\left(\frac{R_{ij}}{R_{0}}-1\right)}}.

:math:`r_c` is the cutoff.

See the :doc:`lj/cut/coul <pair_lj_cut_coul>` styles to add a Coulombic
pairwise interaction and the :doc:`lj/cut/tip4p <pair_lj_cut_tip4p>` styles to
add the TIP4P water model.

Coefficients
""""""""""""

The following coefficients must be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data
file or restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, or by mixing as described below:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* LJ cutoff (distance units)

Note that :math:`\sigma` is defined in the LJ formula as the zero-crossing
distance for the potential, not as the energy minimum at :math:`2^{\frac{1}{6}} \sigma`.

The last coefficient is optional.  If not specified, the global
LJ cutoff specified in the pair_style command are used.

See the :doc:`run_style <run_style>` command for details.

----------

Related commands
""""""""""""""""

* :doc:`pair_coeff <pair_coeff>`

Default
"""""""

none
