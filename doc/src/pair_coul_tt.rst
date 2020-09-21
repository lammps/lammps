.. index:: pair_style coul/tt

pair_style coul/tt command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *coul/tt* 
* args = list of arguments for a particular style

.. parsed-literal::

     *coul/tt* args = n cutoff
       n = degree of polynomial
       cutoff = global cutoff (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style hybrid/overlay ... coul/tt 4 12.0
   pair_coeff 1 2  coul/tt 4.5 1.0
   pair_coeff 1 2  coul/tt 4.0 1.0 4 12.0
   pair_coeff 1 3* coul/tt 4.5 1.0 4


Description
"""""""""""

The *coul/tt* pair styles are meant to be used with force fields that
include explicit polarization through Drude dipoles.

The *coul/tt* pair style should be used as a sub-style within in the
:doc:`pair_style hybrid/overlay <pair_hybrid>` command, in conjunction with a
main pair style including Coulomb interactions and *thole* pair style, 
or with *lj/cut/thole/long* pair style that is equivalent to the combination 
of two preceding.

The *coul/tt* pair styles compute the charge-dipole Coulomb interaction damped 
at short distances by a function

.. math::

  TT_{ij}(r) = 1 - c_{ij} \cdot e^{-b_{ij} r} \sum_{k=0}^n \frac{(b_{ij} r)^k}{k!}

This function results from an adaptation to Coulomb interaction :ref:`(Salanne) 
<Salanne1>` the damping function originally proposed
by :ref:`(Tang Toennies) <TangToennies1>` for van der Waals interactions.

The polynomial takes the degree of 4 for damping Coulomb interaction.
The parameters :math:`b_{ij}` and :math:`c_{ij}` could be determined from 
first-principle calculations for small, mainly mono-atomic, ions :ref:`(Salanne) 
<Salanne1>` or chosen as empirical for large molecules.

The damping function is typically applied to the interactions between a Drude 
charge (:math:`q_{D,i}` on a Drude particle or :math:`-q_{D,i}` on the respective 
Drude core particle bonded to a Drude particle) and a charge of a non-polarizable 
atom, :math:`q_{j}`. The Tang-Toennies function could be used to damp electrostatic 
interaction between two Drude cores acting on the partial charge of the one core 
:math:`q_{i}-q_{D,i}` and a Drude charge of the another one :math:`-q_{D,j}`, and 
the opposite case, respectively. The :math:`b_{ij}` and :math:`c_{ij}` are equal 
to :math:`b_{ji}` and :math:`c_{ji}` in case of core - core interaction.
Therefore, the screening is not applied to the full charge of the Drude core 
:math:`q_i`, but only to the :math:`-q_{D,i}` part of it when it acts as a 
dipole or :math:`q_{i}-q_{D,i}` when it acts as a charge. 

For pair_style *coul/tt*\ , the following coefficients must be defined for
each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command
as in the example above.

* :math:`b_{ij}`
* :math:`c_{ij}`
* degree of polynomial (positive integer)
* cutoff (distance units)

The last two coefficients are optional.  If not specified the global 
degree of polynomial or global cutoff specified in the pair_style
command are used. In order to specify a cutoff (forth argument) a damp
parameter (third argument) must also be specified.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The *coul/tt* pair style does not support mixing.  Thus, coefficients
for all I,J pairs must be specified explicitly.

Restrictions
""""""""""""

These pair styles are part of the USER-DRUDE package. They are only
enabled if LAMMPS was built with that package. See the :doc:`Build package 
<Build_package>` doc page for more info.

This pair_style should currently not be used with the :doc:`charmm dihedral 
style <dihedral_charmm>` if the latter has non-zero 1-4 weighting
factors. This is because the *coul/tt* pair style does not know which
pairs are 1-4 partners of which dihedrals.

Related commands
""""""""""""""""

:doc:`fix drude <fix_drude>`, :doc:`fix langevin/drude <fix_langevin_drude>`, 
:doc:`fix drude/transform <fix_drude_transform>`, 
:doc:`compute temp/drude <compute_temp_drude>`,
:doc:`pair_style thole <pair_thole>`

Default
"""""""

none

----------

.. _Thole1:

**(Thole)** Chem Phys, 59, 341 (1981).

.. _Salanne1:

**(Salanne)** Salanne, Rotenberg, Jahn, Vuilleumier, Simon, Christian and Madden, Theor Chem Acc, 131, 1143 (2012).

.. _TangToennies1:

**(Tang and Toennies)** J Chem Phys, 80, 3726 (1984).
