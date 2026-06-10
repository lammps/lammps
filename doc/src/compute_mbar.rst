.. index:: compute mbar

compute mbar command
====================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID mbar temp attribute args ... keyword value ...

* ID, group-ID are documented in the :doc:`compute <compute>` command
* mbar = name of this compute command
* temp = external temperature (as specified for constant-temperature run)
* one or more attributes with args may be appended
* attribute = *pair* or *atom*

  .. parsed-literal::

       *pair* args = pstyle pparam I J v_name
         pstyle = pair style name (e.g., *lj/cut/soft*)
         pparam = parameter to set
         I,J = type pair(s) to set parameter for (I<=J)
         v_name = :doc:`vector-style variable <variable>` with the values of *pparam* at each state
       *atom* args = aparam I v_name
         aparam = *charge* = parameter to perturb
         I = type to set parameter for
         v_name = :doc:`vector-style variable <variable>` with the charge values at each state

* zero or more keyword/value pairs may be appended
* keyword = *tail*

  .. parsed-literal::

       *tail* value = *no* or *yes*
         *no* = ignore tail correction to pair energies
         *yes* = include tail correction to pair energies

Examples
""""""""

.. code-block:: LAMMPS

   variable lam vector [0.0,0.2,0.4,0.6,0.8,1.0]
   compute 1 all mbar 300 pair lj/cut/soft lambda 1 * v_lam

   variable lam vector [0.0,0.25,0.5,0.75,1.0]
   variable eps vector [0.0,0.05,0.1,0.15,0.2]
   compute 2 all mbar 298 pair lj/cut/soft epsilon 1 1 v_eps tail yes

   variable q vector [-0.8,-0.4,0.0]
   compute 3 all mbar 300 atom charge 1 v_q

Description
"""""""""""

This compute is the multistate analogue of :doc:`compute fep <compute_fep>`.
It applies perturbations to parameters of the interaction potential and
recalculates the potential energy *without* changing the atomic
coordinates from those of the sampled (reference) configuration. Whereas
:doc:`compute fep <compute_fep>` evaluates the potential energy at a single
perturbed state, this compute evaluates it at a series
of states defined by a vector of coupling-parameter (:math:`\lambda`)
values. This is the information needed by the multistate Bennett
acceptance ratio (MBAR) method :ref:`(Shirts) <Shirts>` to estimate free
energy differences along an alchemical transformation.

The MBAR estimator is built from the matrix of *reduced* potential
energies

.. math::

   u_{kn} = \frac{U_k(\mathbf{x}_n)}{k_B T}

where :math:`U_k(\mathbf{x}_n)` is the potential energy of configuration
:math:`n` (sampled from any of the states) evaluated with the
interaction parameters of state :math:`k`. For each invocation, this
compute produces one row of this matrix: the reduced potential energies
of the current configuration evaluated at every state requested.
Accumulating these rows over a trajectory (and over trajectories sampled at
the different states) yields the full :math:`u_{kn}` matrix, which is then
passed to an external MBAR solver such as
`pymbar <https://github.com/choderalab/pymbar>`_ to obtain the free energy
differences and their uncertainties.

Each perturbed parameter carries its *own* grid, supplied as a
:doc:`vector-style variable <variable>` referenced as *v_name* and listed
immediately after the type specification of that perturbation. The grid
holds the **absolute values** that the parameter takes at each of the
states, one entry per state. All perturbation grids must have the same
length, which defines the number of states; this length sets the length of
the global vector produced by the compute. The variables are evaluated once
when the compute is defined and their values are copied, so they must
already be defined and should hold fixed numeric values (e.g.
``variable lam vector [0.0,0.25,0.5,0.75,1.0]``) rather than quantities that
could change during the run.

At each sampling state of the MBAR method, the perturbed parameter is *set*
to the corresponding grid value (its original value is saved and restored
after the compute).

For the *pair* attribute the user decides which parameter to couple: it can
be the explicit activation (:math:`\lambda`) parameter of a soft-core pair
style such as :doc:`pair_style .../soft <pair_fep_soft>`, or any other pair
parameter such as *epsilon* or *sigma* of a Lennard-Jones potential. The
grid then lists the values that parameter takes along the alchemical path;
any non-linear schedule is encoded directly in those values. The grid for a
soft-core activation parameter typically runs from ``0.0`` (interaction off)
to ``1.0`` (fully coupled), but the sequence of values is to be adapted for
each specific parameter.

Electrostatic charges are coupled in the same way using the *atom charge*
attribute. Although the electrostatic charges can also be coupled through a 
soft-core Coulomb pair style, in many applications the van der Waals
interactions are coupled first and the charges are activated afterwards, so
it is convenient to be able to set the charges directly. The grid then
lists the charge of the selected atom type at each state (e.g. from ``0.0``
to the fully-coupled charge), and these values are assigned to all atoms of
that type in the group; the original charges are restored after the compute.

As with :doc:`compute fep <compute_fep>`, the kinetic energy is not taken
into account, so the masses of the particles should not be modified
between states, and internal coordinates (bond lengths, angles, etc.)
are not changed.

----------

The mechanics of the *pair* attribute (the *pstyle* matching, the *pparam*
extraction, and the I,J type specification including the wild-card type
syntax) work as described for :doc:`compute fep <compute_fep>`; see that
page for the supported pair styles and the type syntax. For this compute
*pparam* is *set* to the grid value at each state, as discussed above, so it
is typically the activation (:math:`\lambda`) parameter of a soft-core pair
style, but any extractable pair parameter may be chosen.

----------

The *tail* keyword controls the calculation of the tail correction to
"van der Waals" pair energies beyond the cutoff, if this has been
activated via the :doc:`pair_modify <pair_modify>` command.

----------

Output info
"""""""""""

This compute calculates a global vector of length equal to the number of
states (the common length of the perturbation grids). Element :math:`k`
(1-based as c_ID[k]) is the reduced potential energy
:math:`U_k(\mathbf{x})/k_B T` of the current configuration evaluated with
the parameters of state :math:`k`. The energies include kspace terms when
charges are perturbed; for a perturbation of pair parameters only, the
kspace contribution is the same in every state and cancels in the MBAR free
energy differences, so it is omitted.

These output results can be used by any command that uses a global vector
from a compute as input. See the :doc:`Howto output <Howto_output>` page
for an overview of LAMMPS output options. In a typical workflow the vector
is written to a file with :doc:`fix ave/time <fix_ave_time>` (using the
*mode vector* option) for later post-processing with an MBAR solver.

The values calculated by this compute are "extensive".

Restrictions
""""""""""""

This compute is distributed as the FEP package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

Related commands
""""""""""""""""

:doc:`compute fep <compute_fep>`, :doc:`fix adapt/fep <fix_adapt_fep>`,
:doc:`fix ave/time <fix_ave_time>`, :doc:`pair_style .../soft <pair_fep_soft>`

Default
"""""""

The option default is *tail* = *no*\ .

----------

.. _Shirts:

**(Shirts)** Shirts and Chodera, J Chem Phys, 129, 124105 (2008)
