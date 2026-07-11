.. index:: compute ti

compute ti command
==================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group ti keyword args ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* ti = style name of this compute command
* one or more attribute/arg pairs may be appended
* keyword = pair style (lj/cut, gauss, born, etc.) or *tail* or *kspace*

  .. parsed-literal::

       pair style args = atype v_name1 v_name2
         atype = atom type (see asterisk form below)
         v_name1 = variable with name1 that is energy scale factor and function of lambda
         v_name2 = variable with name2 that is derivative of v_name1 with respect to lambda
       *tail* args = atype v_name1 v_name2
         atype = atom type (see asterisk form below)
         v_name1 = variable with name1 that is energy tail correction scale factor and function of lambda
         v_name2 = variable with name2 that is derivative of v_name1 with respect to lambda
       *kspace* args = atype v_name1 v_name2
         atype = atom type (see asterisk form below)
         v_name1 = variable with name1 that is K-Space scale factor and function of lambda
         v_name2 = variable with name2 that is derivative of v_name1 with respect to lambda

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all ti lj/cut 1 v_lj v_dlj coul/long 2 v_c v_dc kspace 1 v_ks v_dks
   compute 1 all ti lj/cut 1*3 v_lj v_dlj coul/long * v_c v_dc kspace * v_ks v_dks

Description
"""""""""""

Define a computation that calculates the derivative of the interaction
potential with respect to *lambda*, the coupling parameter used in a
thermodynamic integration.  This derivative can be used to infer a
free energy difference resulting from an alchemical simulation, as
described in :ref:`Eike <Eike>`.

Typically this compute will be used in conjunction with the
:doc:`fix adapt <fix_adapt>` command which can perform alchemical
transformations by adjusting the strength of an interaction potential
as a simulation runs, as defined by one or more
:doc:`pair_style <pair_style>` or :doc:`kspace_style <kspace_style>`
commands.  This scaling is done via a prefactor on the energy, forces,
virial calculated by the pair or :math:`k`-space style.  The prefactor is
often a function of a *lambda* parameter which may be adjusted from 0 to 1
(or vice versa) over the course of a :doc:`run <run>`.
The time-dependent adjustment is what the :doc:`fix adapt <fix_adapt>`
command does.

Assume that the unscaled energy of a pair_style or kspace_style is
given by :math:`U`.  Then the scaled energy is

.. math::

   U_s = f(\lambda) U

where :math:`f` is some function of :math:`\lambda`.  What this compute
calculates is

.. math::

   \frac{dU_s}{d\lambda} = U \frac{df(\lambda)}{d\lambda}
     = \frac{U_s}{f(\lambda)} \frac{df(\lambda)}{d\lambda},

which is the derivative of the system's scaled potential energy :math:`U_s`
with respect to :math:`\lambda`.

To perform this calculation, you provide one or more atom types as
*atype*\ .  The variable *atype* can be specified in one of two ways.
An explicit numeric value can be used, as in the first example above, or a
wildcard asterisk can be used in place of or in conjunction with the
*atype* argument to select multiple atom types.  This takes the form
"\*" or "\*n" or "m\*" or "m\*n".  If :math:`N` is the number of atom types,
then an asterisk with no numeric values means all types from 1 to :math:`N`.
A leading asterisk means all types from 1 to n (inclusive).  A trailing
asterisk means all types from m to N (inclusive).  A middle asterisk
means all types from m to n (inclusive).

You also specify two functions, as :doc:`equal-style variables <variable>`.
The first is specified as *v_name1*, where *name1* is the name of the
variable, and is :math:`f(\lambda)` in the notation above.  The second is
specified as *v_name2*, where *name2* is the name of the variable, and is
:math:`df(\lambda)/d\lambda` in the notation above (i.e., it is the analytic
derivative of :math:`f` with respect to :math:`\lambda`).
Note that the *name1* variable is also typically given as an
argument to the :doc:`fix adapt <fix_adapt>` command.

An alchemical simulation may use several pair potentials together,
invoked via the :doc:`pair_style hybrid or hybrid/overlay <pair_hybrid>`
command.  The total :math:`dU_s/d\lambda` for the overall system is calculated
as the sum of each contributing term as listed by the keywords in the
:doc:`compute ti <compute_ti>` command.  Individual pair potentials can be
listed, which will be sub-styles in the hybrid case.  You can also include a
:math:`k`-space term via the *kspace* keyword.  You can also include a pairwise
long-range tail correction to the energy via the *tail* keyword.

For each term, you can specify a different (or the same) scale factor
by the two variables that you list.  Again, these will typically
correspond toe the scale factors applied to these various potentials
and the :math:`k`-space contribution via the :doc:`fix adapt <fix_adapt>`
command.

More details about the exact functional forms for the computation of
:math:`du/dl` can be found in the paper by :ref:`Eike <Eike>`.

----------

Output info
"""""""""""

This compute calculates a global scalar, namely :math:`dU_s/d\lambda`.  This
value can be used by any command that uses a global scalar value from
a compute as input.  See the :doc:`Howto output <Howto_output>` doc page
for an overview of LAMMPS output options.

The scalar value calculated by this compute is "extensive".

The scalar value will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix adapt <fix_adapt>`

Default
"""""""

none

----------

.. _Eike:

**(Eike)** Eike and Maginn, Journal of Chemical Physics, 124, 164503 (2006).
