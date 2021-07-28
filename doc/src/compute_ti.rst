.. index:: compute ti

compute ti command
==================

Syntax
""""""

.. parsed-literal::

   compute ID group ti keyword args ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* ti = style name of this compute command
* one or more attribute/arg pairs may be appended
* keyword = pair style (lj/cut, gauss, born, etc) or *tail* or *kspace*

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

Typically this compute will be used in conjunction with the :doc:`fix adapt <fix_adapt>` command which can perform alchemical
transformations by adjusting the strength of an interaction potential
as a simulation runs, as defined by one or more
:doc:`pair_style <pair_style>` or :doc:`kspace_style <kspace_style>`
commands.  This scaling is done via a prefactor on the energy, forces,
virial calculated by the pair or K-Space style.  The prefactor is
often a function of a *lambda* parameter which may be adjusted from 0
to 1 (or vice versa) over the course of a :doc:`run <run>`.  The
time-dependent adjustment is what the :doc:`fix adapt <fix_adapt>`
command does.

Assume that the unscaled energy of a pair_style or kspace_style is
given by U.  Then the scaled energy is

.. parsed-literal::

   Us = f(lambda) U

where f() is some function of lambda.  What this compute calculates is

.. parsed-literal::

   dUs / d(lambda) = U df(lambda)/dlambda = Us / f(lambda) df(lambda)/dlambda

which is the derivative of the system's scaled potential energy Us
with respect to *lambda*\ .

To perform this calculation, you provide one or more atom types as
*atype*\ .  *Atype* can be specified in one of two ways.  An explicit
numeric values can be used, as in the first example above.  Or a
wildcard asterisk can be used in place of or in conjunction with the
*atype* argument to select multiple atom types.  This takes the form
"\*" or "\*n" or "n\*" or "m\*n".  If N = the number of atom types, then
an asterisk with no numeric values means all types from 1 to N.  A
leading asterisk means all types from 1 to n (inclusive).  A trailing
asterisk means all types from n to N (inclusive).  A middle asterisk
means all types from m to n (inclusive).

You also specify two functions, as :doc:`equal-style variables <variable>`.  The first is specified as *v_name1*, where
*name1* is the name of the variable, and is f(lambda) in the notation
above.  The second is specified as *v_name2*, where *name2* is the
name of the variable, and is df(lambda) / dlambda in the notation
above.  I.e. it is the analytic derivative of f() with respect to
lambda.  Note that the *name1* variable is also typically given as an
argument to the :doc:`fix adapt <fix_adapt>` command.

An alchemical simulation may use several pair potentials together,
invoked via the :doc:`pair_style hybrid or hybrid/overlay <pair_hybrid>`
command.  The total dUs/dlambda for the overall system is calculated
as the sum of each contributing term as listed by the keywords in the
compute ti command.  Individual pair potentials can be listed, which
will be sub-styles in the hybrid case.  You can also include a K-space
term via the *kspace* keyword.  You can also include a pairwise
long-range tail correction to the energy via the *tail* keyword.

For each term you can specify a different (or the same) scale factor
by the two variables that you list.  Again, these will typically
correspond toe the scale factors applied to these various potentials
and the K-Space contribution via the :doc:`fix adapt <fix_adapt>`
command.

More details about the exact functional forms for the computation of
du/dl can be found in the paper by :ref:`Eike <Eike>`.

----------

Output info
"""""""""""

This compute calculates a global scalar, namely dUs/dlambda.  This
value can be used by any command that uses a global scalar value from
a compute as input.  See the :doc:`Howto output <Howto_output>` doc page
for an overview of LAMMPS output options.

The scalar value calculated by this compute is "extensive".

The scalar value will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix adapt <fix_adapt>`

Default
"""""""

none

----------

.. _Eike:

**(Eike)** Eike and Maginn, Journal of Chemical Physics, 124, 164503 (2006).
