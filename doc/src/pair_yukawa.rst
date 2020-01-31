.. index:: pair\_style yukawa

pair\_style yukawa command
==========================

pair\_style yukawa/gpu command
==============================

pair\_style yukawa/omp command
==============================

pair\_style yukawa/kk command
=============================

Syntax
""""""


.. parsed-literal::

   pair_style yukawa kappa cutoff

* kappa = screening length (inverse distance units)
* cutoff = global cutoff for Yukawa interactions (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style yukawa 2.0 2.5
   pair_coeff 1 1 100.0 2.3
   pair_coeff \* \* 100.0

Description
"""""""""""

Style *yukawa* computes pairwise interactions with the formula

.. image:: Eqs/pair_yukawa.jpg
   :align: center

Rc is the cutoff.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* A (energy\*distance units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global yukawa
cutoff is used.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, the A coefficient and cutoff
distance for this pair style can be mixed.  A is an energy value mixed
like a LJ epsilon.  The default mix value is *geometric*\ .  See the
"pair\_modify" command for details.

This pair style supports the :doc:`pair_modify <pair_modify>` shift
option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none


