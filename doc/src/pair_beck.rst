.. index:: pair\_style beck

pair\_style beck command
========================

pair\_style beck/gpu command
============================

pair\_style beck/omp command
============================

Syntax
""""""


.. parsed-literal::

   pair_style beck Rc

* Rc = cutoff for interactions (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style beck 8.0
   pair_coeff \* \* 399.671876712 0.0000867636112694 0.675 4.390 0.0003746
   pair_coeff 1 1 399.671876712 0.0000867636112694 0.675 4.390 0.0003746 6.0

Description
"""""""""""

Style *beck* computes interactions based on the potential by
:ref:`(Beck) <Beck>`, originally designed for simulation of Helium.  It
includes truncation at a cutoff distance Rc.

.. image:: Eqs/pair_beck.jpg
   :align: center

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands.

* A (energy units)
* B (energy-distance\^6 units)
* a (distance units)
* alpha (1/distance units)
* beta  (1/distance\^6 units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global cutoff
Rc is used.


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

For atom type pairs I,J and I != J, coefficients must be specified.
No default mixing rules are used.

This pair style does not support the :doc:`pair_modify <pair_modify>` shift
option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections.

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


----------


.. _Beck:



**(Beck)** Beck, Molecular Physics, 14, 311 (1968).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
