.. index:: pair_style wf

pair_style wf command
===========================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style wf cutoff

* cutoff = cutoff for wf interactions (distance units)

Examples
""""""""


.. code-block:: LAMMPS

   variable         sigma  equal   1.0
   variable         epsilon equal 1.0
   variable         nu   equal 1.0
   variable         mu   equal 1.0
   variable         rc equal 2.0*${sigma}

   pair_style        wf ${rc}
   pair_coeff        1 1 ${epsilon} ${sigma} ${nu} ${mu}  ${rc}	

Description
"""""""""""

The *wf* style computes the potential in :ref:`Wang2020 <Wang2020>`, which is given by:

.. math::

  \phi(r)= \epsilon \alpha \left(\left[{\sigma\over r}\right]^{2\mu} -1 \right)\left(\left[{r_c\over r}\right]^{2\mu}-1\right)^{2\nu}

with

.. math::
  \alpha=2\nu\left(\frac{r_c}{\sigma}\right)^{2\mu}\left[\frac{1+2\nu}{2\nu\left[(r_c/\sigma)^{2\mu}-1\right]}\right]^{2\nu+1}

and

.. math::
  r_{min}=r_c\left[\frac{1+2\nu}{1+2\nu(r_c/\sigma)^{2\nu}}\right]^{1/{2\nu}}

:math:`r_c` is the cutoff. 

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* :math:`\nu` 
* :math:`\mu`
* :math:`r_c` (distance units)

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

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
can be  mixed. Careful is required with the cut-off radius. 
The default mix value is *geometric*\ .  See the "pair\_modify" command
for details.

The :doc:`pair_modify <pair_modify>` tail option is not relevant
for this pair style as it goes to zero at the cut-off radius.


This pair style writes its information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style does not support the use of the *inner*\ , *middle*\ , and *outer*
keywords of the :doc:`run_style respa <run_style>` command.


----------


Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none


----------

.. _Wang2020:

**(Wang2020)** X. Wang, S. Ramírez-Hinestrosa, J. Dobnikar, and D. Frenkel, Phys. Chem. Chem. Phys. 22, 10624 (2020).
