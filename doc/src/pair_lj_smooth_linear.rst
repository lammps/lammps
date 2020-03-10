.. index:: pair_style lj/smooth/linear

pair_style lj/smooth/linear command
===================================

pair_style lj/smooth/linear/omp command
=======================================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style lj/smooth/linear cutoff

* cutoff = global cutoff for Lennard-Jones interactions (distance units)

Examples
""""""""


.. code-block:: LAMMPS

   pair_style lj/smooth/linear 2.5
   pair_coeff * * 1.0 1.0
   pair_coeff 1 1 0.3 3.0 9.0

Description
"""""""""""

Style *lj/smooth/linear* computes a truncated and force-shifted LJ
interaction (aka Shifted Force Lennard-Jones) that combines the
standard 12/6 Lennard-Jones function and subtracts a linear term based
on the cutoff distance, so that both, the potential and the force, go
continuously to zero at the cutoff Rc :ref:`(Toxvaerd) <Toxvaerd>`:

.. math::

   \phi\left(r\right) & =  4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
                       \left(\frac{\sigma}{r}\right)^6 \right] \\
   E\left(r\right) & =  \phi\left(r\right)  - \phi\left(R_c\right) - \left(r - R_c\right) \left.\frac{d\phi}{d r} \right|_{r=R_c}       \qquad r < R_c


The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* cutoff (distance units)

The last coefficient is optional. If not specified, the global
LJ cutoff specified in the pair\_style command is used.


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
and cutoff distance can be mixed. The default mix value is geometric.
See the "pair\_modify" command for details.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift option for the energy of the pair interaction, since it goes
to 0.0 at the cutoff by construction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure, since the energy of the pair interaction is smoothed to 0.0
at the cutoff.

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

:doc:`pair_coeff <pair_coeff>`, :doc:`pair lj/smooth <pair_lj_smooth>`

**Default:** none


----------


.. _Toxvaerd:



**(Toxvaerd)** Toxvaerd, Dyre, J Chem Phys, 134, 081102 (2011).
