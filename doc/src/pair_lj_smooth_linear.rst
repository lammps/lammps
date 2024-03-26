.. index:: pair_style lj/smooth/linear
.. index:: pair_style lj/smooth/linear/omp

pair_style lj/smooth/linear command
===================================

Accelerator Variants: *lj/smooth/linear/omp*

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
continuously to zero at the cutoff :math:`r_c` :ref:`(Toxvaerd) <Toxvaerd>`:

.. math::

   \phi\left(r\right) & =  4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
                       \left(\frac{\sigma}{r}\right)^6 \right] \\
   E\left(r\right) & =  \phi\left(r\right)  - \phi\left(r_c\right) - \left(r - r_c\right) \left.\frac{d\phi}{d r} \right|_{r=r_c}       \qquad r < r_c

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* cutoff (distance units)

The last coefficient is optional. If not specified, the global
LJ cutoff specified in the pair_style command is used.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance can be mixed. The default mix value is geometric.
See the "pair_modify" command for details.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift option for the energy of the pair interaction, since it goes
to 0.0 at the cutoff by construction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure, since the energy of the pair interaction is smoothed to 0.0
at the cutoff.

This pair style writes its information to :doc:`binary restart files <restart>`,
so pair_style and pair_coeff commands do not need to be specified
in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the EXTRA-PAIR package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair lj/smooth <pair_lj_smooth>`

Default
"""""""

none

----------

.. _Toxvaerd:

**(Toxvaerd)** Toxvaerd, Dyre, J Chem Phys, 134, 081102 (2011).
