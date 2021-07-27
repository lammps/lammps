.. index:: pair_style lj/cubic
.. index:: pair_style lj/cubic/gpu
.. index:: pair_style lj/cubic/omp

pair_style lj/cubic command
===========================

Accelerator Variants: *lj/cubic/gpu*, *lj/cubic/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/cubic

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cubic
   pair_coeff * * 1.0 0.8908987

Description
"""""""""""

The *lj/cubic* style computes a truncated LJ interaction potential
whose energy and force are continuous everywhere.  Inside the
inflection point the interaction is identical to the standard 12/6
:doc:`Lennard-Jones <pair_lj>` potential.  The LJ function outside the
inflection point is replaced with a cubic function of distance. The
energy, force, and second derivative are continuous at the inflection
point.  The cubic coefficient A3 is chosen so that both energy and
force go to zero at the cutoff distance.  Outside the cutoff distance
the energy and force are zero.

.. math::

   E & = u_{LJ}(r) \qquad r \leq r_s \\
     & = u_{LJ}(r_s) + (r-r_s) u'_{LJ}(r_s) - \frac{1}{6} A_3 (r-r_s)^3 \qquad r_s < r \leq r_c \\
     & = 0 \qquad r > r_c

The location of the inflection point :math:`r_s` is defined
by the LJ diameter, :math:`r_s/\sigma = (26/7)^{1/6}`. The cutoff distance
is defined by :math:`r_c/r_s = 67/48` or :math:`r_c/\sigma = 1.737...`
The analytic expression for the
the cubic coefficient
:math:`A_3 r_{min}^3/\epsilon = 27.93...` is given in the paper by
Holian and Ravelo :ref:`(Holian) <Holian>`.

This potential is commonly used to study the shock mechanics of FCC
solids, as in Ravelo et al. :ref:`(Ravelo) <Ravelo2>`.

The following coefficients must be defined for each pair of atom types
via the :doc:`pair_coeff <pair_coeff>` command as in the example above,
or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)

Note that :math:`\sigma` is defined in the LJ formula as the
zero-crossing distance for the potential, not as the energy minimum,
which is located at :math:`r_{min} = 2^{\frac{1}{6}} \sigma`. In the
above example, :math:`\sigma = 0.8908987`, so :math:`r_{min} = 1.0`.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/cut pair styles can be mixed.
The default mix value is *geometric*\ .  See the "pair_modify" command
for details.

The lj/cubic pair style does not support the
:doc:`pair_modify <pair_modify>` shift option,
since pair interaction is already smoothed to 0.0 at the
cutoff.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

The lj/cubic pair style does not support the
:doc:`pair_modify <pair_modify>` tail option for adding long-range tail
corrections to energy and pressure, since there are no corrections for
a potential that goes to 0.0 at the cutoff.

The lj/cubic pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do
not need to be specified in an input script that reads a restart file.

The lj/cubic pair style can only be used via the *pair*
keyword of the :doc:`run_style respa <run_style>` command.  It does not
support the *inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the EXTRA-PAIR package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Holian:

.. _Ravelo2:

**(Holian)** Holian and Ravelo, Phys Rev B, 51, 11275 (1995).

**(Ravelo)** Ravelo, Holian, Germann and Lomdahl, Phys Rev B, 70, 014103 (2004).
