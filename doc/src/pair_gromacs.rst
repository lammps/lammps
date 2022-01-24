.. index:: pair_style lj/gromacs
.. index:: pair_style lj/gromacs/gpu
.. index:: pair_style lj/gromacs/kk
.. index:: pair_style lj/gromacs/omp
.. index:: pair_style lj/gromacs/coul/gromacs
.. index:: pair_style lj/gromacs/coul/gromacs/kk
.. index:: pair_style lj/gromacs/coul/gromacs/omp

pair_style lj/gromacs command
=============================

Accelerator Variants: *lj/gromacs/gpu*, *lj/gromacs/kk*, *lj/gromacs/omp*

pair_style lj/gromacs/coul/gromacs command
==========================================

Accelerator Variants: *lj/gromacs/coul/gromacs/kk*, *lj/gromacs/coul/gromacs/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/gromacs* or *lj/gromacs/coul/gromacs*
* args = list of arguments for a particular style

.. parsed-literal::

     *lj/gromacs* args = inner outer
       inner, outer = global switching cutoffs for Lennard Jones
     *lj/gromacs/coul/gromacs* args = inner outer (inner2) (outer2)
       inner, outer = global switching cutoffs for Lennard Jones (and Coulombic if only 2 args)
       inner2, outer2 = global switching cutoffs for Coulombic (optional)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/gromacs 9.0 12.0
   pair_coeff * * 100.0 2.0
   pair_coeff 2 2 100.0 2.0 8.0 10.0

   pair_style lj/gromacs/coul/gromacs 9.0 12.0
   pair_style lj/gromacs/coul/gromacs 8.0 10.0 7.0 9.0
   pair_coeff * * 100.0 2.0

Description
"""""""""""

The *lj/gromacs* styles compute shifted LJ and Coulombic interactions
with an additional switching function S(r) that ramps the energy and force
smoothly to zero between an inner and outer cutoff.  It is a commonly
used potential in the `GROMACS <http://www.gromacs.org>`_ MD code and for
the coarse-grained models of :ref:`(Marrink) <Marrink>`.

.. math::

   E_{LJ} = & 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
            \left(\frac{\sigma}{r}\right)^6 \right] + S_{LJ}(r)
                       \qquad r < r_c \\
   E_C  = & \frac{C q_i q_j}{\epsilon  r} + S_C(r) \qquad r < r_c \\
   S(r) = & C \qquad r < r_1 \\
   S(r) = & \frac{A}{3} (r - r_1)^3 + \frac{B}{4} (r - r_1)^4 + C \qquad  r_1 < r < r_c \\
   A = & (-3 E'(r_c) + (r_c - r_1) E''(r_c))/(r_c - r_1)^2 \\
   B = & (2 E'(r_c) - (r_c - r_1) E''(r_c))/(r_c - r_1)^3 \\
   C = & -E(r_c) + \frac{1}{2} (r_c - r_1) E'(r_c) - \frac{1}{12} (r_c - r_1)^2 E''(r_c)

:math:`r_1` is the inner cutoff; :math:`r_c` is the outer cutoff.  The
coefficients A, B, and C are computed by LAMMPS to perform the shifting
and smoothing.  The function S(r) is actually applied once to each term
of the LJ formula and once to the Coulombic formula, so there are 2 or 3
sets of A,B,C coefficients depending on which pair_style is used.  The
boundary conditions applied to the smoothing function are as follows:
:math:`S'(r_1) = S''(r_1) = 0, S(r_c) = -E(r_c), S'(r_c) = -E'(r_c)`,
and :math:`S''(r_c) = -E''(r_c)`, where E(r) is the corresponding term
in the LJ or Coulombic potential energy function.  Single and double
primes denote first and second derivatives with respect to r,
respectively.

The inner and outer cutoff for the LJ and Coulombic terms can be the
same or different depending on whether 2 or 4 arguments are used in
the pair_style command.  The inner LJ cutoff must be > 0, but the
inner Coulombic cutoff can be >= 0.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* inner (distance units)
* outer (distance units)

Note that sigma is defined in the LJ formula as the zero-crossing
distance for the potential, not as the energy minimum at :math:`2^{1/6} \sigma`.

The last 2 coefficients are optional inner and outer cutoffs for style
*lj/gromacs*\ .  If not specified, the global *inner* and *outer* values
are used.

The last 2 coefficients cannot be used with style
*lj/gromacs/coul/gromacs* because this force field does not allow
varying cutoffs for individual atom pairs; all pairs use the global
cutoff(s) specified in the pair_style command.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/cut pair styles can be mixed.
The default mix value is *geometric*\ .  See the "pair_modify" command
for details.

None of the GROMACS pair styles support the
:doc:`pair_modify <pair_modify>` shift option, since the Lennard-Jones
portion of the pair interaction is already smoothed to 0.0 at the
cutoff.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

None of the GROMACS pair styles support the
:doc:`pair_modify <pair_modify>` tail option for adding long-range tail
corrections to energy and pressure, since there are no corrections for
a potential that goes to 0.0 at the cutoff.

All of the GROMACS pair styles write their information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do
not need to be specified in an input script that reads a restart file.

All of the GROMACS pair styles can only be used via the *pair*
keyword of the :doc:`run_style respa <run_style>` command.  They do not
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

.. _Marrink:

**(Marrink)** Marrink, de Vries, Mark, J Phys Chem B, 108, 750-760 (2004).
