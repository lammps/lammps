.. index:: pair_style zbl

pair_style zbl command
======================

pair_style zbl/gpu command
==========================

pair_style zbl/kk command
=========================

pair_style zbl/omp command
==========================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style zbl inner outer

* inner = distance where switching function begins
* outer = global cutoff for ZBL interaction

Examples
""""""""


.. code-block:: LAMMPS

   pair_style zbl 3.0 4.0
   pair_coeff * * 73.0 73.0
   pair_coeff 1 1 14.0 14.0

Description
"""""""""""

Style *zbl* computes the Ziegler-Biersack-Littmark (ZBL) screened nuclear
repulsion for describing high-energy collisions between atoms.
:ref:`(Ziegler) <Ziegler>`. It includes an additional switching function
that ramps the energy, force, and curvature smoothly to zero
between an inner and outer cutoff. The potential
energy due to a pair of atoms at a distance r\_ij is given by:

.. math::

   E^{ZBL}_{ij} & = \frac{1}{4\pi\epsilon_0} \frac{Z_i Z_j \,e^2}{r_{ij}} \phi(r_{ij}/a)+ S(r_{ij})\\
   a & =  \frac{0.46850}{Z_{i}^{0.23} + Z_{j}^{0.23}}\\
   \phi(x) & =  0.18175e^{-3.19980x} + 0.50986e^{-0.94229x} + 0.28022e^{-0.40290x} + 0.02817e^{-0.20162x}\\

where *e* is the electron charge, :math:`\epsilon_0` is the electrical
permittivity of vacuum, and :math:`Z_i` and :math:`Z_j` are the nuclear
charges of the
two atoms.  The switching function :math:`S(r)` is identical to that used by
:doc:`pair_style lj/gromacs <pair_gromacs>`.  Here, the inner and outer
cutoff are the same for all pairs of atom types.

The following coefficients must be defined for each pair of atom types
via the :doc:`pair_coeff <pair_coeff>` command as in the examples above,
or in the LAMMPS data file.

* :math:`Z_i` (atomic number for first atom type, e.g. 13.0 for aluminum)

* :math:`Z_j` (ditto for second atom type)

The values of :math:`Z_i` and :math:`Z_j` are normally equal to the atomic
numbers of the two atom types. Thus, the user may optionally
specify only the coefficients for each :math:`i==i` pair, and rely
on the obvious mixing rule for cross interactions (see below).
Note that when :math:`i==i` it is required that :math:`Z_i == Z_j`.
When used with :doc:`hybrid/overlay <pair_hybrid>` and pairs are
assigned
to more than one sub-style, the mixing rule is not used and
each pair of types interacting with the ZBL sub-style must
be included in a pair\_coeff command.

.. note::

   The numerical values of the exponential decay constants in the
   screening function depend on the unit of distance. In the above
   equation they are given for units of angstroms. LAMMPS will
   automatically convert these values to the distance unit of the
   specified LAMMPS :doc:`units <units>` setting.  The values of Z should
   always be given as multiples of a proton's charge, e.g. 29.0 for
   copper.


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

For atom type pairs *i,j* and :math:`i \neq i`, the :math:`Z_i` and
:math:`Z_j` coefficients
can be mixed by taking :math:`Z_i` and :math:`Z_j` from the values
specified for
:math:`i == i` and :math:`j == j` cases. When used
with :doc:`hybrid/overlay <pair_hybrid>` and pairs are assigned
to more than one sub-style, the mixing rule is not used and
each pair of types interacting with the ZBL sub-style
must be included in a pair\_coeff command.
The :doc:`pair_modify <pair_modify>` mix option has no effect on
the mixing behavior

The ZBL pair style does not support the :doc:`pair_modify <pair_modify>`
shift option, since the ZBL interaction is already smoothed to 0.0 at
the cutoff.

The :doc:`pair_modify <pair_modify>` table option is not relevant for
this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure, since there are no corrections for a potential that goes to
0.0 at the cutoff.

This pair style does not write information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands must be
specified in an input script that reads a restart file.

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


.. _Ziegler:



**(Ziegler)** J.F. Ziegler, J. P. Biersack and U. Littmark, "The
Stopping and Range of Ions in Matter," Volume 1, Pergamon, 1985.
