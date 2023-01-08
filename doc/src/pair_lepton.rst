.. index:: pair_style lepton
.. index:: pair_style lepton/omp
.. index:: pair_style lepton/coul
.. index:: pair_style lepton/coul/omp

pair_style lepton command
=========================

Accelerator Variants: *lepton/omp*, *lepton/coul/comp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lepton* or *lepton/coul*
* args = list of arguments for a particular style

.. parsed-literal::

    *lepton* args = cutoff
      cutoff = global cutoff for the interactions (distance units)
    *lepton/coul* args = cutoff keyword
      cutoff = global cutoff for the interactions (distance units)
      zero or more keywords may be appended
      keyword = *ewald* or *pppm* or *msm* or *dispersion* or *tip4p*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lepton 2.5

   pair_coeff  * *  "k*((r-r0)^2*step(r0-r)); k=200; r0=1.5" 2.0
   pair_coeff  1 2  "4.0*eps*((sig/r)^12 - (sig/r)^6);eps=1.0;sig=1.0" 1.12246204830937
   pair_coeff  2 2  "eps*(2.0*(sig/r)^9 - 3.0*(sig/r)^6);eps=1.0;sig=1.0"
   pair_coeff  1 3  "zbl(13,6,r)"
   pair_coeff  3 3  "(1.0-switch)*zbl(6,6,r)-switch*4.0*eps*((sig/r)^6);switch=0.5*(tanh(10.0*(r-sig))+1.0);eps=0.05;sig=3.20723"

   pair_style lepton/coul 2.5
   pair_coeff 1 1 "qi*qj/r" 4.0
   pair_coeff 1 2 "lj+coul; lj=4.0*eps*((sig/r)^12 - (sig/r)^6); eps=1.0; sig=1.0; coul=qi*qj/r"

   pair_style lepton/coul 2.5 pppm
   kspace_style pppm 1.0e-4
   pair_coeff 1 1 "qi*qj/r*erfc(alpha*r); alpha=1.067"

Description
"""""""""""

.. versionadded:: TBD

Pair styles *lepton* and *lepton/coul* compute pairwise interactions
between particles which depend solely on the distance and have a cutoff.
The potential function must be provided as an expression string using
"r" as the distance variable.  With pair style *lepton/coul* one may
additionally reference the charges of the two atoms of the pair with
"qi" and "qj", respectively.  Note that further constants in the
expression can be defined in the same string as additional expressions
separated by semi-colons as shown in the examples above.

The expression `"200.0*(r-1.5)^2"` represents a harmonic potential
around the pairwise distance :math:`r_0` of 1.5 distance units and a
force constant *K* of 200.0 energy units:

.. math::

   U_{ij} = K (r-r_0)^2

The expression `"qi*qj/r"` represents a regular Coulombic potential with cutoff:

.. math::

   U_{ij} = \frac{C q_i q_j}{\epsilon  r} \qquad r < r_c

The `Lepton library <https://simtk.org/projects/lepton>`_, that the
*lepton* pair style interfaces with, evaluates this expression string at
run time to compute the pairwise energy.  It also creates an analytical
representation of the first derivative of this expression with respect
to "r" and then uses that to compute the force between the pairs of
particles within the given cutoff.

The following coefficients must be defined for each pair of atoms types
via the :doc:`pair_coeff <pair_coeff>` command as in the examples above,
or in the data file or restart files read by the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands:

* Lepton expression (energy units)
* cutoff (distance units)

The Lepton expression must be either enclosed in quotes or must not
contain any whitespace so that LAMMPS recognizes it as a single keyword.
More on valid Lepton expressions below.  The last coefficient is
optional; it allows to set the cutoff for a pair of atom types to a
different value than the global cutoff.

For pair style *lepton* only the "lj" value of the :doc:`special_bonds <special_bonds>`
settings apply in case the interacting pair is also connected with a bond.
The potential energy will *only* be added to the "evdwl" property.

For pair style *lepton/coul* only the "coul" value of the :doc:`special_bonds <special_bonds>`
settings apply in case the interacting pair is also connected with a bond.
The potential energy will *only* be added to the "ecoul" property.

In addition to the functions listed below, both pair styles support in
addition a custom "zbl(zi,zj,r)" function which computes the
Ziegler-Biersack-Littmark (ZBL) screened nuclear repulsion for
describing high-energy collisions between atoms.  For details of the
function please see the documentation for :doc:`pair style zbl
<pair_zbl>`. The arguments of the function are the atomic numbers of
atom i (zi), atom j (zj) and the distance r.  Please see the examples
above.

----------

.. include:: lepton_expression.rst

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Pair styles *lepton* and *lepton/coul* do not support mixing.  Thus,
expressions for *all* I,J pairs must be specified explicitly.

Only pair style *lepton* supports the :doc:`pair_modify shift <pair_modify>`
option for shifting the energy of the pair interaction so that it is
0 at the cutoff, pair style *lepton/coul* does *not*.

The :doc:`pair_modify table <pair_modify>` options are not relevant for
the these pair styles.

These pair styles do not support the :doc:`pair_modify tail
<pair_modify>` option for adding long-range tail corrections to energy
and pressure.

These pair styles write its information to :doc:`binary restart files
<restart>`, so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

These pair styles are part of the LEPTON package and only enabled if
LAMMPS was built with this package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style python <pair_python>`,
:doc:`pair_style table <pair_table>`, :doc:`pair_write <pair_write>`

Default
"""""""

none
