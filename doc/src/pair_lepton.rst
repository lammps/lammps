.. index:: pair_style lepton
.. index:: pair_style lepton/omp

pair_style lepton command
=========================

Accelerator Variants: *lepton/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lepton*
* args = list of arguments for a particular style

.. parsed-literal::

    *lepton* args = cutoff
      cutoff = global cutoff for the interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lepton 2.5
   pair_coeff * * "k*((r-r0)^2*step(r0-r)); k=200; r0=1.5" 2.0
   pair_coeff 1 2 "4.0*eps*((sig/r)^12 - (sig/r)^6);eps=1.0;sig=1.0" 1.12246204830937
   pair_coeff 2 2 "eps*(2.0*(sig/r)^9 - 3.0*(sig/r)^6);eps=1.0;sig=1.0"

Description
"""""""""""

.. versionadded:: TBD

Pair style *lepton* computes spherical pairwise interactions based on
evaluating strings.  The potential function must be provided as an
expression string using "r" as the distance variable, for example
`"200.0*(r-1.5)^2"` represents a harmonic potential with equilibrium
distance :math:`r_0` of 1.5 distance units and a force constant *K* of
200.0 energy units:

.. math::

   U_{ij} = K (r-r_0)^2

The `Lepton library <https://simtk.org/projects/lepton>`_, that the
*lepton* pair style interfaces with, evaluates this expression string at
run time to compute the pairwise energy.  It also creates an
analytical representation of the differentiation of this expression with
respect to "r" and then uses that to compute the force between the pairs
of particles within the given cutoff.

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

----------

.. include:: lepton_expression.rst

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Pair style *lepton* does **not** support mixing.  Thus, expressions for
all I,J pairs must be specified explicitly.

This pair style supports the :doc:`pair_modify <pair_modify>`
shift option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table options are not relevant for
the this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files
<restart>`, so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

The *lepton* pair style is part of the LEPTON package and only enabled if
LAMMPS was built with this package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style python <pair_python>`,
:doc:`pair_style table <pair_table>`, :doc:`pair_write <pair_write>`

Default
"""""""

none
