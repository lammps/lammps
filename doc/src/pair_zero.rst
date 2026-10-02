.. index:: pair_style zero
.. index:: pair_style zero/coul

pair_style zero command
=======================

pair_style zero/coul command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style cutoff [nocoeff] [full]

* style = *zero* or *zero/coul*
* cutoff = global cutoff (distance units)
* nocoeff = ignore all pair_coeff parameters (optional)
* full = build full neighbor list (optional)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style zero 10.0
   pair_style zero 5.0 nocoeff
   pair_coeff * *
   pair_coeff 1 2*4 3.0
   pair_style zero/coul 12.0

Description
"""""""""""

Define a global or per-type cutoff length for the purpose of
building a neighbor list and acquiring ghost atoms, but do
not compute any pairwise forces or energies.

This can be useful for fixes or computes which require a neighbor list
to enumerate pairs of atoms within some cutoff distance, but when
pairwise forces are not otherwise needed.  Examples are the :doc:`fix bond/create <fix_bond_create>`, :doc:`compute rdf <compute_rdf>`,
:doc:`compute voronoi/atom <compute_voronoi_atom>` commands.

Note that the :doc:`comm_modify cutoff <comm_modify>` command can be
used to ensure communication of ghost atoms even when a pair style is
not defined, but it will not trigger neighbor list generation.

The optional *nocoeff* flag allows to read data files with a PairCoeff
section for any pair style. Similarly, any pair_coeff commands
will only be checked for the atom type numbers and the rest ignored.
In this case, only the global cutoff will be used.

.. versionadded:: 3Nov2022

The optional *full* flag builds a full neighbor list instead of the default
half neighbor list.

.. versionadded:: 2Sep2026

Pair style *zero/coul* behaves exactly like pair style *zero*, but
additionally presents itself as a Coulombic pair style: it declares
compatibility with :doc:`kspace styles <kspace_style>` and provides its
cutoff as the real-space Coulomb cutoff, while still computing no
pairwise interactions.  This allows computing only the k-space
contribution to forces and energies, or satisfying commands that require
the presence of a Coulombic pair style, for example for testing and
debugging purposes.  Note that the k-space calculation then misses the
compensating real-space contribution, so the resulting forces and
energies do not correspond to a complete Coulomb interaction.  Since the
compatibility checks between kspace and pair styles are applied in both
directions, pair style *zero/coul* cannot be used with kspace styles
that require a pair style providing dispersion or TIP4P support (e.g.
*pppm/disp* or *pppm/tip4p*).

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* cutoff (distance units)

This coefficient is optional.  If not specified, the global cutoff
specified in the pair_style command is used. If the pair_style has
been specified with the optional *nocoeff* flag, then a cutoff
pair coefficient is ignored.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The cutoff distance for this pair style can be mixed.  The default mix
value is *geometric*\ .  See the "pair_modify" command for details.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style supports the use of the *inner*, *middle*,
and *outer* keywords of the :doc:`run_style respa <run_style>` command.

----------

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`pair_style none <pair_none>`

Default
"""""""

none
