.. index:: pair_style zero

pair_style zero command
=======================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style zero cutoff [nocoeff]

* zero = style name of this pair style
* cutoff = global cutoff (distance units)
* nocoeff = ignore all pair\_coeff parameters (optional)

Examples
""""""""


.. code-block:: LAMMPS

   pair_style zero 10.0
   pair_style zero 5.0 nocoeff
   pair_coeff * *
   pair_coeff 1 2*4 3.0

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
used to insure communication of ghost atoms even when a pair style is
not defined, but it will not trigger neighbor list generation.

The optional *nocoeff* flag allows to read data files with a PairCoeff
section for any pair style. Similarly, any pair\_coeff commands
will only be checked for the atom type numbers and the rest ignored.
In this case, only the global cutoff will be used.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* cutoff (distance units)

This coefficient is optional.  If not specified, the global cutoff
specified in the pair\_style command is used. If the pair\_style has
been specified with the optional *nocoeff* flag, then a cutoff
pair coefficient is ignored.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

The cutoff distance for this pair style can be mixed.  The default mix
value is *geometric*\ .  See the "pair\_modify" command for details.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style supports the use of the *inner*\ , *middle*\ ,
and *outer* keywords of the :doc:`run_style respa <run_style>` command.


----------


Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`pair_style none <pair_none>`

**Default:** none
