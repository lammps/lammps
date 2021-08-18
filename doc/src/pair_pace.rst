.. index:: pair_style pace

pair_style pace command
========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style pace ... keyword values ...

* an optional keyword may be appended
* keyword = *product* or *recursive*

  .. parsed-literal::

       *product* = use product algorithm for basis functions
       *recursive* = use recursive algorithm for basis functions

Examples
""""""""

.. code-block:: LAMMPS

   pair_style pace
   pair_style pace product
   pair_coeff * * Cu-PBE-core-rep.ace Cu

Description
"""""""""""

Pair style *pace* computes interactions using the Atomic Cluster
Expansion (ACE), which is a general expansion of the atomic energy in
multi-body basis functions. :ref:`(Drautz) <Drautz20191>`.
The *pace* pair style
provides an efficient implementation that
is described in this paper :ref:`(Lysogorskiy) <Lysogorskiy20211>`.

In ACE, the total energy is decomposed into a sum over
atomic energies. The energy of atom *i* is expressed as a
linear or non-linear function of one or more density functions.
By projecting the
density onto a local atomic base, the lowest order contributions
to the energy can be expressed as a set of scalar polynomials in
basis function contributions summed over neighbor atoms.

Only a single pair_coeff command is used with the *pace* style which
specifies an ACE coefficient file followed by N additional arguments
specifying the mapping of ACE elements to LAMMPS atom types,
where N is the number of LAMMPS atom types:

* ACE coefficient file
* N element names = mapping of ACE elements to atom types

Only a single pair_coeff command is used with the *pace* style which
specifies an ACE file that fully defines the potential.
Note that unlike for other potentials, cutoffs are
not set in the pair_style or pair_coeff command; they are specified in
the ACE file.

The pair_style *pace* command may be followed by an optional keyword
*product* or *recursive*, which determines which of two algorithms
is used for the calculation of basis functions and derivatives.
The default is *recursive*.

See the :doc:`pair_coeff <pair_coeff>` page for alternate ways
to specify the path for the ACE coefficient file.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS with
user-specifiable parameters as described above.  You never need to
specify a pair_coeff command with I != J arguments for this style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the ML-PACE package.  It is only enabled if LAMMPS
was built with that package.
See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_style snap  <pair_snap>`

Default
"""""""

recursive

.. _Drautz20191:

**(Drautz)** Drautz, Phys Rev B, 99, 014104 (2019).

.. _Lysogorskiy20211:

**(Lysogorskiy)** Lysogorskiy, van der Oord, Bochkarev, Menon, Rinaldi, Hammerschmidt, Mrovec, Thompson, Csanyi, Ortner, Drautz, TBD (2021).
