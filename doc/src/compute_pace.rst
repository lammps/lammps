.. index:: compute pace

compute pace command
========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID pace ace_potential_filename ... keyword values ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* pace = style name of this compute command
* ace_potential_filename = file name (in the .yace or .ace format from :doc:`pace pair_style <pair_pace>`) including ACE hyper-parameters, bonds, and generalized coupling coefficients
* keyword = *bikflag* or *dgradflag*

  .. parsed-literal::

       *bikflag* value = *0* or *1*
          *0* = descriptors are summed over atoms of each type
          *1* = descriptors are listed separately for each atom
       *dgradflag* value = *0* or *1*
          *0* = descriptor gradients are summed over atoms of each type
          *1* = descriptor gradients are listed separately for each atom pair

Examples
""""""""

.. code-block:: LAMMPS

   compute pace all pace coupling_coefficients.yace
   compute pace all pace coupling_coefficients.yace 0 1
   compute pace all pace coupling_coefficients.yace 1 1

Description
"""""""""""

.. versionadded:: TBD

This compute calculates a set of quantities related to the atomic
cluster expansion (ACE) descriptors of the atoms in a group.  ACE
descriptors are highly general atomic descriptors, encoding the radial
and angular distribution of neighbor atoms, up to arbitrary bond order
(rank).  The detailed mathematical definition is given in the paper by
:ref:`(Drautz) <Drautz19>`.  These descriptors are used in the
:doc:`pace pair_style <pair_pace>`.  Quantities obtained from `compute
pace` are related to those used in :doc:`pace pair_style <pair_pace>` to
evaluate atomic energies, forces, and stresses for linear ACE models.

For example, the energy for a linear ACE model is calculated as:
:math:`E=\sum_i^{N\_atoms} \sum_{\boldsymbol{\nu}} c_{\boldsymbol{\nu}}
B_{i,\boldsymbol{\boldsymbol{\nu}}}`.  The ACE descriptors for atom `i`
:math:`B_{i,\boldsymbol{\nu}}`, and :math:`c_{\nu}` are linear model
parameters.  The detailed definition and indexing convention for ACE
descriptors is given in :ref:`(Drautz) <Drautz19>`.  In short, body
order :math:`N`, angular character, radial character, and chemical
elements in the *N-body* descriptor are encoded by :math:`\nu`.  In the
:doc:`pace pair_style <pair_pace>`, the linear model parameters and the
ACE descriptors are combined for efficient evaluation of energies and
forces.  The details and benefits of this efficient implementation are
given in :ref:`(Lysogorskiy) <Lysogorskiy21>`, but the combined
descriptors and linear model parameters for the purposes of `compute
pace` may be expressed in terms of the ACE descriptors mentioned above.

:math:`c_{\boldsymbol{\nu}} B_{i,\boldsymbol{\nu}}= \sum_{\boldsymbol{\nu}' \in \boldsymbol{\nu} } \big[ c_{\boldsymbol{\nu}} C(\boldsymbol{\nu}') \big] A_{i,\boldsymbol{\nu}'}`

where the bracketed terms on the right-hand side are the combined functions
with linear model parameters typically provided in the `<name>.yace` potential
file for `pace pair_style`. When these bracketed terms are multiplied by the
products of the atomic base from :ref:`(Drautz) <Drautz19>`,
:math:`A_{i,\boldsymbol{\nu'}}`, the ACE descriptors are recovered but they
are also scaled by linear model parameters. The generalized coupling coefficients,
written in short-hand here as :math:`C(\boldsymbol{\nu}')`, are the generalized
Clebsch-Gordan or generalized Wigner symbols. It may be desirable to reverse the
combination of these descriptors and the linear model parameters so that the
ACE descriptors themselves may be used. The ACE descriptors and their gradients
are often used when training ACE models, performing custom data analysis,
generalizing ACE model forms, and other tasks that involve direct computation of
descriptors. The key utility of `compute pace` is that it can compute the ACE
descriptors and gradients so that these tasks can be performed during a LAMMPS
simulation or so that LAMMPS can be used as a driver for tasks like ACE model
parameterization. To see how this command can be used within a Python workflow
to train ACE potentials, see the examples in
`FitSNAP <https://github.com/FitSNAP/FitSNAP>`_. Examples on using outputs from
this compute to construct general ACE potential forms are demonstrated in
:ref:`(Goff) <Goff23>`. The various keywords and inputs to `compute pace`
determine what ACE descriptors and related quantities are returned in a compute
array.

The coefficient file, `<name>.yace`, ultimately defines the number of ACE
descriptors to be computed, their maximum body-order, the degree of angular
character they have, the degree of radial character they have, the chemical
character (which element-element interactions are encoded by descriptors),
and other hyper-parameters defined in :ref:`(Drautz) <Drautz19>`. These may
be modeled after the potential files in :doc:`pace pair_style <pair_pace>`,
and have the same format. Details on how to generate the coefficient files
to train ACE models may be found in `FitSNAP <https://github.com/FitSNAP/FitSNAP>`_.

The keyword *bikflag* determines whether or not to list the descriptors of
each atom separately, or sum them together and list in a single row. If
*bikflag* is set to *0* then a single descriptor row is used, which contains
the per-atom ACE descriptors :math:`B_{i,\boldsymbol{\nu}}` summed over all
atoms *i* to produce :math:`B_{\boldsymbol{\nu}}`. If *bikflag* is set to
*1* this is replaced by a separate per-atom ACE descriptor row for each atom.
In this case, the entries in the final column for these rows are set to zero.

The keyword *dgradflag* determines whether to sum atom gradients or list
them separately. If *dgradflag* is set to 0, the ACE
descriptor gradients w.r.t. atom *j* are summed over all atoms *i'*
of, which may be useful when training linear ACE models on atomic forces.
If *dgradflag* is set to 1, gradients are listed separately for each pair of atoms.
Each row corresponds
to a single term :math:`\frac{\partial {B_{i,\boldsymbol{\nu}}}}{\partial {r}^a_j}`
where :math:`{r}^a_j` is the *a-th* position coordinate of the atom with global
index *j*. This also changes the number of columns to be equal to the number of
ACE descriptors, with 3 additional columns representing the indices :math:`i`,
:math:`j`, and :math:`a`, as explained more in the Output info section below.
The option *dgradflag=1* requires that *bikflag=1*.

.. note::

    It is noted here that in contrast to :doc:`pace pair_style <pair_pace>`,
    the *.yace* file for `compute pace` typically should not contain linear
    parameters for an ACE potential. If :math:`c_{\nu}` are included,
    the value of the descriptor will not be returned in the `compute` array,
    but instead, the energy contribution from that descriptor will be returned.
    Do not do this unless it is the desired behavior.
    *In short, you should not plug in a '.yace' for a pace potential into this
    compute to evaluate descriptors.*

.. note::

    *Generalized Clebsch-Gordan or Generalized Wigner symbols (with appropriate
    factors) must be used to evaluate ACE descriptors with this compute.* There
    are multiple ways to define the generalized coupling coefficients. Because
    of this, this compute will not revert your potential file to a coupling
    coefficient file. Instead this compute allows the user to supply coupling
    coefficients that follow any convention.

.. note::

   Using *dgradflag* = 1 produces a global array with :math:`N + 3N^2 + 1` rows
   which becomes expensive for systems with more than 1000 atoms.

.. note::

   If you have a bonded system, then the settings of :doc:`special_bonds
   <special_bonds>` command can remove pairwise interactions between
   atoms in the same bond, angle, or dihedral.  This is the default
   setting for the :doc:`special_bonds <special_bonds>` command, and
   means those pairwise interactions do not appear in the neighbor list.
   Because this fix uses the neighbor list, it also means those pairs
   will not be included in the calculation.  One way to get around this,
   is to write a dump file, and use the :doc:`rerun <rerun>` command to
   compute the ACE descriptors for snapshots in the dump file.
   The rerun script can use a :doc:`special_bonds <special_bonds>`
   command that includes all pairs in the neighbor list.

----------

Output info
"""""""""""

Compute *pace* evaluates a global array.  The columns are arranged into
*ntypes* blocks, listed in order of atom type *I*\ . Each block contains
one column for each ACE descriptor, the same as for compute
*sna/atom*\ in :doc:`compute snap <compute_sna_atom>`. A final column contains the corresponding energy, force
component on an atom, or virial stress component. The rows of the array
appear in the following order:

* 1 row: *pace* average descriptor values for all atoms of type *I*
* 3\*\ *n* force rows: quantities, with derivatives w.r.t. x, y, and z coordinate of atom *i* appearing in consecutive rows. The atoms are sorted based on atom ID and run up to the total number of atoms, *n*.
* 6 rows: *virial* quantities summed for all atoms of type *I*

For example, if :math:`\# \; B_{i, \boldsymbol{\nu}}` =30 and ntypes=1, the number of columns in the
The number of columns in the global array generated by *pace* are 31, and
931, respectively, while the number of rows is 1+3\*\ *n*\ +6, where *n*
is the total number of atoms.

If the *bik* keyword is set to 1, the structure of the pace array is expanded.
The first :math:`N` rows of the pace array
correspond to :math:`\# \; B_{i,\boldsymbol{\nu}}` instead of a single row summed over atoms :math:`i`.
In this case, the entries in the final column for these rows
are set to zero. Also, each row contains only non-zero entries for the
columns corresponding to the type of that atom. This is not true in the case
of *dgradflag* keyword = 1 (see below).

If the *dgradflag* keyword is set to 1, this changes the structure of the
global array completely.
Here the per-atom quantities are replaced with rows corresponding to
descriptor gradient components on single atoms:

.. math::

  \frac{\partial {B_{i,\boldsymbol{\nu}}  }}{\partial {r}^a_j}

where :math:`{r}^a_j` is the *a-th* position coordinate of the atom with global
index *j*. The rows are
organized in chunks, where each chunk corresponds to an atom with global index
:math:`j`. The rows in an atom :math:`j` chunk correspond to
atoms with global index :math:`i`. The total number of rows for
these descriptor gradients is therefore :math:`3N^2`.
The number of columns is equal to the number of ACE descriptors,
plus 3 additional left-most columns representing the global atom indices
:math:`i`, :math:`j`,
and Cartesian direction :math:`a`  (0, 1, 2, for x, y, z).
The first 3 columns of the first :math:`N` rows belong to the reference
potential force components. The remaining K columns contain the
:math:`B_{i,\boldsymbol{\nu}}` per-atom descriptors corresponding to the non-zero entries
obtained when *bikflag* = 1.
The first column of the last row, after the first
:math:`N + 3N^2` rows, contains the reference potential
energy. The virial components are not used with this option. The total number of
rows is therefore :math:`N + 3N^2 + 1` and the number of columns is :math:`K + 3`.

These values can be accessed by any command that uses global values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

Restrictions
""""""""""""

These computes are part of the ML-PACE package.  They are only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_style pace <pair_pace>`
:doc:`pair_style snap <pair_snap>`
:doc:`compute snap <compute_sna_atom>`

Default
"""""""

The optional keyword defaults are *bikflag* = 0,
*dgradflag* = 0

----------

.. _Drautz19:

**(Drautz)** Drautz, Phys Rev B, 99, 014104 (2019).

.. _Lysogorskiy21:

**(Lysogorskiy)** Lysogorskiy, van der Oord, Bochkarev, Menon, Rinaldi, Hammerschmidt, Mrovec, Thompson, Csanyi, Ortner, Drautz, npj Comp Mat, 7, 97 (2021).

.. _Goff23:

**(Goff)** Goff, Zhang, Negre, Rohskopf, Niklasson, Journal of Chemical Theory and Computation 19, no. 13 (2023).
