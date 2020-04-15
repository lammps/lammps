.. index:: compute orientorder/atom

compute orientorder/atom command
================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID orientorder/atom keyword values ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* orientorder/atom = style name of this compute command
* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *cutoff* or *nnn* or *degrees* or *components*
       *cutoff* value = distance cutoff
       *nnn* value = number of nearest neighbors
       *degrees* values = nlvalues, l1, l2,...
       *wl* value = yes or no
       *wl/hat* value = yes or no
       *components* value = ldegree

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all orientorder/atom
   compute 1 all orientorder/atom degrees 5 4 6 8 10 12 nnn NULL cutoff 1.5
   compute 1 all orientorder/atom wl/hat yes
   compute 1 all orientorder/atom components 6

Description
"""""""""""

Define a computation that calculates a set of bond-orientational
order parameters :math:`Q_l` for each atom in a group. These order parameters
were introduced by :ref:`Steinhardt et al. <Steinhardt>` as a way to
characterize the local orientational order in atomic structures.
For each atom, :math:`Q_l` is a real number defined as follows:

.. math::

   \bar{Y}_{lm} = & \frac{1}{nnn}\sum_{j = 1}^{nnn} Y_{lm}( \theta( {\bf r}_{ij} ), \phi( {\bf r}_{ij} ) ) \\
   Q_l = & \sqrt{\frac{4 \pi}{2 l + 1} \sum_{m = -l}^{m = l} \bar{Y}_{lm} \bar{Y}^*_{lm}}

The first equation defines the spherical harmonic order parameters.
These are complex number components of the 3D analog of the 2D order
parameter :math:`q_n`, which is implemented as LAMMPS compute
:doc:`hexorder/atom <compute_hexorder_atom>`.
The summation is over the *nnn* nearest
neighbors of the central atom.
The angles theta and phi are the standard spherical polar angles
defining the direction of the bond vector :math:`r_{ij}`.
The second equation defines :math:`Q_l`, which is a
rotationally invariant non-negative amplitude obtained by summing
over all the components of degree *l*\ .

The optional keyword *cutoff* defines the distance cutoff
used when searching for neighbors. The default value, also
the maximum allowable value, is the cutoff specified
by the pair style.

The optional keyword *nnn* defines the number of nearest
neighbors used to calculate :math:`Q_l`. The default value is 12.
If the value is NULL, then all neighbors up to the
specified distance cutoff are used.

The optional keyword *degrees* defines the list of order parameters to
be computed.  The first argument *nlvalues* is the number of order
parameters. This is followed by that number of non-negative integers giving the
degree of each order parameter. Because :math:`Q_2` and all odd-degree order
parameters are zero for atoms in cubic crystals (see
:ref:`Steinhardt <Steinhardt>`), the default order parameters are :math:`Q_4`,
:math:`Q_6`, :math:`Q_8`, :math:`Q_{10}`, and :math:`Q_{12}`. For the FCC
crystal with *nnn* =12, :math:`Q_4 = \sqrt{\frac{7}{192}} = 0.19094...`.
The numerical values of all order
parameters up to :math:`Q_12` for a range of commonly encountered
high-symmetry structures are given in Table I of :ref:`Mickel et al. <Mickel>`,
and these can be reproduced with this compute.

The optional keyword *wl* will output the third-order invariants :math:`W_l`
(see Eq. 1.4 in :ref:`Steinhardt <Steinhardt>`) for the same degrees as
for the :math:`Q_l` parameters. For the FCC crystal with *nnn* =12,
:math:`W_4` = -sqrt(14/143).(49/4096)/Pi\^1.5 = -0.0006722136...

The optional keyword *wl/hat* will output the normalized third-order
invariants :math:`\hat{W}_l` (see Eq. 2.2 in :ref:`Steinhardt <Steinhardt>`)
for the same degrees as for the :math:`Q_l` parameters. For the FCC crystal
with *nnn* =12, :math:`\hat{W}_4 = -\frac{7}{3} \sqrt{\frac{2}{429}} = -0.159317...`
The numerical
values of :math:`\hat{W}_l` for a range of commonly encountered high-symmetry
structures are given in Table I of :ref:`Steinhardt <Steinhardt>`, and these
can be reproduced with this keyword.

The optional keyword *components* will output the components of the
normalized complex vector :math:`\bar{Y}_{lm}` of degree *ldegree*\ , which must be
explicitly included in the keyword *degrees*\ . This option can be used
in conjunction with :doc:`compute coord_atom <compute_coord_atom>` to
calculate the ten Wolde's criterion to identify crystal-like
particles, as discussed in :ref:`ten Wolde <tenWolde2>`.

The value of :math:`Q_l` is set to zero for atoms not in the
specified compute group, as well as for atoms that have less than
*nnn* neighbors within the distance cutoff, unless *nnn* is NULL.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e. each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently.

.. note::

   If you have a bonded system, then the settings of
   :doc:`special_bonds <special_bonds>` command can remove pairwise
   interactions between atoms in the same bond, angle, or dihedral.  This
   is the default setting for the :doc:`special_bonds <special_bonds>`
   command, and means those pairwise interactions do not appear in the
   neighbor list.  Because this fix uses the neighbor list, it also means
   those pairs will not be included in the order parameter.  This
   difficulty can be circumvented by writing a dump file, and using the
   :doc:`rerun <rerun>` command to compute the order parameter for
   snapshots in the dump file.  The rerun script can use a
   :doc:`special_bonds <special_bonds>` command that includes all pairs in
   the neighbor list.

**Output info:**

This compute calculates a per-atom array with *nlvalues* columns,
giving the :math:`Q_l` values for each atom, which are real numbers on the
range :math:`0 <= Q_l <= 1`.

If the keyword *wl* is set to yes, then the :math:`W_l` values for each
atom will be added to the output array, which are real numbers.

If the keyword *wl/hat* is set to yes, then the :math:`\hat{W}_l`
values for each atom will be added to the output array, which are real numbers.

If the keyword *components* is set, then the real and imaginary parts
of each component of (normalized) :math:`\bar{Y}_{lm}` will be added to the
output array in the following order: :math:`Re(\bar{Y}_{-m}) Im(\bar{Y}_{-m})
Re(\bar{Y}_{-m+1}) Im(\bar{Y}_{-m+1}) ... Re(\bar{Y}_m) Im(\bar{Y}_m)`.  This
way, the per-atom array will have a total of *nlvalues*\ +2\*(2\ *l*\ +1)
columns.

These values can be accessed by any command that uses per-atom values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute coord/atom <compute_coord_atom>`, :doc:`compute centro/atom <compute_centro_atom>`, :doc:`compute hexorder/atom <compute_hexorder_atom>`

Default
"""""""

The option defaults are *cutoff* = pair style cutoff, *nnn* = 12,
*degrees* = 5 4 6 8 10 12 i.e. :math:`Q_4`, :math:`Q_6`, :math:`Q_8`, :math:`Q_{10}`, and :math:`Q_{12}`,
*wl* = no, *wl/hat* = no, and *components* off

----------

.. _Steinhardt:

**(Steinhardt)** P. Steinhardt, D. Nelson, and M. Ronchetti,
Phys. Rev. B 28, 784 (1983).

.. _Mickel:

**(Mickel)** W. Mickel, S. C. Kapfer, G. E. Schroeder-Turkand, K. Mecke,
J. Chem. Phys. 138, 044501 (2013).

.. _tenWolde2:

**(tenWolde)** P. R. ten Wolde, M. J. Ruiz-Montero, D. Frenkel,
J. Chem. Phys. 104, 9932 (1996).
