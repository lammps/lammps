.. index:: pair\_style snap

pair\_style snap command
========================

pair\_style snap/kk command
===========================

Syntax
""""""


.. parsed-literal::

   pair_style snap

Examples
""""""""


.. parsed-literal::

   pair_style snap
   pair_coeff \* \* InP.snapcoeff InP.snapparam In In P P

Description
"""""""""""

Pair style *snap* computes interactions using the spectral
neighbor analysis potential (SNAP) :ref:`(Thompson) <Thompson20142>`.
Like the GAP framework of Bartok et al. :ref:`(Bartok2010) <Bartok20102>`,
:ref:`(Bartok2013) <Bartok2013>` which uses bispectrum components
to characterize the local neighborhood of each atom
in a very general way. The mathematical definition of the
bispectrum calculation used by SNAP is identical
to that used by :doc:`compute sna/atom <compute_sna_atom>`.
In SNAP, the total energy is decomposed into a sum over
atom energies. The energy of atom *i* is
expressed as a weighted sum over bispectrum components.

.. math::

   E^i_{SNAP}(B_1^i,...,B_K^i) = \beta^{\alpha_i}_0 + \sum_{k=1}^K \beta_k^{\alpha_i} B_k^i


where :math:`B_k^i` is the *k*\ -th bispectrum component of atom *i*\ ,
and :math:`\beta_k^{\alpha_i}` is the corresponding linear coefficient
that depends on :math:\alpha_i`, the SNAP element of atom *i*\ . The
number of bispectrum components used and their definitions
depend on the value of *twojmax*
defined in the SNAP parameter file described below.
The bispectrum calculation is described in more detail
in :doc:`compute sna/atom <compute_sna_atom>`.

Note that unlike for other potentials, cutoffs for SNAP potentials are
not set in the pair\_style or pair\_coeff command; they are specified in
the SNAP potential files themselves.

Only a single pair\_coeff command is used with the *snap* style which
specifies a SNAP coefficient file followed by a SNAP parameter file
and then N additional arguments specifying the mapping of SNAP
elements to LAMMPS atom types, where N is the number of
LAMMPS atom types:

* SNAP coefficient file
* SNAP parameter file
* N element names = mapping of SNAP elements to atom types

As an example, if a LAMMPS indium phosphide simulation has 4 atoms
types, with the first two being indium and the 3rd and 4th being
phophorous, the pair\_coeff command would look like this:


.. parsed-literal::

   pair_coeff \* \* snap InP.snapcoeff InP.snapparam In In P P

The 1st 2 arguments must be \* \* so as to span all LAMMPS atom types.
The two filenames are for the coefficient and parameter files, respectively.
The two trailing 'In' arguments map LAMMPS atom types 1 and 2 to the
SNAP 'In' element. The two trailing 'P' arguments map LAMMPS atom types
3 and 4 to the SNAP 'P' element.

If a SNAP mapping value is
specified as NULL, the mapping is not performed.
This can be used when a *snap* potential is used as part of the
*hybrid* pair style.  The NULL values are placeholders for atom types
that will be used with other potentials.

The name of the SNAP coefficient file usually ends in the
".snapcoeff" extension. It may contain coefficients
for many SNAP elements. The only requirement is that it
contain at least those element names appearing in the
LAMMPS mapping list.
The name of the SNAP parameter file usually ends in the ".snapparam"
extension. It contains a small number
of parameters that define the overall form of the SNAP potential.
See the :doc:`pair_coeff <pair_coeff>` doc page for alternate ways
to specify the path for these files.

Quite commonly,
SNAP potentials are combined with one or more other LAMMPS pair styles
using the *hybrid/overlay* pair style. As an example, the SNAP
tantalum potential provided in the LAMMPS potentials directory
combines the *snap* and *zbl* pair styles. It is invoked
by the following commands:


.. parsed-literal::

           variable zblcutinner equal 4
           variable zblcutouter equal 4.8
           variable zblz equal 73
           pair_style hybrid/overlay &
           zbl ${zblcutinner} ${zblcutouter} snap
           pair_coeff \* \* zbl 0.0
           pair_coeff 1 1 zbl ${zblz}
           pair_coeff \* \* snap Ta06A.snapcoeff Ta06A.snapparam Ta

It is convenient to keep these commands in a separate file that can
be inserted in any LAMMPS input script using the :doc:`include <include>`
command.

The top of the SNAP coefficient file can contain any number of blank and comment lines (start with #), but follows a strict
format after that. The first non-blank non-comment
line must contain two integers:

* nelem  = Number of elements
* ncoeff = Number of coefficients

This is followed by one block for each of the *nelem* elements.
The first line of each block contains three entries:

* Element symbol (text string)
* R = Element radius (distance units)
* w = Element weight (dimensionless)

This line is followed by *ncoeff* coefficients, one per line.

The SNAP parameter file can contain blank and comment lines (start
with #) anywhere. Each non-blank non-comment line must contain one
keyword/value pair. The required keywords are *rcutfac* and
*twojmax*\ . Optional keywords are *rfac0*\ , *rmin0*\ ,
*switchflag*\ , and *bzeroflag*\ .

The default values for these keywords are

* *rfac0* = 0.99363
* *rmin0* = 0.0
* *switchflag* = 0
* *bzeroflag* = 1
* *quadraticflag* = 1

Detailed definitions for all the keywords are given on the :doc:`compute sna/atom <compute_sna_atom>` doc page.
If *quadraticflag* is set to 1, then the SNAP energy expression includes the quadratic term,
0.5\*B\^t.alpha.B, where alpha is a symmetric *K* by *K* matrix.
The SNAP element file should contain *K*\ (\ *K*\ +1)/2 additional coefficients
for each element, the upper-triangular elements of alpha.

.. note::

   The previously used *diagonalstyle* keyword was removed in 2019,
   since all known SNAP potentials use the default value of 3.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS with
user-specifiable parameters as described above.  You never need to
specify a pair\_coeff command with I != J arguments for this style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair\_style and pair\_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


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


Restrictions
""""""""""""


This style is part of the SNAP package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`compute sna/atom <compute_sna_atom>`,
:doc:`compute snad/atom <compute_sna_atom>`,
:doc:`compute snav/atom <compute_sna_atom>`

**Default:** none


----------


.. _Thompson20142:



**(Thompson)** Thompson, Swiler, Trott, Foiles, Tucker, J Comp Phys, 285, 316 (2015).

.. _Bartok20102:



**(Bartok2010)** Bartok, Payne, Risi, Csanyi, Phys Rev Lett, 104, 136403 (2010).

.. _Bartok2013:



**(Bartok2013)** Bartok, Gillan, Manby, Csanyi, Phys Rev B 87, 184115 (2013).
