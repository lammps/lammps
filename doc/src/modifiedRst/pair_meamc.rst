.. index:: pair\_style meam/c

pair\_style meam/c command
==========================

Syntax
""""""


.. parsed-literal::

   pair_style style

style = *meam/c*

Examples
""""""""


.. parsed-literal::

   pair_style meam/c
   pair_coeff \* \* ../potentials/library.meam Si ../potentials/si.meam Si
   pair_coeff \* \* ../potentials/library.meam Ni Al NULL Ni Al Ni Ni

Description
"""""""""""

.. note::

   The behavior of the MEAM potential for alloy systems has changed
   as of November 2010; see description below of the mixture\_ref\_t
   parameter

Style *meam/c* computes pairwise interactions for a variety of materials
using modified embedded-atom method (MEAM) potentials
:ref:`(Baskes) <Baskes>`.  Conceptually, it is an extension to the original
:doc:`EAM potentials <pair_eam>` which adds angular forces.  It is
thus suitable for modeling metals and alloys with fcc, bcc, hcp and
diamond cubic structures, as well as covalently bonded materials like
silicon and carbon. Style *meam/c* is a translation of the (now obsolete)
*meam* code from Fortran to C++. It is functionally equivalent to *meam*
but more efficient, and thus *meam* has been removed from LAMMPS after
the 12 December 2018 release.

In the MEAM formulation, the total energy E of a system of atoms is
given by:

.. math::

  E = \sum_i \left\{ F_i(\bar{\rho}_i)
      + \frac{1}{2} \sum_{i \neq j} \phi_{ij} (r_{ij}) \right\}


where F is the embedding energy which is a function of the atomic
electron density rho, and phi is a pair potential interaction.  The
pair interaction is summed over all neighbors J of atom I within the
cutoff distance.  As with EAM, the multi-body nature of the MEAM
potential is a result of the embedding energy term.  Details of the
computation of the embedding and pair energies, as implemented in
LAMMPS, are given in :ref:`(Gullet) <Gullet>` and references therein.

The various parameters in the MEAM formulas are listed in two files
which are specified by the :doc:`pair\_coeff <pair_coeff>` command.
These are ASCII text files in a format consistent with other MD codes
that implement MEAM potentials, such as the serial DYNAMO code and
Warp.  Several MEAM potential files with parameters for different
materials are included in the "potentials" directory of the LAMMPS
distribution with a ".meam" suffix.  All of these are parameterized in
terms of LAMMPS :doc:`metal units <units>`.

Note that unlike for other potentials, cutoffs for MEAM potentials are
not set in the pair\_style or pair\_coeff command; they are specified in
the MEAM potential files themselves.

Only a single pair\_coeff command is used with the *meam* style which
specifies two MEAM files and the element(s) to extract information
for.  The MEAM elements are mapped to LAMMPS atom types by specifying
N additional arguments after the 2nd filename in the pair\_coeff
command, where N is the number of LAMMPS atom types:

* MEAM library file
* Elem1, Elem2, ...
* MEAM parameter file
* N element names = mapping of MEAM elements to atom types

See the :doc:`pair\_coeff <pair_coeff>` doc page for alternate ways
to specify the path for the potential files.

As an example, the potentials/library.meam file has generic MEAM
settings for a variety of elements.  The potentials/SiC.meam file has
specific parameter settings for a Si and C alloy system.  If your
LAMMPS simulation has 4 atoms types and you want the 1st 3 to be Si,
and the 4th to be C, you would use the following pair\_coeff command:


.. parsed-literal::

   pair_coeff \* \* library.meam Si C sic.meam Si Si Si C

The 1st 2 arguments must be \* \* so as to span all LAMMPS atom types.
The two filenames are for the library and parameter file respectively.
The Si and C arguments (between the file names) are the two elements
for which info will be extracted from the library file.  The first
three trailing Si arguments map LAMMPS atom types 1,2,3 to the MEAM Si
element.  The final C argument maps LAMMPS atom type 4 to the MEAM C
element.

If the 2nd filename is specified as NULL, no parameter file is read,
which simply means the generic parameters in the library file are
used.  Use of the NULL specification for the parameter file is
discouraged for systems with more than a single element type
(e.g. alloys), since the parameter file is expected to set element
interaction terms that are not captured by the information in the
library file.

If a mapping value is specified as NULL, the mapping is not performed.
This can be used when a *meam* potential is used as part of the
*hybrid* pair style.  The NULL values are placeholders for atom types
that will be used with other potentials.

.. note::

   If the 2nd filename is NULL, the element names between the two
   filenames can appear in any order, e.g. "Si C" or "C Si" in the
   example above.  However, if the 2nd filename is not NULL (as in the
   example above), it contains settings that are Fortran-indexed for the
   elements that preceed it.  Thus you need to insure you list the
   elements between the filenames in an order consistent with how the
   values in the 2nd filename are indexed.  See details below on the
   syntax for settings in the 2nd file.

The MEAM library file provided with LAMMPS has the name
potentials/library.meam.  It is the "meamf" file used by other MD
codes.  Aside from blank and comment lines (start with #) which can
appear anywhere, it is formatted as a series of entries, each of which
has 19 parameters and can span multiple lines:

elt, lat, z, ielement, atwt, alpha, b0, b1, b2, b3, alat, esub, asub,
t0, t1, t2, t3, rozero, ibar

The "elt" and "lat" parameters are text strings, such as elt = Si or
Cu and lat = dia or fcc.  Because the library file is used by Fortran
MD codes, these strings may be enclosed in single quotes, but this is
not required.  The other numeric parameters match values in the
formulas above.  The value of the "elt" string is what is used in the
pair\_coeff command to identify which settings from the library file
you wish to read in.  There can be multiple entries in the library
file with the same "elt" value; LAMMPS reads the 1st matching entry it
finds and ignores the rest.

Other parameters in the MEAM library file correspond to single-element
potential parameters:


.. parsed-literal::

   lat      = lattice structure of reference configuration
   z        = number of nearest neighbors in the reference structure
   ielement = atomic number
   atwt     = atomic weight
   alat     = lattice constant of reference structure
   esub     = energy per atom (eV) in the reference structure at equilibrium
   asub     = "A" parameter for MEAM (see e.g. :ref:`(Baskes) <Baskes>`)

The alpha, b0, b1, b2, b3, t0, t1, t2, t3 parameters correspond to the
standard MEAM parameters in the literature :ref:`(Baskes) <Baskes>` (the b
parameters are the standard beta parameters).  The rozero parameter is
an element-dependent density scaling that weights the reference
background density (see e.g. equation 4.5 in :ref:`(Gullet) <Gullet>`) and
is typically 1.0 for single-element systems.  The ibar parameter
selects the form of the function G(Gamma) used to compute the electron
density; options are


.. parsed-literal::

      0 => G = sqrt(1+Gamma)
      1 => G = exp(Gamma/2)
      2 => not implemented
      3 => G = 2/(1+exp(-Gamma))
      4 => G = sqrt(1+Gamma)
     -5 => G = +-sqrt(abs(1+Gamma))

If used, the MEAM parameter file contains settings that override or
complement the library file settings.  Examples of such parameter
files are in the potentials directory with a ".meam" suffix.  Their
format is the same as is read by other Fortran MD codes.  Aside from
blank and comment lines (start with #) which can appear anywhere, each
line has one of the following forms.  Each line can also have a
trailing comment (starting with #) which is ignored.


.. parsed-literal::

   keyword = value
   keyword(I) = value
   keyword(I,J) = value
   keyword(I,J,K) = value

The indices I, J, K correspond to the elements selected from the
MEAM library file numbered in the order of how those elements were
selected starting from 1. Thus for the example given below


.. parsed-literal::

   pair_coeff \* \* library.meam Si C sic.meam Si Si Si C

an index of 1 would refer to Si and an index of 2 to C.

The recognized keywords for the parameter file are as follows:

Ec, alpha, rho0, delta, lattce, attrac, repuls, nn2, Cmin, Cmax, rc, delr,
augt1, gsmooth\_factor, re

where


.. parsed-literal::

   rc          = cutoff radius for cutoff function; default = 4.0
   delr        = length of smoothing distance for cutoff function; default = 0.1
   rho0(I)     = relative density for element I (overwrites value
                 read from meamf file)
   Ec(I,J)     = cohesive energy of reference structure for I-J mixture
   delta(I,J)  = heat of formation for I-J alloy; if Ec_IJ is input as
                 zero, then LAMMPS sets Ec_IJ = (Ec_II + Ec_JJ)/2 - delta_IJ
   alpha(I,J)  = alpha parameter for pair potential between I and J (can
                 be computed from bulk modulus of reference structure
   re(I,J)     = equilibrium distance between I and J in the reference
                 structure
   Cmax(I,J,K) = Cmax screening parameter when I-J pair is screened
                 by K (I<=J); default = 2.8
   Cmin(I,J,K) = Cmin screening parameter when I-J pair is screened
                 by K (I<=J); default = 2.0
   lattce(I,J) = lattice structure of I-J reference structure:
                   dia = diamond (interlaced fcc for alloy)
                   fcc = face centered cubic
                   bcc = body centered cubic
                   dim = dimer
                   b1  = rock salt (NaCl structure)
                   hcp = hexagonal close-packed
                   c11 = MoSi2 structure
                   l12 = Cu3Au structure (lower case L, followed by 12)
                   b2  = CsCl structure (interpenetrating simple cubic)
   nn2(I,J)    = turn on second-nearest neighbor MEAM formulation for
                 I-J pair (see for example :ref:`(Lee) <Lee>`).
                   0 = second-nearest neighbor formulation off
                   1 = second-nearest neighbor formulation on
                   default = 0
   attrac(I,J) = additional cubic attraction term in Rose energy I-J pair potential
                   default = 0
   repuls(I,J) = additional cubic repulsive term in Rose energy I-J pair potential
                   default = 0
   zbl(I,J)    = blend the MEAM I-J pair potential with the ZBL potential for small
                 atom separations :ref:`(ZBL) <ZBL>`
                   default = 1
   gsmooth_factor  = factor determining the length of the G-function smoothing
                     region; only significant for ibar=0 or ibar=4.
                         99.0 = short smoothing region, sharp step
                         0.5  = long smoothing region, smooth step
                         default = 99.0
   augt1           = integer flag for whether to augment t1 parameter by
                     3/5\*t3 to account for old vs. new meam formulations;
                       0 = don't augment t1
                       1 = augment t1
                       default = 1
   ialloy          = integer flag to use alternative averaging rule for t parameters,
                     for comparison with the DYNAMO MEAM code
                       0 = standard averaging (matches ialloy=0 in DYNAMO)
                       1 = alternative averaging (matches ialloy=1 in DYNAMO)
                       2 = no averaging of t (use single-element values)
                       default = 0
   mixture_ref_t   = integer flag to use mixture average of t to compute the background
                     reference density for alloys, instead of the single-element values
                     (see description and warning elsewhere in this doc page)
                       0 = do not use mixture averaging for t in the reference density
                       1 = use mixture averaging for t in the reference density
                       default = 0
   erose_form      = integer value to select the form of the Rose energy function
                     (see description below).
                       default = 0
   emb_lin_neg     = integer value to select embedding function for negative densities
                       0 = F(rho)=0
                       1 = F(rho) = -asub\*esub\*rho (linear in rho, matches DYNAMO)
                       default = 0
   bkgd_dyn        = integer value to select background density formula
                       0 = rho_bkgd = rho_ref_meam(a) (as in the reference structure)
                       1 = rho_bkgd = rho0_meam(a)\*Z_meam(a) (matches DYNAMO)
                       default = 0

Rc, delr, re are in distance units (Angstroms in the case of metal
units).  Ec and delta are in energy units (eV in the case of metal
units).

Each keyword represents a quantity which is either a scalar, vector,
2d array, or 3d array and must be specified with the correct
corresponding array syntax.  The indices I,J,K each run from 1 to N
where N is the number of MEAM elements being used.

Thus these lines


.. parsed-literal::

   rho0(2) = 2.25
   alpha(1,2) = 4.37

set rho0 for the 2nd element to the value 2.25 and set alpha for the
alloy interaction between elements 1 and 2 to 4.37.

The augt1 parameter is related to modifications in the MEAM
formulation of the partial electron density function.  In recent
literature, an extra term is included in the expression for the
third-order density in order to make the densities orthogonal (see for
example :ref:`(Wang) <Wang2>`, equation 3d); this term is included in the
MEAM implementation in lammps.  However, in earlier published work
this term was not included when deriving parameters, including most of
those provided in the library.meam file included with lammps, and to
account for this difference the parameter t1 must be augmented by
3/5\*t3.  If augt1=1, the default, this augmentation is done
automatically.  When parameter values are fit using the modified
density function, as in more recent literature, augt1 should be set to
0.

The mixture\_ref\_t parameter is available to match results with those
of previous versions of lammps (before January 2011).  Newer versions
of lammps, by default, use the single-element values of the t
parameters to compute the background reference density.  This is the
proper way to compute these parameters.  Earlier versions of lammps
used an alloy mixture averaged value of t to compute the background
reference density.  Setting mixture\_ref\_t=1 gives the old behavior.
WARNING: using mixture\_ref\_t=1 will give results that are demonstrably
incorrect for second-neighbor MEAM, and non-standard for
first-neighbor MEAM; this option is included only for matching with
previous versions of lammps and should be avoided if possible.

The parameters attrac and repuls, along with the integer selection
parameter erose\_form, can be used to modify the Rose energy function
used to compute the pair potential.  This function gives the energy of
the reference state as a function of interatomic spacing.  The form of
this function is:


.. parsed-literal::

   astar = alpha \* (r/re - 1.d0)
   if erose_form = 0: erose = -Ec\*(1+astar+a3\*(astar\*\*3)/(r/re))\*exp(-astar)
   if erose_form = 1: erose = -Ec\*(1+astar+(-attrac+repuls/r)\*(astar\*\*3))\*exp(-astar)
   if erose_form = 2: erose = -Ec\*(1 +astar + a3\*(astar\*\*3))\*exp(-astar)
   a3 = repuls, astar < 0
   a3 = attrac, astar >= 0

Most published MEAM parameter sets use the default values attrac=repulse=0.
Setting repuls=attrac=delta corresponds to the form used in several
recent published MEAM parameter sets, such as :ref:`(Valone) <Valone>`

.. note::

   The default form of the erose expression in LAMMPS was corrected
   in March 2009.  The current version is correct, but may show different
   behavior compared with earlier versions of lammps with the attrac
   and/or repuls parameters are non-zero.  To obtain the previous default
   form, use erose\_form = 1 (this form does not seem to appear in the
   literature).  An alternative form (see e.g. :ref:`(Lee2) <Lee2>`) is
   available using erose\_form = 2.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS with
user-specifiable parameters as described above.  You never need to
specify a pair\_coeff command with I != J arguments for this style.

This pair style does not support the :doc:`pair\_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair\_style and pair\_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run\_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


The *meam/c* style is provided in the USER-MEAMC package. It is
only enabled if LAMMPS was built with that package.
See the :doc:`Build package <Build_package>` doc page for more info.

The maximum number of elements, that can be read from the MEAM
library file, is determined at compile time. The default is 5.
If you need support for more elements, you have to change the
define for the constant 'maxelt' at the beginning of the file
src/USER-MEAMC/meam.h and update/recompile LAMMPS. There is no
limit on the number of atoms types.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`, :doc:`pair\_style eam <pair_eam>`,
:doc:`pair\_style meam/spline <pair_meam_spline>`

**Default:** none


----------


.. _Baskes:



**(Baskes)** Baskes, Phys Rev B, 46, 2727-2742 (1992).

.. _Gullet:



**(Gullet)** Gullet, Wagner, Slepoy, SANDIA Report 2003-8782 (2003).
This report may be accessed on-line via `this link <sandreport_>`_.

.. _sandreport: http://infoserve.sandia.gov/sand\_doc/2003/038782.pdf



.. _Lee:



**(Lee)** Lee, Baskes, Phys. Rev. B, 62, 8564-8567 (2000).

.. _Lee2:



**(Lee2)** Lee, Baskes, Kim, Cho.  Phys. Rev. B, 64, 184102 (2001).

.. _Valone:



**(Valone)** Valone, Baskes, Martin, Phys. Rev. B, 73, 214209 (2006).

.. _Wang2:



**(Wang)** Wang, Van Hove, Ross, Baskes, J. Chem. Phys., 121, 5410 (2004).

.. _ZBL:



**(ZBL)** J.F. Ziegler, J.P. Biersack, U. Littmark, "Stopping and Ranges
of Ions in Matter", Vol 1, 1985, Pergamon Press.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
