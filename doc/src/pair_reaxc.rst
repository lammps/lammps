.. index:: pair\_style reax/c

pair\_style reax/c command
==========================

pair\_style reax/c/kk command
=============================

pair\_style reax/c/omp command
==============================

Syntax
""""""


.. parsed-literal::

   pair_style reax/c cfile keyword value

* cfile = NULL or name of a control file
* zero or more keyword/value pairs may be appended
  
  .. parsed-literal::
  
     keyword = *checkqeq* or *lgvdw* or *safezone* or *mincap*
       *checkqeq* value = *yes* or *no* = whether or not to require qeq/reax fix
       *enobonds* value = *yes* or *no* = whether or not to tally energy of atoms with no bonds
       *lgvdw* value = *yes* or *no* = whether or not to use a low gradient vdW correction
       *safezone* = factor used for array allocation
       *mincap* = minimum size for array allocation



Examples
""""""""


.. parsed-literal::

   pair_style reax/c NULL
   pair_style reax/c controlfile checkqeq no
   pair_style reax/c NULL lgvdw yes
   pair_style reax/c NULL safezone 1.6 mincap 100
   pair_coeff \* \* ffield.reax C H O N

Description
"""""""""""

Style *reax/c* computes the ReaxFF potential of van Duin, Goddard and
co-workers.  ReaxFF uses distance-dependent bond-order functions to
represent the contributions of chemical bonding to the potential
energy. There is more than one version of ReaxFF. The version
implemented in LAMMPS uses the functional forms documented in the
supplemental information of the following paper: :ref:`(Chenoweth et al., 2008) <Chenoweth_20082>`.  The version integrated into LAMMPS matches
the most up-to-date version of ReaxFF as of summer 2010.  For more
technical details about the pair reax/c implementation of ReaxFF, see
the :ref:`(Aktulga) <Aktulga>` paper. The *reax/c* style was initially
implemented as a stand-alone C code and is now integrated into LAMMPS
as a package.

The *reax/c/kk* style is a Kokkos version of the ReaxFF potential that
is derived from the *reax/c* style. The Kokkos version can run on GPUs
and can also use OpenMP multithreading. For more information about the
Kokkos package, see :doc:`Packages details <Packages_details>` and
:doc:`Speed kokkos <Speed_kokkos>` doc pages.  One important
consideration when using the *reax/c/kk* style is the choice of either
half or full neighbor lists. This setting can be changed using the
Kokkos :doc:`package <package>` command.

The *reax/c* style differs from the (obsolete) "pair\_style reax"
command in the implementation details.  The *reax* style was a
Fortran library, linked to LAMMPS.  The *reax* style has been removed
from LAMMPS after the 12 December 2018 version.

LAMMPS provides several different versions of ffield.reax in its
potentials dir, each called potentials/ffield.reax.label.  These are
documented in potentials/README.reax.  The default ffield.reax
contains parameterizations for the following elements: C, H, O, N.

The format of these files is identical to that used originally by van
Duin.  We have tested the accuracy of *pair\_style reax/c* potential
against the original ReaxFF code for the systems mentioned above.  You
can use other ffield files for specific chemical systems that may be
available elsewhere (but note that their accuracy may not have been
tested).

.. note::

   We do not distribute a wide variety of ReaxFF force field files
   with LAMMPS.  Adri van Duin's group at PSU is the central repository
   for this kind of data as they are continuously deriving and updating
   parameterizations for different classes of materials.  You can submit
   a contact request at the Materials Computation Center (MCC) website
   `https://www.mri.psu.edu/materials-computation-center/connect-mcc <https://www.mri.psu.edu/materials-computation-center/connect-mcc>`_,
   describing the material(s) you are interested in modeling with ReaxFF.
   They can tell you what is currently available or what it would take to
   create a suitable ReaxFF parameterization.

The *cfile* setting can be specified as NULL, in which case default
settings are used. A control file can be specified which defines
values of control variables. Some control variables are
global parameters for the ReaxFF potential. Others define certain
performance and output settings.
Each line in the control file specifies the value for
a control variable.  The format of the control file is described
below.

.. note::

   The LAMMPS default values for the ReaxFF global parameters
   correspond to those used by Adri van Duin's stand-alone serial
   code. If these are changed by setting control variables in the control
   file, the results from LAMMPS and the serial code will not agree.

Examples using *pair\_style reax/c* are provided in the examples/reax
sub-directory.

Use of this pair style requires that a charge be defined for every
atom.  See the :doc:`atom_style <atom_style>` and
:doc:`read_data <read_data>` commands for details on how to specify
charges.

The ReaxFF parameter files provided were created using a charge
equilibration (QEq) model for handling the electrostatic interactions.
Therefore, by default, LAMMPS requires that the :doc:`fix qeq/reax <fix_qeq_reax>` command be used with *pair\_style reax/c*
when simulating a ReaxFF model, to equilibrate charge each timestep.
Using the keyword *checkqeq* with the value *no*
turns off the check for *fix qeq/reax*\ ,
allowing a simulation to be run without charge equilibration.
In this case, the static charges you
assign to each atom will be used for computing the electrostatic
interactions in the system.
See the :doc:`fix qeq/reax <fix_qeq_reax>` command for details.

Using the optional keyword *lgvdw* with the value *yes* turns on the
low-gradient correction of the ReaxFF/C for long-range London
Dispersion, as described in the :ref:`(Liu) <Liu_2011>` paper. Force field
file *ffield.reax.lg* is designed for this correction, and is trained
for several energetic materials (see "Liu"). When using lg-correction,
recommended value for parameter *thb* is 0.01, which can be set in the
control file.  Note: Force field files are different for the original
or lg corrected pair styles, using wrong ffield file generates an
error message.

Using the optional keyword *enobonds* with the value *yes*\ , the energy
of atoms with no bonds (i.e. isolated atoms) is included in the total
potential energy and the per-atom energy of that atom.  If the value
*no* is specified then the energy of atoms with no bonds is set to
zero.  The latter behavior is usual not desired, as it causes
discontinuities in the potential energy when the bonding of an atom
drops to zero.

Optional keywords *safezone* and *mincap* are used for allocating
reax/c arrays.  Increasing these values can avoid memory problems,
such as segmentation faults and bondchk failed errors, that could
occur under certain conditions. These keywords aren't used by the
Kokkos version, which instead uses a more robust memory allocation
scheme that checks if the sizes of the arrays have been exceeded and
automatically allocates more memory.

The thermo variable *evdwl* stores the sum of all the ReaxFF potential
energy contributions, with the exception of the Coulombic and charge
equilibration contributions which are stored in the thermo variable
*ecoul*\ .  The output of these quantities is controlled by the
:doc:`thermo <thermo>` command.

This pair style tallies a breakdown of the total ReaxFF potential
energy into sub-categories, which can be accessed via the :doc:`compute pair <compute_pair>` command as a vector of values of length 14.
The 14 values correspond to the following sub-categories (the variable
names in italics match those used in the original FORTRAN ReaxFF
code):

1. *eb* = bond energy
2. *ea* = atom energy
3. *elp* = lone-pair energy
4. *emol* = molecule energy (always 0.0)
5. *ev* = valence angle energy
6. *epen* = double-bond valence angle penalty
7. *ecoa* = valence angle conjugation energy
8. *ehb* = hydrogen bond energy
9. *et* = torsion energy
10. *eco* = conjugation energy
11. *ew* = van der Waals energy
12. *ep* = Coulomb energy
13. *efi* = electric field energy (always 0.0)
14. *eqeq* = charge equilibration energy

To print these quantities to the log file (with descriptive column
headings) the following commands could be included in an input script:


.. parsed-literal::

   compute reax all pair reax/c
   variable eb      equal c_reax[1]
   variable ea      equal c_reax[2]
   [...]
   variable eqeq    equal c_reax[14]
   thermo_style custom step temp epair v_eb v_ea [...] v_eqeq

Only a single pair\_coeff command is used with the *reax/c* style which
specifies a ReaxFF potential file with parameters for all needed
elements.  These are mapped to LAMMPS atom types by specifying N
additional arguments after the filename in the pair\_coeff command,
where N is the number of LAMMPS atom types:

* filename
* N indices = ReaxFF elements

The filename is the ReaxFF potential file.

In the ReaxFF potential file, near the top, after the general
parameters, is the atomic parameters section that contains element
names, each with a couple dozen numeric parameters.  If there are M
elements specified in the *ffield* file, think of these as numbered 1
to M. Each of the N indices you specify for the N atom types of LAMMPS
atoms must be an integer from 1 to M.  Atoms with LAMMPS type 1 will
be mapped to whatever element you specify as the first index value,
etc.  If a mapping value is specified as NULL, the mapping is not
performed.  This can be used when the *reax/c* style is used as part
of the *hybrid* pair style.  The NULL values are placeholders for atom
types that will be used with other potentials.

As an example, say your LAMMPS simulation has 4 atom types and the
elements are ordered as C, H, O, N in the *ffield* file.  If you want
the LAMMPS atom type 1 and 2 to be C, type 3 to be N, and type 4 to be
H, you would use the following pair\_coeff command:


.. parsed-literal::

   pair_coeff \* \* ffield.reax C C N H


----------


The format of a line in the control file is as follows:


.. parsed-literal::

   variable_name value

and it may be followed by an "!" character and a trailing comment.

If the value of a control variable is not specified, then default
values are used.  What follows is the list of variables along with a
brief description of their use and default values.

simulation\_name: Output files produced by *pair\_style reax/c* carry
this name + extensions specific to their contents.  Partial energies
are reported with a ".pot" extension, while the trajectory file has
".trj" extension.

tabulate\_long\_range: To improve performance, long range interactions
can optionally be tabulated (0 means no tabulation). Value of this
variable denotes the size of the long range interaction table.  The
range from 0 to long range cutoff (defined in the *ffield* file) is
divided into *tabulate\_long\_range* points.  Then at the start of
simulation, we fill in the entries of the long range interaction table
by computing the energies and forces resulting from van der Waals and
Coulomb interactions between every possible atom type pairs present in
the input system.  During the simulation we consult to the long range
interaction table to estimate the energy and forces between a pair of
atoms. Linear interpolation is used for estimation. (default value =
0)

energy\_update\_freq: Denotes the frequency (in number of steps) of
writes into the partial energies file. (default value = 0)

nbrhood\_cutoff: Denotes the near neighbors cutoff (in Angstroms)
regarding the bonded interactions. (default value = 5.0)

hbond\_cutoff: Denotes the cutoff distance (in Angstroms) for hydrogen
bond interactions.(default value = 7.5. A value of 0.0 turns off
hydrogen bonds)

bond\_graph\_cutoff: is the threshold used in determining what is a
physical bond, what is not. Bonds and angles reported in the
trajectory file rely on this cutoff. (default value = 0.3)

thb\_cutoff: cutoff value for the strength of bonds to be considered in
three body interactions. (default value = 0.001)

thb\_cutoff\_sq: cutoff value for the strength of bond order products
to be considered in three body interactions. (default value = 0.00001)

write\_freq: Frequency of writes into the trajectory file. (default
value = 0)

traj\_title: Title of the trajectory - not the name of the trajectory
file.

atom\_info: 1 means print only atomic positions + charge (default = 0)

atom\_forces: 1 adds net forces to atom lines in the trajectory file
(default = 0)

atom\_velocities: 1 adds atomic velocities to atoms line (default = 0)

bond\_info: 1 prints bonds in the trajectory file (default = 0)

angle\_info: 1 prints angles in the trajectory file (default = 0)


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair style does not support the :doc:`pair_modify <pair_modify>`
mix, shift, table, and tail options.

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


This pair style is part of the USER-REAXC package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

The ReaxFF potential files provided with LAMMPS in the potentials
directory are parameterized for real :doc:`units <units>`.  You can use
the ReaxFF potential with any LAMMPS units, but you would need to
create your own potential file with coefficients listed in the
appropriate units if your simulation doesn't use "real" units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`fix qeq/reax <fix_qeq_reax>`, :doc:`fix reax/c/bonds <fix_reaxc_bonds>`, :doc:`fix reax/c/species <fix_reaxc_species>`

Default
"""""""

The keyword defaults are checkqeq = yes, enobonds = yes, lgvdw = no,
safezone = 1.2, mincap = 50.


----------


.. _Chenoweth\_20082:



**(Chenoweth\_2008)** Chenoweth, van Duin and Goddard,
Journal of Physical Chemistry A, 112, 1040-1053 (2008).

.. _Aktulga:



(Aktulga) Aktulga, Fogarty, Pandit, Grama, Parallel Computing, 38,
245-259 (2012).

.. _Liu\_2011:



**(Liu)** L. Liu, Y. Liu, S. V. Zybin, H. Sun and W. A. Goddard, Journal
of Physical Chemistry A, 115, 11016-11022 (2011).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
