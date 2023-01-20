Package details
===============

Here is a brief description of all packages in LAMMPS.  It lists authors
(if applicable) and summarizes the package contents.  It has specific
instructions on how to install the package, including, if necessary,
info on how to download or build any extra library it requires.  It also
gives links to documentation, example scripts, and pictures/movies (if
available) that illustrate use of the package.

The majority of packages can be included in a LAMMPS build with a
single setting (``-D PKG_<NAME>=on`` for CMake) or command
(``make yes-<name>`` for make).  See the :doc:`Build package <Build_package>`
page for more info.  A few packages may require additional steps;
this is indicated in the descriptions below.  The :doc:`Build extras <Build_extras>`
page gives those details.

.. note::

   To see the complete list of commands a package adds to LAMMPS,
   you can examine the files in its src directory, e.g. "ls
   src/GRANULAR".  Files with names that start with fix, compute, atom,
   pair, bond, angle, etc correspond to commands with the same style name
   as contained in the file name.

.. table_from_list::
   :columns: 6

   * :ref:`ADIOS <PKG-ADIOS>`
   * :ref:`AMOEBA <PKG-AMOEBA>`
   * :ref:`ASPHERE <PKG-ASPHERE>`
   * :ref:`ATC <PKG-ATC>`
   * :ref:`AWPMD <PKG-AWPMD>`
   * :ref:`BOCS <PKG-BOCS>`
   * :ref:`BODY <PKG-BODY>`
   * :ref:`BPM <PKG-BPM>`
   * :ref:`BROWNIAN <PKG-BROWNIAN>`
   * :ref:`CG-DNA <PKG-CG-DNA>`
   * :ref:`CG-SPICA <PKG-CG-SPICA>`
   * :ref:`CLASS2 <PKG-CLASS2>`
   * :ref:`COLLOID <PKG-COLLOID>`
   * :ref:`COLVARS <PKG-COLVARS>`
   * :ref:`COMPRESS <PKG-COMPRESS>`
   * :ref:`CORESHELL <PKG-CORESHELL>`
   * :ref:`DIELECTRIC <PKG-DIELECTRIC>`
   * :ref:`DIFFRACTION <PKG-DIFFRACTION>`
   * :ref:`DIPOLE <PKG-DIPOLE>`
   * :ref:`DPD-BASIC <PKG-DPD-BASIC>`
   * :ref:`DPD-MESO <PKG-DPD-MESO>`
   * :ref:`DPD-REACT <PKG-DPD-REACT>`
   * :ref:`DPD-SMOOTH <PKG-DPD-SMOOTH>`
   * :ref:`DRUDE <PKG-DRUDE>`
   * :ref:`EFF <PKG-EFF>`
   * :ref:`ELECTRODE <PKG-ELECTRODE>`
   * :ref:`EXTRA-COMPUTE <PKG-EXTRA-COMPUTE>`
   * :ref:`EXTRA-DUMP <PKG-EXTRA-DUMP>`
   * :ref:`EXTRA-FIX <PKG-EXTRA-FIX>`
   * :ref:`EXTRA-MOLECULE <PKG-EXTRA-MOLECULE>`
   * :ref:`EXTRA-PAIR <PKG-EXTRA-PAIR>`
   * :ref:`FEP <PKG-FEP>`
   * :ref:`GPU <PKG-GPU>`
   * :ref:`GRANULAR <PKG-GRANULAR>`
   * :ref:`H5MD <PKG-H5MD>`
   * :ref:`INTEL <PKG-INTEL>`
   * :ref:`INTERLAYER <PKG-INTERLAYER>`
   * :ref:`KIM <PKG-KIM>`
   * :ref:`KOKKOS <PKG-KOKKOS>`
   * :ref:`KSPACE <PKG-KSPACE>`
   * :ref:`LATBOLTZ <PKG-LATBOLTZ>`
   * :ref:`LATTE <PKG-LATTE>`
   * :ref:`LEPTON <PKG-LEPTON>`
   * :ref:`MACHDYN <PKG-MACHDYN>`
   * :ref:`MANIFOLD <PKG-MANIFOLD>`
   * :ref:`MANYBODY <PKG-MANYBODY>`
   * :ref:`MC <PKG-MC>`
   * :ref:`MDI <PKG-MDI>`
   * :ref:`MEAM <PKG-MEAM>`
   * :ref:`MESONT <PKG-MESONT>`
   * :ref:`MGPT <PKG-MGPT>`
   * :ref:`MISC <PKG-MISC>`
   * :ref:`ML-HDNNP <PKG-ML-HDNNP>`
   * :ref:`ML-IAP <PKG-ML-IAP>`
   * :ref:`ML-PACE <PKG-ML-PACE>`
   * :ref:`ML-POD <PKG-ML-POD>`
   * :ref:`ML-QUIP <PKG-ML-QUIP>`
   * :ref:`ML-RANN <PKG-ML-RANN>`
   * :ref:`ML-SNAP <PKG-ML-SNAP>`
   * :ref:`MOFFF <PKG-MOFFF>`
   * :ref:`MOLECULE <PKG-MOLECULE>`
   * :ref:`MOLFILE <PKG-MOLFILE>`
   * :ref:`MPIIO <PKG-MPIIO>`
   * :ref:`MSCG <PKG-MSCG>`
   * :ref:`NETCDF <PKG-NETCDF>`
   * :ref:`OPENMP <PKG-OPENMP>`
   * :ref:`OPT <PKG-OPT>`
   * :ref:`ORIENT <PKG-ORIENT>`
   * :ref:`PERI <PKG-PERI>`
   * :ref:`PHONON <PKG-PHONON>`
   * :ref:`PLUGIN <PKG-PLUGIN>`
   * :ref:`PLUMED <PKG-PLUMED>`
   * :ref:`POEMS <PKG-POEMS>`
   * :ref:`PTM <PKG-PTM>`
   * :ref:`PYTHON <PKG-PYTHON>`
   * :ref:`QEQ <PKG-QEQ>`
   * :ref:`QMMM <PKG-QMMM>`
   * :ref:`QTB <PKG-QTB>`
   * :ref:`REACTION <PKG-REACTION>`
   * :ref:`REAXFF <PKG-REAXFF>`
   * :ref:`REPLICA <PKG-REPLICA>`
   * :ref:`RIGID <PKG-RIGID>`
   * :ref:`SCAFACOS <PKG-SCAFACOS>`
   * :ref:`SHOCK <PKG-SHOCK>`
   * :ref:`SMTBQ <PKG-SMTBQ>`
   * :ref:`SPH <PKG-SPH>`
   * :ref:`SPIN <PKG-SPIN>`
   * :ref:`SRD <PKG-SRD>`
   * :ref:`TALLY <PKG-TALLY>`
   * :ref:`UEF <PKG-UEF>`
   * :ref:`VORONOI <PKG-VORONOI>`
   * :ref:`VTK <PKG-VTK>`
   * :ref:`YAFF <PKG-YAFF>`

----------

.. _PKG-ADIOS:

ADIOS package
------------------

**Contents:**

ADIOS is a high-performance I/O library. This package implements the
:doc:`dump atom/adios <dump_adios>`, :doc:`dump custom/adios <dump_adios>` and
:doc:`read_dump ... format adios <read_dump>`
commands to write and read data using the ADIOS library.

**Authors:** Norbert Podhorszki (ORNL) from the ADIOS developer team.

.. versionadded:: 28Feb2019

**Install:**

This package has :ref:`specific installation instructions <adios>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/ADIOS: filenames -> commands
* src/ADIOS/README
* examples/PACKAGES/adios
* https://github.com/ornladios/ADIOS2
* :doc:`dump atom/adios <dump_adios>`
* :doc:`dump custom/adios <dump_adios>`
* :doc:`read_dump <read_dump>`

----------

.. _PKG-AMOEBA:

AMOEBA package
---------------

**Contents:**

Implementation of the AMOEBA and HIPPO polarized force fields
originally developed by Jay Ponder's group at the U Washington at St
Louis.  The LAMMPS implementation is based on Fortran 90 code
provided by the Ponder group in their
`Tinker MD software <https://dasher.wustl.edu/tinker/>`_.

**Authors:** Josh Rackers and Steve Plimpton (Sandia), Trung Nguyen (U
 Chicago)

**Supporting info:**

* src/AMOEBA: filenames -> commands
* :doc:`AMOEBA and HIPPO howto <Howto_amoeba>`
* :doc:`pair_style amoeba <pair_amoeba>`
* :doc:`pair_style hippo <pair_amoeba>`
* :doc:`atom_style amoeba <atom_style>`
* :doc:`angle_style amoeba <angle_amoeba>`
* :doc:`improper_style amoeba <improper_amoeba>`
* :doc:`fix amoeba/bitorsion <fix_amoeba_bitorsion>`
* :doc:`fix amoeba/pitorsion <fix_amoeba_pitorsion>`
* tools/tinker/tinker2lmp.py
* examples/amoeba

----------

.. _PKG-ASPHERE:

ASPHERE package
---------------

**Contents:**

Computes, time-integration fixes, and pair styles for aspherical
particle models including ellipsoids, 2d lines, and 3d triangles.

**Supporting info:**

* src/ASPHERE: filenames -> commands
* :doc:`Howto spherical <Howto_spherical>`
* :doc:`pair_style gayberne <pair_gayberne>`
* :doc:`pair_style resquared <pair_resquared>`
* :doc:`pair_style ylz <pair_ylz>`
* `doc/PDF/pair_gayberne_extra.pdf <PDF/pair_gayberne_extra.pdf>`_
* `doc/PDF/pair_resquared_extra.pdf <PDF/pair_resquared_extra.pdf>`_
* examples/ASPHERE
* examples/ellipse
* https://www.lammps.org/movies.html#line
* https://www.lammps.org/movies.html#tri

----------

.. _PKG-ATC:

ATC package
----------------

**Contents:**

ATC stands for atoms-to-continuum.  This package implements a
:doc:`fix atc <fix_atc>` command to either couple molecular dynamics
with continuum finite element equations or perform on-the-fly
conversion of atomic information to continuum fields.

**Authors:** Reese Jones, Jeremy Templeton, Jon Zimmerman (Sandia).

**Install:**

This package has :ref:`specific installation instructions <atc>` on the :doc:`Build extras <Build_extras>` page.
The ATC package requires that also the `MANYBODY <PKG-MANYBODY>`_ package is installed.

**Supporting info:**

* src/ATC: filenames -> commands
* src/ATC/README
* :doc:`fix atc <fix_atc>`
* examples/PACKAGES/atc
* https://www.lammps.org/pictures.html#atc

----------

.. _PKG-AWPMD:

AWPMD package
------------------

**Contents:**

AWPMD stands for Antisymmetrized Wave Packet Molecular Dynamics.  This
package implements an atom, pair, and fix style which allows electrons
to be treated as explicit particles in a classical molecular dynamics
model.

**Author:** Ilya Valuev (JIHT, Russia).

**Install:**

This package has :ref:`specific installation instructions <awpmd>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/AWPMD: filenames -> commands
* src/AWPMD/README
* :doc:`pair_style awpmd/cut <pair_awpmd>`
* examples/PACKAGES/awpmd

----------

.. _PKG-BOCS:

BOCS package
-----------------

**Contents:**

This package provides :doc:`fix bocs <fix_bocs>`, a modified version
of :doc:`fix npt <fix_nh>` which includes the pressure correction to
the barostat as outlined in:

N. J. H. Dunn and W. G. Noid, "Bottom-up coarse-grained models that
accurately describe the structure, pressure, and compressibility of
molecular liquids", J. Chem. Phys. 143, 243148 (2015).

**Authors:** Nicholas J. H. Dunn and Michael R. DeLyser (The
Pennsylvania State University)

**Supporting info:**

The BOCS package for LAMMPS is part of the BOCS software package:
`https://github.com/noid-group/BOCS <https://github.com/noid-group/BOCS>`_

See the following reference for information about the entire package:

Dunn, NJH; Lebold, KM; DeLyser, MR; Rudzinski, JF; Noid, WG.
"BOCS: Bottom-Up Open-Source Coarse-Graining Software."
J. Phys. Chem. B. 122, 13, 3363-3377 (2018).

Example inputs are in the examples/PACKAGES/bocs folder.

----------

.. _PKG-BODY:

BODY package
------------

**Contents:**

Body-style particles with internal structure.  Computes,
time-integration fixes, pair styles, as well as the body styles
themselves.  See the :doc:`Howto body <Howto_body>` page for an
overview.

**Supporting info:**

* src/BODY filenames -> commands
* :doc:`Howto_body <Howto_body>`
* :doc:`atom_style body <atom_style>`
* :doc:`fix nve/body <fix_nve_body>`
* :doc:`pair_style body/nparticle <pair_body_nparticle>`
* examples/body

----------

.. _PKG-BPM:

BPM package
------------

**Contents:**

Pair styles, bond styles, fixes, and computes for bonded particle
models for mesoscale simulations of solids and fracture.  See the
:doc:`Howto bpm <Howto_bpm>` page for an overview.

**Authors:** Joel T. Clemmer (Sandia National Labs)

.. versionadded:: 4May2022

**Supporting info:**

* src/BPM filenames -> commands
* :doc:`Howto_bpm <Howto_bpm>`
* :doc:`atom_style bpm/sphere <atom_style>`
* :doc:`bond_style bpm/rotational <bond_bpm_rotational>`
* :doc:`bond_style bpm/spring <bond_bpm_spring>`
* :doc:`compute nbond/atom <compute_nbond_atom>`
* :doc:`fix nve/bpm/sphere <fix_nve_bpm_sphere>`
* :doc:`pair_style bpm/spring <pair_bpm_spring>`
* examples/bpm

----------

.. _PKG-BROWNIAN:

BROWNIAN package
---------------------

**Contents:**

This package provides :doc:`fix brownian, fix brownian/sphere, and
fix brownian/asphere <fix_brownian>` as well as
:doc:`fix propel/self <fix_propel_self>` which allow to do Brownian
Dynamics time integration of point, spherical and aspherical particles
and also support self-propelled particles.

**Authors:** Sam Cameron (University of Bristol),
Stefan Paquay (while at Brandeis University) (initial version of fix propel/self)

.. versionadded:: 14May2021

Example inputs are in the examples/PACKAGES/brownian folder.

----------

.. _PKG-CG-DNA:

CG-DNA package
------------------

**Contents:**

Several pair styles, bond styles, and integration fixes for coarse-grained
modelling of single- and double-stranded DNA and RNA based on the oxDNA and
oxRNA model of Doye, Louis and Ouldridge. The package includes Langevin-type
rigid-body integrators with improved stability.

**Author:** Oliver Henrich (University of Strathclyde, Glasgow).

**Install:**

The CG-DNA package requires that also the `MOLECULE <PKG-MOLECULE>`_ and
`ASPHERE <PKG-ASPHERE>`_ packages are installed.

**Supporting info:**

* src/CG-DNA: filenames -> commands
* /src/CG-DNA/README
* :doc:`pair_style oxdna/\* <pair_oxdna>`
* :doc:`pair_style oxdna2/\* <pair_oxdna2>`
* :doc:`pair_style oxrna2/\* <pair_oxrna2>`
* :doc:`bond_style oxdna/\* <bond_oxdna>`
* :doc:`bond_style oxdna2/\* <bond_oxdna>`
* :doc:`bond_style oxrna2/\* <bond_oxdna>`
* :doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`

----------

.. _PKG-CG-SPICA:

CG-SPICA package
------------------

**Contents:**

Several pair styles and an angle style which implement the
coarse-grained SPICA (formerly called SDK) model which enables
simulation of biological or soft material systems.

**Original Author:** Axel Kohlmeyer (Temple U).

**Maintainers:** Yusuke Miyazaki and Wataru Shinoda (Okayama U).

**Supporting info:**

* src/CG-SPICA: filenames -> commands
* src/CG-SPICA/README
* :doc:`pair_style lj/spica/\* <pair_spica>`
* :doc:`angle_style spica <angle_spica>`
* examples/PACKAGES/cgspica
* https://www.lammps.org/pictures.html#cg
* https://www.spica-ff.org/

----------

.. _PKG-CLASS2:

CLASS2 package
--------------

**Contents:**

Bond, angle, dihedral, improper, and pair styles for the COMPASS
CLASS2 molecular force field.

**Supporting info:**

* src/CLASS2: filenames -> commands
* :doc:`bond_style class2 <bond_class2>`
* :doc:`angle_style class2 <angle_class2>`
* :doc:`dihedral_style class2 <dihedral_class2>`
* :doc:`improper_style class2 <improper_class2>`
* :doc:`pair_style lj/class2 <pair_class2>`

----------

.. _PKG-COLLOID:

COLLOID package
---------------

**Contents:**

Coarse-grained finite-size colloidal particles.  Pair styles and fix
wall styles for colloidal interactions.  Includes the Fast Lubrication
Dynamics (FLD) method for hydrodynamic interactions, which is a
simplified approximation to Stokesian dynamics.

**Authors:** This package includes Fast Lubrication Dynamics pair styles
which were created by Amit Kumar and Michael Bybee from Jonathan
Higdon's group at UIUC.

**Supporting info:**

* src/COLLOID: filenames -> commands
* :doc:`fix wall/colloid <fix_wall>`
* :doc:`pair_style colloid <pair_colloid>`
* :doc:`pair_style yukawa/colloid <pair_yukawa_colloid>`
* :doc:`pair_style brownian <pair_brownian>`
* :doc:`pair_style lubricate <pair_lubricate>`
* :doc:`pair_style lubricateU <pair_lubricateU>`
* examples/colloid
* examples/srd

----------

.. _PKG-COLVARS:

COLVARS package
--------------------

**Contents:**

Colvars stands for collective variables, which can be used to implement
various enhanced sampling methods, including Adaptive Biasing Force,
Metadynamics, Steered MD, Umbrella Sampling and Restraints.  A :doc:`fix
colvars <fix_colvars>` command is implemented which wraps a COLVARS
library, which implements these methods.  simulations.

**Authors:** The COLVARS library is written and maintained by Giacomo
Fiorin (NIH, Bethesda, MD, USA) and Jerome Henin (CNRS, Paris, France),
originally for the NAMD MD code, but with portability in mind.  Axel
Kohlmeyer (Temple U) provided the interface to LAMMPS.

**Install:**

This package has :ref:`specific installation instructions <colvar>` on
the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/COLVARS: filenames -> commands
* `doc/PDF/colvars-refman-lammps.pdf <PDF/colvars-refman-lammps.pdf>`_
* src/COLVARS/README
* lib/colvars/README
* :doc:`fix colvars <fix_colvars>`
* :doc:`group2ndx <group2ndx>`
* :doc:`ndx2group <group2ndx>`
* examples/PACKAGES/colvars

----------

.. _PKG-COMPRESS:

COMPRESS package
----------------

**Contents:**

Compressed output of dump files via the zlib compression library,
using dump styles with a "gz" in their style name.

To use this package you must have the zlib compression library
available on your system.

**Author:** Axel Kohlmeyer (Temple U).

**Install:**

This package has :ref:`specific installation instructions <compress>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/COMPRESS: filenames -> commands
* src/COMPRESS/README
* lib/compress/README
* :doc:`dump atom/gz <dump>`
* :doc:`dump cfg/gz <dump>`
* :doc:`dump custom/gz <dump>`
* :doc:`dump xyz/gz <dump>`

----------

.. _PKG-CORESHELL:

CORESHELL package
-----------------

**Contents:**

Compute and pair styles that implement the adiabatic core/shell model
for polarizability.  The pair styles augment Born, Buckingham, and
Lennard-Jones styles with core/shell capabilities.  The :doc:`compute
temp/cs <compute_temp_cs>` command calculates the temperature of a
system with core/shell particles.  See the :doc:`Howto coreshell
<Howto_coreshell>` page for an overview of how to use this package.

**Author:** Hendrik Heenen (Technical U of Munich).

**Supporting info:**

* src/CORESHELL: filenames -> commands
* :doc:`Howto coreshell <Howto_coreshell>`
* :doc:`Howto polarizable <Howto_polarizable>`
* :doc:`compute temp/cs <compute_temp_cs>`
* :doc:`pair_style born/coul/long/cs <pair_cs>`
* :doc:`pair_style buck/coul/long/cs <pair_cs>`
* :doc:`pair_style lj/cut/coul/long/cs <pair_lj>`
* examples/coreshell

----------

.. _PKG-DIELECTRIC:

DIELECTRIC package
------------------------

**Contents:**

An atom style, multiple pair styles, several fixes, Kspace styles and a
compute for simulating systems using boundary element solvers for
computing the induced charges at the interface between two media with
different dielectric constants.

**Install:**

To use this package, also the :ref:`KSPACE <PKG-KSPACE>` and
:ref:`EXTRA-PAIR <PKG-EXTRA-PAIR>` packages need to be installed.

**Author:** Trung Nguyen and Monica Olvera de la Cruz (Northwestern U)

.. versionadded:: 2Jul2021

**Supporting info:**

* src/DIELECTRIC: filenames -> commands
* :doc:`atom_style dielectric <atom_style>`
* :doc:`pair_style coul/cut/dielectric <pair_dielectric>`
* :doc:`pair_style coul/long/dielectric <pair_dielectric>`
* :doc:`pair_style lj/cut/coul/cut/dielectric <pair_dielectric>`
* :doc:`pair_style lj/cut/coul/debye/dielectric <pair_dielectric>`
* :doc:`pair_style lj/cut/coul/long/dielectric <pair_dielectric>`
* :doc:`pair_style lj/cut/coul/msm/dielectric <pair_dielectric>`
* :doc:`pair_style pppm/dielectric <kspace_style>`
* :doc:`pair_style pppm/disp/dielectric <kspace_style>`
* :doc:`pair_style msm/dielectric <kspace_style>`
* :doc:`fix_style polarize/bem/icc <fix_polarize>`
* :doc:`fix_style polarize/bem/gmres <fix_polarize>`
* :doc:`fix_style polarize/functional <fix_polarize>`
* :doc:`compute efield/atom  <compute_efield_atom>`
* examples/PACKAGES/dielectric

----------

.. _PKG-DIFFRACTION:

DIFFRACTION package
------------------------

**Contents:**

Two computes and a fix for calculating x-ray and electron diffraction
intensities based on kinematic diffraction theory.

**Author:** Shawn Coleman while at the U Arkansas.

**Supporting info:**

* src/DIFFRACTION: filenames -> commands
* :doc:`compute saed <compute_saed>`
* :doc:`compute xrd <compute_xrd>`
* :doc:`fix saed/vtk <fix_saed_vtk>`
* examples/PACKAGES/diffraction

----------

.. _PKG-DIPOLE:

DIPOLE package
--------------

**Contents:**

An atom style and several pair styles for point dipole models with
short-range or long-range interactions.

**Supporting info:**

* src/DIPOLE: filenames -> commands
* :doc:`atom_style dipole <atom_style>`
* :doc:`pair_style lj/cut/dipole/cut <pair_dipole>`
* :doc:`pair_style lj/cut/dipole/long <pair_dipole>`
* :doc:`pair_style lj/long/dipole/long <pair_dipole>`
* :doc:`angle_style dipole <angle_dipole>`
* examples/dipole

----------

.. _PKG-DPD-BASIC:

DPD-BASIC package
--------------------

**Contents:**

Pair styles for the basic dissipative particle dynamics (DPD) method
and DPD thermostatting.

**Author:** Kurt Smith (U Pittsburgh), Martin Svoboda, Martin Lisal (ICPF and UJEP)

**Supporting info:**

* src/DPD-BASIC: filenames -> commands
* :doc:`pair_style dpd <pair_dpd>`
* :doc:`pair_style dpd/tstat <pair_dpd>`
* :doc:`pair_style dpd/ext <pair_dpd_ext>`
* :doc:`pair_style dpd/ext/tstat <pair_dpd_ext>`
* examples/PACKAGES/dpd-basic

----------

.. _PKG-DPD-MESO:

DPD-MESO package
--------------------

**Contents:**

Several extensions of the dissipative particle dynamics (DPD)
method.  Specifically, energy-conserving DPD (eDPD) that can model
non-isothermal processes, many-body DPD (mDPD) for simulating
vapor-liquid coexistence, and transport DPD (tDPD) for modeling
advection-diffusion-reaction systems. The equations of motion of these
DPD extensions are integrated through a modified velocity-Verlet (MVV)
algorithm.

**Author:** Zhen Li (Department of Mechanical Engineering, Clemson University)

**Supporting info:**

* src/DPD-MESO: filenames -> commands
* src/DPD-MESO/README
* :doc:`atom_style edpd <atom_style>`
* :doc:`pair_style edpd <pair_mesodpd>`
* :doc:`pair_style mdpd <pair_mesodpd>`
* :doc:`pair_style tdpd <pair_mesodpd>`
* :doc:`fix mvv/dpd <fix_mvv_dpd>`
* examples/PACKAGES/mesodpd
* https://www.lammps.org/movies.html#mesodpd

----------

.. _PKG-DPD-REACT:

DPD-REACT package
-----------------

**Contents:**

DPD stands for dissipative particle dynamics.  This package implements
coarse-grained DPD-based models for energetic, reactive molecular
crystalline materials.  It includes many pair styles specific to these
systems, including for reactive DPD, where each particle has internal
state for multiple species and a coupled set of chemical reaction ODEs
are integrated each timestep.  Highly accurate time integrators for
isothermal, isoenergetic, isobaric and isenthalpic conditions are
included.  These enable long timesteps via the Shardlow splitting
algorithm.

**Authors:** Jim Larentzos (ARL), Tim Mattox (Engility Corp), and John
Brennan (ARL).

**Supporting info:**

* src/DPD-REACT: filenames -> commands
* /src/DPD-REACT/README
* :doc:`compute dpd <compute_dpd>`
* :doc:`compute dpd/atom <compute_dpd_atom>`
* :doc:`fix eos/cv <fix_eos_table>`
* :doc:`fix eos/table <fix_eos_table>`
* :doc:`fix eos/table/rx <fix_eos_table_rx>`
* :doc:`fix shardlow <fix_shardlow>`
* :doc:`fix rx <fix_rx>`
* :doc:`pair_style table/rx <pair_table_rx>`
* :doc:`pair_style dpd/fdt <pair_dpd_fdt>`
* :doc:`pair_style dpd/fdt/energy <pair_dpd_fdt>`
* :doc:`pair_style exp6/rx <pair_exp6_rx>`
* :doc:`pair_style multi/lucy <pair_multi_lucy>`
* :doc:`pair_style multi/lucy/rx <pair_multi_lucy_rx>`
* examples/PACKAGES/dpd-react

----------

.. _PKG-DPD-SMOOTH:

DPD-SMOOTH package
------------------

**Contents:**

A pair style for smoothed dissipative particle dynamics (SDPD), which
is an extension of smoothed particle hydrodynamics (SPH) to mesoscale
where thermal fluctuations are important (see the
:ref:`SPH package <PKG-SPH>`).
Also two fixes for moving and rigid body integration of SPH/SDPD particles
(particles of atom_style meso).

**Author:** Morteza Jalalvand (Institute for Advanced Studies in Basic
Sciences, Iran).

**Supporting info:**

* src/DPD-SMOOTH: filenames -> commands
* src/DPD-SMOOTH/README
* :doc:`pair_style sdpd/taitwater/isothermal <pair_sdpd_taitwater_isothermal>`
* :doc:`fix meso/move <fix_meso_move>`
* :doc:`fix rigid/meso <fix_rigid_meso>`
* examples/PACKAGES/dpd-smooth

----------

.. _PKG-DRUDE:

DRUDE package
------------------

**Contents:**

Fixes, pair styles, and a compute to simulate thermalized Drude
oscillators as a model of polarization.  See the :doc:`Howto drude <Howto_drude>` and :doc:`Howto drude2 <Howto_drude2>` pages
for an overview of how to use the package.  There are auxiliary tools
for using this package in tools/drude.

**Authors:** Alain Dequidt (U Clermont Auvergne), Julien
Devemy (CNRS), and Agilio Padua (ENS de Lyon).

**Supporting info:**

* src/DRUDE: filenames -> commands
* :doc:`Howto drude <Howto_drude>`
* :doc:`Howto drude2 <Howto_drude2>`
* :doc:`Howto polarizable <Howto_polarizable>`
* src/DRUDE/README
* :doc:`fix drude <fix_drude>`
* :doc:`fix drude/transform/\* <fix_drude_transform>`
* :doc:`compute temp/drude <compute_temp_drude>`
* :doc:`pair_style thole <pair_thole>`
* :doc:`pair_style lj/cut/thole/long <pair_thole>`
* examples/PACKAGES/drude
* tools/drude

----------

.. _PKG-EFF:

EFF package
----------------

**Contents:**

EFF stands for electron force field which allows a classical MD code
to model electrons as particles of variable radius.  This package
contains atom, pair, fix and compute styles which implement the eFF as
described in A. Jaramillo-Botero, J. Su, Q. An, and W.A. Goddard III,
JCC, 2010.  The eFF potential was first introduced by Su and Goddard,
in 2007.  There are auxiliary tools for using this package in
tools/eff; see its README file.

**Author:** Andres Jaramillo-Botero (CalTech).

**Supporting info:**

* src/EFF: filenames -> commands
* src/EFF/README
* :doc:`atom_style electron <atom_style>`
* :doc:`fix nve/eff <fix_nve_eff>`
* :doc:`fix nvt/eff <fix_nh_eff>`
* :doc:`fix npt/eff <fix_nh_eff>`
* :doc:`fix langevin/eff <fix_langevin_eff>`
* :doc:`compute temp/eff <compute_temp_eff>`
* :doc:`pair_style eff/cut <pair_eff>`
* :doc:`pair_style eff/inline <pair_eff>`
* examples/PACKAGES/eff
* tools/eff/README
* tools/eff
* https://www.lammps.org/movies.html#eff

-------------------

.. _PKG-ELECTRODE:

ELECTRODE package
-----------------

**Contents:**

The ELECTRODE package allows the user to enforce a constant potential method for
groups of atoms that interact with the remaining atoms as electrolyte.

**Authors:** The ELECTRODE package is written and maintained by Ludwig
Ahrens-Iwers (TUHH, Hamburg, Germany), Shern Tee (UQ, Brisbane, Australia) and
Robert Meissner (TUHH, Hamburg, Germany).

.. versionadded:: 4May2022

**Install:**

This package has :ref:`specific installation instructions <electrode>` on the
:doc:`Build extras <Build_extras>` page.

**Supporting info:**

* :doc:`fix electrode/conp <fix_electrode>`
* :doc:`fix electrode/conq <fix_electrode>`
* :doc:`fix electrode/thermo <fix_electrode>`

----------

.. _PKG-EXTRA-COMPUTE:

EXTRA-COMPUTE package
---------------------

**Contents:**

Additional compute styles that are less commonly used.

**Supporting info:**

* src/EXTRA-COMPUTE: filenames -> commands
* :doc:`compute <compute>`

----------

.. _PKG-EXTRA-DUMP:

EXTRA-DUMP package
------------------

**Contents:**

Additional dump styles that are less commonly used.

**Supporting info:**

* src/EXTRA-DUMP: filenames -> commands
* :doc:`dump <dump>`

----------

.. _PKG-EXTRA-FIX:

EXTRA-FIX package
-----------------

**Contents:**

Additional fix styles that are less commonly used.

**Supporting info:**

* src/EXTRA-FIX: filenames -> commands
* :doc:`fix <fix>`

----------

.. _PKG-EXTRA-MOLECULE:

EXTRA-MOLECULE package
----------------------

**Contents:**

Additional bond, angle, dihedral, and improper styles that are less commonly used.

**Install:**

To use this package, also the :ref:`MOLECULE <PKG-MOLECULE>` package needs to be installed.

**Supporting info:**

* src/EXTRA-MOLECULE: filenames -> commands
* :doc:`molecular styles <Commands_bond>`

----------

.. _PKG-EXTRA-PAIR:

EXTRA-PAIR package
------------------

**Contents:**

Additional pair styles that are less commonly used.

**Supporting info:**

* src/EXTRA-PAIR: filenames -> commands
* :doc:`pair_style <pair_style>`

----------

.. _PKG-FEP:

FEP package
----------------

**Contents:**

FEP stands for free energy perturbation.  This package provides methods
for performing FEP simulations by using a :doc:`fix adapt/fep
<fix_adapt_fep>` command with soft-core pair potentials, which have a
"soft" in their style name.  There are auxiliary tools for using this
package in tools/fep; see its README file.

**Author:** Agilio Padua (ENS de Lyon)

**Supporting info:**

* src/FEP: filenames -> commands
* src/FEP/README
* :doc:`fix adapt/fep <fix_adapt_fep>`
* :doc:`compute fep <compute_fep>`
* :doc:`pair_style \*/soft <pair_fep_soft>`
* examples/PACKAGES/fep
* tools/fep/README
* tools/fep

----------

.. _PKG-GPU:

GPU package
-----------

**Contents:**

Dozens of pair styles and a version of the PPPM long-range Coulombic
solver optimized for GPUs.  All such styles have a "gpu" as a suffix
in their style name. The GPU code can be compiled with either CUDA or
OpenCL, however the OpenCL variants are no longer actively maintained
and only the CUDA versions are regularly tested.  The
:doc:`Speed_gpu` page gives details of what hardware and GPU
software is required on your system, and details on how to build and
use this package.  Its styles can be invoked at run time via the "-sf
gpu" or "-suffix gpu" :doc:`command-line switches <Run_options>`.  See
also the :ref:`KOKKOS <PKG-KOKKOS>` package, which has GPU-enabled styles.

**Authors:** Mike Brown (Intel) while at Sandia and ORNL and Trung Nguyen
(Northwestern U) while at ORNL and later. AMD HIP support by Evgeny
Kuznetsov, Vladimir Stegailov, and Vsevolod Nikolskiy (HSE University).

**Install:**

This package has :ref:`specific installation instructions <gpu>` on the
:doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/GPU: filenames -> commands
* src/GPU/README
* lib/gpu/README
* :doc:`Accelerator packages <Speed_packages>`
* :doc:`GPU package <Speed_gpu>`
* :doc:`Section 2.6 -sf gpu <Run_options>`
* :doc:`Section 2.6 -pk gpu <Run_options>`
* :doc:`package gpu <package>`
* :doc:`Commands <Commands_all>` pages (:doc:`pair <Commands_pair>`, :doc:`kspace <Commands_kspace>`)
  for styles followed by (g)
* `Benchmarks page <https://www.lammps.org/bench.html>`_ of website

----------

.. _PKG-GRANULAR:

GRANULAR package
----------------

**Contents:**

Pair styles and fixes for finite-size granular particles, which
interact with each other and boundaries via frictional and dissipative
potentials.

**Supporting info:**

* src/GRANULAR: filenames -> commands
* :doc:`Howto granular <Howto_granular>`
* :doc:`fix pour <fix_pour>`
* :doc:`fix wall/gran <fix_wall_gran>`
* :doc:`pair_style gran/hooke <pair_gran>`
* :doc:`pair_style gran/hertz/history <pair_gran>`
* examples/granregion
* examples/pour
* bench/in.chute
* https://www.lammps.org/pictures.html#jamming
* https://www.lammps.org/movies.html#hopper
* https://www.lammps.org/movies.html#dem
* https://www.lammps.org/movies.html#brazil
* https://www.lammps.org/movies.html#granregion

----------

.. _PKG-H5MD:

H5MD package
-----------------

**Contents:**

H5MD stands for HDF5 for MD.  `HDF5 <HDF5_>`_ is a portable, binary,
self-describing file format, used by many scientific simulations.
H5MD is a format for molecular simulations, built on top of HDF5.
This package implements a :doc:`dump h5md <dump_h5md>` command to output
LAMMPS snapshots in this format.

.. _HDF5: https://www.hdfgroup.org/solutions/hdf5

To use this package you must have the HDF5 library available on your
system.

**Author:** Pierre de Buyl (KU Leuven) created both the package and the
H5MD format.

**Install:**

This package has :ref:`specific installation instructions <h5md>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/H5MD: filenames -> commands
* src/H5MD/README
* lib/h5md/README
* :doc:`dump h5md <dump_h5md>`

----------

.. _PKG-INTEL:

INTEL package
------------------

**Contents:**

Dozens of pair, fix, bond, angle, dihedral, improper, and kspace
styles which are optimized for Intel CPUs and KNLs (Knights Landing).
All of them have an "intel" in their style name.  The
:doc:`INTEL package <Speed_intel>` page gives details of what hardware and
compilers are required on your system, and how to build and use this
package.  Its styles can be invoked at run time via the "-sf intel" or
"-suffix intel" :doc:`command-line switches <Run_options>`.  Also see
the :ref:`KOKKOS <PKG-KOKKOS>`, :ref:`OPT <PKG-OPT>`, and :ref:`OPENMP <PKG-OPENMP>` packages,
which have styles optimized for CPUs and KNLs.

You need to have an Intel compiler, version 14 or higher to take full
advantage of this package. While compilation with GNU compilers is
supported, performance will be sub-optimal.

.. note::

   the INTEL package contains styles that require using the
   -restrict flag, when compiling with Intel compilers.

**Author:** Mike Brown (Intel).

**Install:**

This package has :ref:`specific installation instructions <intel>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/INTEL: filenames -> commands
* src/INTEL/README
* :doc:`Accelerator packages <Speed_packages>`
* :doc:`INTEL package <Speed_intel>`
* :doc:`Section 2.6 -sf intel <Run_options>`
* :doc:`Section 2.6 -pk intel <Run_options>`
* :doc:`package intel <package>`
* Search the :doc:`commands <Commands_all>` pages (:doc:`fix <Commands_fix>`, :doc:`compute <Commands_compute>`,
  :doc:`pair <Commands_pair>`, :doc:`bond, angle, dihedral, improper <Commands_bond>`, :doc:`kspace <Commands_kspace>`) for styles followed by (i)
* src/INTEL/TEST
* `Benchmarks page <https://www.lammps.org/bench.html>`_ of website

----------

.. _PKG-INTERLAYER:

INTERLAYER package
------------------

**Contents:**

A collection of pair styles specifically to be used for modeling layered
materials, most commonly graphene sheets (or equivalents).

**Supporting info:**

* src/INTERLAYER: filenames -> commands
* :doc:`Pair style <Commands_pair>` page
* examples/PACKAGES/interlayer

----------

.. _PKG-KIM:

KIM package
-----------

**Contents:**

This package contains a command with a set of sub-commands that serve as a
wrapper on the
`Open Knowledgebase of Interatomic Models (OpenKIM) <https://openkim.org>`_
repository of interatomic models (IMs) enabling compatible ones to be used in
LAMMPS simulations.


This includes :doc:`kim init <kim_commands>`, and
:doc:`kim interactions <kim_commands>` commands to select, initialize and
instantiate the IM, a :doc:`kim query <kim_commands>` command to perform web
queries for material property predictions of OpenKIM IMs, a
:doc:`kim param <kim_commands>` command to access KIM Model Parameters from
LAMMPS, and a :doc:`kim property <kim_commands>` command to write material
properties computed in LAMMPS to standard KIM property instance format.

Support for KIM IMs that conform to the
`KIM Application Programming Interface (API) <https://openkim.org/kim-api/>`_
is provided by the :doc:`pair_style kim <pair_kim>` command.

.. note::

   The command *pair_style kim* is called by *kim interactions* and is not
   recommended to be directly used in input scripts.

To use this package you must have the KIM API library available on your
system. The KIM API is available for download on the
`OpenKIM website <https://openkim.org/kim-api/>`_.
When installing LAMMPS from binary, the kim-api package
is a dependency that is automatically downloaded and installed.

Information about the KIM project can be found at its website:
`https://openkim.org <https://openkim.org>`_.
The KIM project is led by Ellad Tadmor and Ryan Elliott (U Minnesota)
and is funded by the `National Science Foundation <https://www.nsf.gov/>`_.

**Authors:** Ryan Elliott (U Minnesota) is the main developer for the KIM
API and the *pair_style kim* command. Yaser Afshar (U Minnesota),
Axel Kohlmeyer (Temple U), Ellad Tadmor (U Minnesota), and
Daniel Karls (U Minnesota) contributed to the
:doc:`kim command <kim_commands>` interface in close collaboration with
Ryan Elliott.

**Install:**

This package has :ref:`specific installation instructions <kim>` on the
:doc:`Build extras <Build_extras>` page.

**Supporting info:**

* :doc:`kim command <kim_commands>`
* :doc:`pair_style kim <pair_kim>`
* src/KIM: filenames -> commands
* src/KIM/README
* lib/kim/README
* examples/kim

----------

.. _PKG-KOKKOS:

KOKKOS package
--------------

**Contents:**

Dozens of atom, pair, bond, angle, dihedral, improper, fix, compute
styles adapted to compile using the Kokkos library which can convert
them to OpenMP or CUDA code so that they run efficiently on multicore
CPUs, KNLs, or GPUs.  All the styles have a "kk" as a suffix in their
style name.  The :doc:`KOKKOS package <Speed_kokkos>` page gives
details of what hardware and software is required on your system, and
how to build and use this package.  Its styles can be invoked at run
time via the "-sf kk" or "-suffix kk" :doc:`command-line switches <Run_options>`.  Also see the :ref:`GPU <PKG-GPU>`, :ref:`OPT <PKG-OPT>`,
:ref:`INTEL <PKG-INTEL>`, and :ref:`OPENMP <PKG-OPENMP>` packages, which
have styles optimized for CPUs, KNLs, and GPUs.

You must have a C++14 compatible compiler to use this package.
KOKKOS makes extensive use of advanced C++ features, which can
expose compiler bugs, especially when compiling for maximum
performance at high optimization levels. Please see the file
lib/kokkos/README for a list of compilers and their respective
platforms, that are known to work.

**Authors:** The KOKKOS package was created primarily by Christian Trott
and Stan Moore (Sandia), with contributions from other folks as well.
It uses the open-source `Kokkos library <https://github.com/kokkos>`_
which was developed by Carter Edwards, Christian Trott, and others at
Sandia, and which is included in the LAMMPS distribution in
lib/kokkos.

**Install:**

This package has :ref:`specific installation instructions <kokkos>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/KOKKOS: filenames -> commands
* src/KOKKOS/README
* lib/kokkos/README
* :doc:`Accelerator packages <Speed_packages>`
* :doc:`KOKKOS package <Speed_kokkos>`
* :doc:`Section 2.6 -k on ... <Run_options>`
* :doc:`Section 2.6 -sf kk <Run_options>`
* :doc:`Section 2.6 -pk kokkos <Run_options>`
* :doc:`package kokkos <package>`
* Search the :doc:`commands <Commands_all>` pages (:doc:`fix <Commands_fix>`, :doc:`compute <Commands_compute>`,
  :doc:`pair <Commands_pair>`, :doc:`bond, angle, dihedral, improper <Commands_bond>`,
  :doc:`kspace <Commands_kspace>`) for styles followed by (k)
* `Benchmarks page <https://www.lammps.org/bench.html>`_ of website

----------

.. _PKG-KSPACE:

KSPACE package
--------------

**Contents:**

A variety of long-range Coulombic solvers, as well as pair styles
which compute the corresponding short-range pairwise Coulombic
interactions.  These include Ewald, particle-particle particle-mesh
(PPPM), and multilevel summation method (MSM) solvers.

**Install:**

Building with this package requires a 1d FFT library be present on
your system for use by the PPPM solvers.  This can be the KISS FFT
library provided with LAMMPS, third party libraries like FFTW, or a
vendor-supplied FFT library.  See the :doc:`Build settings <Build_settings>` page for details on how to select
different FFT options for your LAMPMS build.

**Supporting info:**

* src/KSPACE: filenames -> commands
* :doc:`kspace_style <kspace_style>`
* `doc/PDF/kspace.pdf <PDF/kspace.pdf>`_
* :doc:`Howto tip3p <Howto_tip3p>`
* :doc:`Howto tip4p <Howto_tip4p>`
* :doc:`Howto spc <Howto_spc>`
* :doc:`pair_style coul <pair_coul>`
* Search the :doc:`pair style <Commands_pair>` page for styles with "long" or "msm" in name
* examples/peptide
* bench/in.rhodo

----------

.. _PKG-LATBOLTZ:

LATBOLTZ package
----------------

**Contents:**

Fixes which implement a background Lattice-Boltzmann (LB) fluid, which
can be used to model MD particles influenced by hydrodynamic forces.

**Authors:** Frances Mackay and Colin Denniston (University of Western
Ontario).

**Install:**

The LATBOLTZ package requires that LAMMPS is build in :ref:`MPI parallel mode <serial>`.

**Supporting info:**

* src/LATBOLTZ: filenames -> commands
* src/LATBOLTZ/README
* :doc:`fix lb/fluid <fix_lb_fluid>`
* :doc:`fix lb/momentum <fix_lb_momentum>`
* :doc:`fix lb/viscous <fix_lb_viscous>`
* examples/PACKAGES/latboltz

----------

.. _PKG-LATTE:

LATTE package
-------------

**Contents:**

A fix command which wraps the LATTE DFTB code, so that molecular
dynamics can be run with LAMMPS using density-functional tight-binding
quantum forces calculated by LATTE.

More information on LATTE can be found at this website:
`https://github.com/lanl/LATTE <latte-home_>`_.  A brief technical
description is given with the :doc:`fix latte <fix_latte>` command.

.. _latte-home: https://github.com/lanl/LATTE

**Authors:** Christian Negre (LANL) and Steve Plimpton (Sandia).  LATTE
itself is developed at Los Alamos National Laboratory by Marc
Cawkwell, Anders Niklasson, and Christian Negre.

**Install:**

This package has :ref:`specific installation instructions <latte>` on
the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/LATTE: filenames -> commands
* src/LATTE/README
* lib/latte/README
* :doc:`fix latte <fix_latte>`
* examples/latte
* `LAMMPS-LATTE tutorial <https://github.com/lanl/LATTE/wiki/Using-LATTE-through-LAMMPS>`_

----------

.. _PKG-LEPTON:

LEPTON package
--------------

**Contents:**

Styles for pair, bond, and angle forces that evaluate the potential
function from a string using the `Lepton mathematical expression parser
<https://simtk.org/projects/lepton>`_.  Lepton is a C++ library that is
bundled with `OpenMM <https://openmm.org/>`_ and can be used for
parsing, evaluating, differentiating, and analyzing mathematical
expressions.  This is a more lightweight and efficient alternative for
evaluating custom potential function to an embedded Python interpreter
as used in the :ref:`PYTHON package <PKG-PYTHON>`.  On the other hand,
since the potentials are evaluated form analytical expressions, they are
more precise than what can be done with :ref:`tabulated potentials
<tabulate>`.

**Authors:** Axel Kohlmeyer (Temple U).  Lepton itself is developed
by Peter Eastman at Stanford University.

.. versionadded:: TBD

**Install:**

This package has :ref:`specific installation instructions <lepton>` on
the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/LEPTON: filenames -> commands
* lib/lepton/README.md
* :doc:`pair_style lepton <pair_lepton>`
* :doc:`bond_style lepton <bond_lepton>`
* :doc:`angle_style lepton <angle_lepton>`
* :doc:`dihedral_style lepton <dihedral_lepton>`

----------

.. _PKG-MACHDYN:

MACHDYN package
----------------

**Contents:**

An atom style, fixes, computes, and several pair styles which
implements smoothed Mach dynamics (SMD) for solids, which is a model
related to smoothed particle hydrodynamics (SPH) for liquids (see the
:ref:`SPH package <PKG-SPH>`).

This package solves solids mechanics problems via a state of the art
stabilized meshless method with hourglass control.  It can specify
hydrostatic interactions independently from material strength models,
i.e. pressure and deviatoric stresses are separated.  It provides many
material models (Johnson-Cook, plasticity with hardening,
Mie-Grueneisen, Polynomial EOS) and allows new material models to be
added.  It implements rigid boundary conditions (walls) which can be
specified as surface geometries from \*.STL files.

**Author:** Georg Ganzenmuller (Fraunhofer-Institute for High-Speed
Dynamics, Ernst Mach Institute, Germany).

**Install:**

This package has :ref:`specific installation instructions <machdyn>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/MACHDYN: filenames -> commands
* src/MACHDYN/README
* `doc/PDF/MACHDYN_LAMMPS_userguide.pdf <PDF/MACHDYN_LAMMPS_userguide.pdf>`_
* examples/PACKAGES/machdyn
* https://www.lammps.org/movies.html#smd

----------

.. _PKG-MANIFOLD:

MANIFOLD package
---------------------

**Contents:**

Several fixes and a "manifold" class which enable simulations of
particles constrained to a manifold (a 2D surface within the 3D
simulation box).  This is done by applying the RATTLE constraint
algorithm to formulate single-particle constraint functions
g(xi,yi,zi) = 0 and their derivative (i.e. the normal of the manifold)
n = grad(g).

**Author:** Stefan Paquay (until 2017: Eindhoven University of
Technology (TU/e), The Netherlands; since 2017: Brandeis University,
Waltham, MA, USA)

**Supporting info:**

* src/MANIFOLD: filenames -> commands
* src/MANIFOLD/README
* :doc:`Howto manifold <Howto_manifold>`
* :doc:`fix manifoldforce <fix_manifoldforce>`
* :doc:`fix nve/manifold/rattle <fix_nve_manifold_rattle>`
* :doc:`fix nvt/manifold/rattle <fix_nvt_manifold_rattle>`
* examples/PACKAGES/manifold
* https://www.lammps.org/movies.html#manifold

----------

.. _PKG-MANYBODY:

MANYBODY package
----------------

**Contents:**

A variety of many-body and bond-order potentials.  These include
(AI)REBO, BOP, EAM, EIM, Stillinger-Weber, and Tersoff potentials.

**Supporting info:**

* src/MANYBODY: filenames -> commands
* :doc:`Pair style <Commands_pair>` page
* examples/comb
* examples/eim
* examples/nb3d
* examples/shear
* examples/streitz
* examples/vashishta
* bench/in.eam

----------

.. _PKG-MC:

MC package
----------

**Contents:**

Several fixes and a pair style that have Monte Carlo (MC) or MC-like
attributes.  These include fixes for creating, breaking, and swapping
bonds, for performing atomic swaps, and performing grand canonical
MC (GCMC), semi-grand canonical MC (SGCMC), or similar processes in
conjunction with molecular dynamics (MD).

**Supporting info:**

* src/MC: filenames -> commands
* :doc:`fix atom/swap <fix_atom_swap>`
* :doc:`fix bond/break <fix_bond_break>`
* :doc:`fix bond/create <fix_bond_create>`
* :doc:`fix bond/create/angle <fix_bond_create>`
* :doc:`fix bond/swap <fix_bond_swap>`
* :doc:`fix charge/regulation <fix_charge_regulation>`
* :doc:`fix gcmc <fix_gcmc>`
* :doc:`fix sgcmc <fix_sgcmc>`
* :doc:`fix tfmc <fix_tfmc>`
* :doc:`fix widom <fix_widom>`
* :doc:`pair_style dsmc <pair_dsmc>`
* https://www.lammps.org/movies.html#gcmc

----------

.. _PKG-MDI:

MDI package
----------------

**Contents:**

A LAMMPS command and fixes to allow client-server coupling of LAMMPS
to other atomic or molecular simulation codes or materials modeling
workflows via the `MolSSI Driver Interface
(MDI) library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`_.

**Author:** Taylor Barnes - MolSSI, taylor.a.barnes at gmail.com

.. versionadded:: 14May2021

**Install:**

This package has :ref:`specific installation instructions <mdi>` on
the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/MDI/README
* lib/mdi/README
* :doc:`Howto MDI <Howto_mdi>`
* :doc:`mdi <mdi>`
* :doc:`fix mdi/qm <fix_mdi_qm>`
* examples/PACKAGES/mdi

----------

.. _PKG-MEAM:

MEAM package
------------------

**Contents:**

A pair style for the modified embedded atom (MEAM) potential
translated from the Fortran version in the (obsolete) MEAM package
to plain C++. The MEAM fully replaces the MEAM package, which
has been removed from LAMMPS after the 12 December 2018 version.

**Author:** Sebastian Huetter, (Otto-von-Guericke University Magdeburg)
based on the Fortran version of Greg Wagner (Northwestern U) while at
Sandia.

**Supporting info:**

* src/MEAM: filenames -> commands
* src/MEAM/README
* :doc:`pair_style meam <pair_meam>`
* examples/meam

----------

.. _PKG-MESONT:

MESONT package
-------------------

**Contents:**

MESONT is a LAMMPS package for simulation of nanomechanics of nanotubes
(NTs). The model is based on a coarse-grained representation of NTs as
"flexible cylinders" consisting of a variable number of
segments. Internal interactions within a NT and the van der Waals
interaction between the tubes are described by a mesoscopic force field
designed and parameterized based on the results of atomic-level
molecular dynamics simulations. The description of the force field is
provided in the papers listed below.

This package contains two independent implementations of this model:
:doc:`pair_style mesont/tpm <pair_mesont_tpm>` is the original
implementation of the model based on a Fortran library in the
``lib/mesont`` folder. The second implementation is provided by the
mesocnt styles (:doc:`bond_style mesocnt <bond_mesocnt>`,
:doc:`angle_style mesocnt <angle_mesocnt>` and :doc:`pair_style mesocnt
<pair_mesocnt>`).  The mesocnt implementation has the same features as
the original implementation with the addition of friction, but is
directly implemented in C++, interfaces more cleanly with general LAMMPS
functionality, and is typically faster. It also does not require its own
atom style and can be installed without any external libraries.

**Download of potential files:**

The potential files for these pair styles are *very* large and thus are
not included in the regular downloaded packages of LAMMPS or the git
repositories.  Instead, they will be automatically downloaded from a web
server when the package is installed for the first time.

**Authors of the *mesont* styles:**

Maxim V. Shugaev (University of Virginia), Alexey N. Volkov (University
of Alabama), Leonid V. Zhigilei (University of Virginia)

**Author of the *mesocnt* styles:**
Philipp Kloza (U Cambridge)

.. versionadded:: 15Jun2020

**Supporting info:**

* src/MESONT: filenames -> commands
* src/MESONT/README
* :doc:`atom_style mesont <atom_style>`
* :doc:`pair_style mesont/tpm <pair_mesont_tpm>`
* :doc:`compute mesont <compute_mesont>`
* :doc:`bond_style mesocnt <bond_mesocnt>`
* :doc:`angle_style mesocnt <angle_mesocnt>`
* :doc:`pair_style mesocnt <pair_mesocnt>`
* examples/PACKAGES/mesont
* tools/mesont

----------

.. _PKG-MGPT:

MGPT package
-----------------

**Contents:**

A pair style which provides a fast implementation of the quantum-based
MGPT multi-ion potentials.  The MGPT or model GPT method derives from
first-principles DFT-based generalized pseudopotential theory (GPT)
through a series of systematic approximations valid for mid-period
transition metals with nearly half-filled d bands.  The MGPT method
was originally developed by John Moriarty at LLNL.  The pair style in
this package calculates forces and energies using an optimized
matrix-MGPT algorithm due to Tomas Oppelstrup at LLNL.

**Authors:** Tomas Oppelstrup and John Moriarty (LLNL).

**Supporting info:**

* src/MGPT: filenames -> commands
* src/MGPT/README
* :doc:`pair_style mgpt <pair_mgpt>`
* examples/PACKAGES/mgpt

----------

.. _PKG-MISC:

MISC package
------------

**Contents:**

A variety of compute, fix, pair, bond styles with specialized
capabilities that don't align with other packages.  Do a directory
listing, "ls src/MISC", to see the list of commands.

.. note::

   the MISC package contains styles that require using the
   -restrict flag, when compiling with Intel compilers.

**Supporting info:**

* src/MISC: filenames -> commands
* :doc:`bond_style special <bond_special>`
* :doc:`compute viscosity/cos <compute_viscosity_cos>`
* :doc:`fix accelerate/cos <fix_accelerate_cos>`
* :doc:`fix imd <fix_imd>`
* :doc:`fix ipi <fix_ipi>`
* :doc:`pair_style agni <pair_agni>`
* :doc:`pair_style list <pair_list>`
* :doc:`pair_style srp <pair_srp>`
* :doc:`pair_style tracker <pair_tracker>`

----------

.. _PKG-ML-HDNNP:

ML-HDNNP package
------------------

**Contents:**

A :doc:`pair_style hdnnp <pair_hdnnp>` command which allows to use
high-dimensional neural network potentials (HDNNPs), a form of machine learning
potentials. HDNNPs must be carefully trained prior to their application in a
molecular dynamics simulation.

.. _n2p2: https://github.com/CompPhysVienna/n2p2

To use this package you must have the `n2p2 <n2p2_>`_ library installed and
compiled on your system.

**Author:** Andreas Singraber

.. versionadded:: 27May2021

**Install:**

This package has :ref:`specific installation instructions <ml-hdnnp>` on the
:doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/ML-HDNNP: filenames -> commands
* src/ML-HDNNP/README
* lib/hdnnp/README
* :doc:`pair_style hdnnp <pair_hdnnp>`
* examples/PACKAGES/hdnnp

----------

.. _PKG-ML-IAP:

ML-IAP package
--------------

**Contents:**

A general interface for machine-learning interatomic potentials, including PyTorch.

**Install:**

To use this package, also the :ref:`ML-SNAP <PKG-ML-SNAP>` package needs
to be installed.  To make the *mliappy* model available, also the
:ref:`PYTHON <PKG-PYTHON>` package needs to be installed, the version
of Python must be 3.6 or later, and the `cython <https://cython.org/>`_ software
must be installed.

**Author:** Aidan Thompson (Sandia), Nicholas Lubbers (LANL).

.. versionadded:: 30Jun2020

**Supporting info:**

* src/ML-IAP: filenames -> commands
* src/ML-IAP/README.md
* :doc:`pair_style mliap <pair_mliap>`
* :doc:`compute_style mliap <compute_mliap>`
* examples/mliap (see README)

When built with the *mliappy* model this package includes an extension for
coupling with Python models, including PyTorch. In this case, the Python
interpreter linked to LAMMPS will need the ``cython`` and ``numpy`` modules
installed.  The provided examples build models with PyTorch, which would
therefore also needs to be installed to run those examples.

----------

.. _PKG-ML-PACE:

ML-PACE package
-------------------

**Contents:**

A pair style for the Atomic Cluster Expansion potential (ACE).
ACE is a methodology for deriving a highly accurate classical potential
fit to a large archive of quantum mechanical (DFT) data. The ML-PACE
package provides an efficient implementation for running simulations
with ACE potentials.

**Authors:**

This package was written by Yury Lysogorskiy^1,
Cas van der Oord^2, Anton Bochkarev^1,
Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1.

 ^1: Ruhr-University Bochum, Bochum, Germany

 ^2: University of Cambridge, Cambridge, United Kingdom

 ^3: Sandia National Laboratories, Albuquerque, New Mexico, USA

 ^4: University of British Columbia, Vancouver, BC, Canada

.. versionadded:: 14May2021

**Install:**

This package has :ref:`specific installation instructions <ml-pace>` on the
:doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/ML-PACE: filenames -> commands
* :doc:`pair_style pace <pair_pace>`
* examples/PACKAGES/pace

----------

.. _PKG-ML-POD:

ML-POD package
-------------------

**Contents:**

A pair style and fitpod style for Proper Orthogonal Descriptors
(POD). POD is a methodology for deriving descriptors based on the proper
orthogonal decomposition. The ML-POD package provides an efficient
implementation for running simulations with POD potentials, along with
fitting the potentials natively in LAMMPS.

**Authors:**

Ngoc Cuong Nguyen (MIT), Andrew Rohskopf (Sandia)

.. versionadded:: 22Dec2022

**Install:**

This package has :ref:`specific installation instructions <ml-pod>` on the
:doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/ML-POD: filenames -> commands
* :doc:`pair_style pod <pair_pod>`
* :doc:`command_style fitpod <fitpod_command>`
* examples/PACKAGES/pod

----------

.. _PKG-ML-QUIP:

ML-QUIP package
-----------------

**Contents:**

A :doc:`pair_style quip <pair_quip>` command which wraps the `QUIP
libAtoms library <quip_>`_, which includes a variety of interatomic
potentials, including Gaussian Approximation Potential (GAP) models
developed by the Cambridge University group.

.. _quip: https://github.com/libAtoms/QUIP

To use this package you must have the QUIP libAtoms library available
on your system.

**Author:** Albert Bartok (Cambridge University)

**Install:**

This package has :ref:`specific installation instructions <ml-quip>` on the
:doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/ML-QUIP: filenames -> commands
* src/ML-QUIP/README
* :doc:`pair_style quip <pair_quip>`
* examples/PACKAGES/quip

----------

.. _PKG-ML-RANN:

ML-RANN package
-----------------

**Contents:**

A pair style for using rapid atomistic neural network (RANN) potentials.
These neural network potentials work by first generating a series of symmetry
functions from the neighbor list and then using these values as the input layer
of a neural network.

**Authors:**

This package was written by Christopher Barrett
with contributions by Doyl Dickel, Mississippi State University.

.. versionadded:: 27May2021

**Supporting info:**

* src/ML-RANN: filenames -> commands
* :doc:`pair_style rann <pair_rann>`
* examples/PACKAGES/rann

----------

.. _PKG-ML-SNAP:

ML-SNAP package
---------------

**Contents:**

A pair style for the spectral neighbor analysis potential (SNAP).
SNAP is methodology for deriving a highly accurate classical potential
fit to a large archive of quantum mechanical (DFT) data. Also several
computes which analyze attributes of the potential.

**Author:** Aidan Thompson (Sandia).

**Supporting info:**

* src/ML-SNAP: filenames -> commands
* :doc:`pair_style snap <pair_snap>`
* :doc:`compute sna/atom <compute_sna_atom>`
* :doc:`compute sna/grid <compute_sna_atom>`
* :doc:`compute sna/grid/local <compute_sna_atom>`
* :doc:`compute snad/atom <compute_sna_atom>`
* :doc:`compute snav/atom <compute_sna_atom>`
* examples/snap

----------

.. _PKG-MOFFF:

MOFFF package
------------------

**Contents:**

Pair, angle and improper styles needed to employ the MOF-FF
force field by Schmid and coworkers with LAMMPS.
MOF-FF is a first principles derived force field with the primary aim
to simulate MOFs and related porous framework materials, using spherical
Gaussian charges. It is described in S. Bureekaew et al., Phys. Stat. Sol. B
2013, 250, 1128-1141.
For the usage of MOF-FF see the example in the example directory as
well as the `MOF+ <MOFplus_>`_ website.

.. _MOFplus: https://www.mofplus.org/content/show/MOF-FF

**Author:** Hendrik Heenen (Technical U of Munich),
Rochus Schmid (Ruhr-University Bochum).

**Supporting info:**

* src/MOFFF: filenames -> commands
* src/MOFFF/README
* :doc:`pair_style buck6d/coul/gauss <pair_buck6d_coul_gauss>`
* :doc:`angle_style class2 <angle_class2>`
* :doc:`angle_style cosine/buck6d <angle_cosine_buck6d>`
* :doc:`improper_style inversion/harmonic <improper_inversion_harmonic>`
* examples/PACKAGES/mofff

----------

.. _PKG-MOLECULE:

MOLECULE package
----------------

**Contents:**

A large number of atom, pair, bond, angle, dihedral, improper styles
that are used to model molecular systems with fixed covalent bonds.
The pair styles include the Dreiding (hydrogen-bonding) and CHARMM
force fields, and a TIP4P water model.

**Supporting info:**

* src/MOLECULE: filenames -> commands
* :doc:`atom_style <atom_style>`
* :doc:`bond_style <bond_style>`
* :doc:`angle_style <angle_style>`
* :doc:`dihedral_style <dihedral_style>`
* :doc:`improper_style <improper_style>`
* :doc:`pair_style hbond/dreiding/lj <pair_hbond_dreiding>`
* :doc:`pair_style lj/charmm/coul/charmm <pair_charmm>`
* :doc:`Howto bioFF <Howto_bioFF>`
* examples/cmap
* examples/dreiding
* examples/micelle,
* examples/peptide
* bench/in.chain
* bench/in.rhodo

----------

.. _PKG-MOLFILE:

MOLFILE package
--------------------

**Contents:**

A :doc:`dump molfile <dump_molfile>` command which uses molfile plugins
that are bundled with the `VMD <vmd-home_>`_
molecular visualization and analysis program, to enable LAMMPS to dump
snapshots in formats compatible with various molecular simulation
tools.

To use this package you must have the desired VMD plugins available on
your system.

Note that this package only provides the interface code, not the
plugins themselves, which will be accessed when requesting a specific
plugin via the :doc:`dump molfile <dump_molfile>` command.  Plugins can
be obtained from a VMD installation which has to match the platform
that you are using to compile LAMMPS for. By adding plugins to VMD,
support for new file formats can be added to LAMMPS (or VMD or other
programs that use them) without having to re-compile the application
itself.  More information about the VMD molfile plugins can be found
at
`https://www.ks.uiuc.edu/Research/vmd/plugins/molfile <https://www.ks.uiuc.edu/Research/vmd/plugins/molfile>`_.

**Author:** Axel Kohlmeyer (Temple U).

**Install:**

This package has :ref:`specific installation instructions <molfile>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/MOLFILE: filenames -> commands
* src/MOLFILE/README
* lib/molfile/README
* :doc:`dump molfile <dump_molfile>`

----------

.. _PKG-MPIIO:

MPIIO package
-------------

**Contents:**

Support for parallel output/input of dump and restart files via the
MPIIO library.  It adds :doc:`dump styles <dump>` with a "mpiio" in
their style name.  Restart files with an ".mpiio" suffix are also
written and read in parallel.

.. warning::

   The MPIIO package is currently unmaintained and has become
   unreliable. Use with caution.


**Install:**

The MPIIO package requires that LAMMPS is build in :ref:`MPI parallel mode <serial>`.

**Supporting info:**

* src/MPIIO: filenames -> commands
* :doc:`dump <dump>`
* :doc:`restart <restart>`
* :doc:`write_restart <write_restart>`
* :doc:`read_restart <read_restart>`

----------

.. _PKG-MSCG:

MSCG package
------------

**Contents:**

A :doc:`fix mscg <fix_mscg>` command which can parameterize a
Multi-Scale Coarse-Graining (MSCG) model using the open-source `MS-CG library <mscg-home_>`_.

.. _mscg-home: https://github.com/uchicago-voth/MSCG-release

To use this package you must have the MS-CG library available on your
system.

**Authors:** The fix was written by Lauren Abbott (Sandia).  The MS-CG
library was developed by Jacob Wagner in Greg Voth's group at the
University of Chicago.

**Install:**

This package has :ref:`specific installation instructions <mscg>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/MSCG: filenames -> commands
* src/MSCG/README
* lib/mscg/README
* examples/mscg

----------

.. _PKG-NETCDF:

NETCDF package
-------------------

**Contents:**

Dump styles for writing NetCDF formatted dump files.  NetCDF is a
portable, binary, self-describing file format developed on top of
HDF5. The file contents follow the AMBER NetCDF trajectory conventions
(https://ambermd.org/netcdf/nctraj.xhtml), but include extensions.

To use this package you must have the NetCDF library available on your
system.

Note that NetCDF files can be directly visualized with the following
tools:

* `Ovito <ovito_>`_ (Ovito supports the AMBER convention and the extensions mentioned above)
* `VMD <vmd-home_>`_

.. _ovito: https://www.ovito.org

.. _vmd-home: https://www.ks.uiuc.edu/Research/vmd/

**Author:** Lars Pastewka (Karlsruhe Institute of Technology).

**Install:**

This package has :ref:`specific installation instructions <netcdf>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/NETCDF: filenames -> commands
* src/NETCDF/README
* lib/netcdf/README
* :doc:`dump netcdf <dump_netcdf>`

----------

.. _PKG-OPENMP:

OPENMP package
----------------

**Contents:**

Hundreds of pair, fix, compute, bond, angle, dihedral, improper, and
kspace styles which are altered to enable threading on many-core CPUs
via OpenMP directives.  All of them have an "omp" in their style name.
The :doc:`OPENMP package <Speed_omp>` page gives details of what hardware
and compilers are required on your system, and how to build and use
this package.  Its styles can be invoked at run time via the "-sf omp"
or "-suffix omp" :doc:`command-line switches <Run_options>`.  Also see
the :ref:`KOKKOS <PKG-KOKKOS>`, :ref:`OPT <PKG-OPT>`, and :ref:`INTEL <PKG-INTEL>`
packages, which have styles optimized for CPUs.

**Author:** Axel Kohlmeyer (Temple U).

.. note::

   To enable multi-threading support the compile flag "-fopenmp"
   and the link flag "-fopenmp" (for GNU compilers, you have to look up
   the equivalent flags for other compilers) must be used to build LAMMPS.
   When using Intel compilers, also the "-restrict" flag is required.
   The OPENMP package can be compiled without enabling OpenMP; then
   all code will be compiled as serial and the only improvement over the
   regular styles are some data access optimization. These flags should
   be added to the CCFLAGS and LINKFLAGS lines of your Makefile.machine.
   See src/MAKE/OPTIONS/Makefile.omp for an example.

Once you have an appropriate Makefile.machine, you can
install/un-install the package and build LAMMPS in the usual manner:

**Install:**

This package has :ref:`specific installation instructions <openmp>` on
the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/OPENMP: filenames -> commands
* src/OPENMP/README
* :doc:`Accelerator packages <Speed_packages>`
* :doc:`OPENMP package <Speed_omp>`
* :doc:`Command line option -suffix/-sf omp <Run_options>`
* :doc:`Command line option -package/-pk omp <Run_options>`
* :doc:`package omp <package>`
* Search the :doc:`commands <Commands_all>` pages (:doc:`fix <Commands_fix>`, :doc:`compute <Commands_compute>`,
  :doc:`pair <Commands_pair>`, :doc:`bond, angle, dihedral, improper <Commands_bond>`,
  :doc:`kspace <Commands_kspace>`) for styles followed by (o)
* `Benchmarks page <https://www.lammps.org/bench.html>`_ of website

----------

.. _PKG-OPT:

OPT package
-----------

**Contents:**

A handful of pair styles which are optimized for improved CPU
performance on single or multiple cores.  These include EAM, LJ,
CHARMM, and Morse potentials.  The styles have an "opt" suffix in
their style name.  The :doc:`OPT package <Speed_opt>` page gives
details of how to build and use this package.  Its styles can be
invoked at run time via the "-sf opt" or "-suffix opt" :doc:`command-line switches <Run_options>`.  See also the :ref:`KOKKOS <PKG-KOKKOS>`,
:ref:`INTEL <PKG-INTEL>`, and :ref:`OPENMP <PKG-OPENMP>` packages, which
have styles optimized for CPU performance.

**Authors:** James Fischer (High Performance Technologies), David Richie,
and Vincent Natoli (Stone Ridge Technology).

**Install:**

This package has :ref:`specific installation instructions <opt>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/OPT: filenames -> commands
* :doc:`Accelerator packages <Speed_packages>`
* :doc:`OPT package <Speed_opt>`
* :doc:`Section 2.6 -sf opt <Run_options>`
* Search the :doc:`pair style <Commands_pair>` page for styles followed by (t)
* `Benchmarks page <https://www.lammps.org/bench.html>`_ of website

.. _PKG-ORIENT:

ORIENT package
--------------

**Contents:**

A few fixes that apply orientation dependent forces for studying
grain boundary migration.

**Supporting info:**

* src/ORIENT: filenames -> commands
* :doc:`fix orient/bcc <fix_orient>`
* :doc:`fix orient/fcc <fix_orient>`
* :doc:`fix orient/eco <fix_orient_eco>`

----------

.. _PKG-PERI:

PERI package
------------

**Contents:**

An atom style, several pair styles which implement different
Peridynamics materials models, and several computes which calculate
diagnostics.  Peridynamics is a particle-based meshless continuum
model.

**Authors:** The original package was created by Mike Parks (Sandia).
Additional Peridynamics models were added by Rezwanur Rahman and John
Foster (UTSA).

**Supporting info:**

* src/PERI: filenames -> commands
* :doc:`Peridynamics Howto <Howto_peri>`
* `doc/PDF/PDLammps_overview.pdf <PDF/PDLammps_overview.pdf>`_
* `doc/PDF/PDLammps_EPS.pdf <PDF/PDLammps_EPS.pdf>`_
* `doc/PDF/PDLammps_VES.pdf <PDF/PDLammps_VES.pdf>`_
* :doc:`atom_style peri <atom_style>`
* :doc:`pair_style peri/\* <pair_peri>`
* :doc:`compute damage/atom <compute_damage_atom>`
* :doc:`compute plasticity/atom <compute_plasticity_atom>`
* examples/peri
* https://www.lammps.org/movies.html#peri

----------

.. _PKG-PHONON:

PHONON package
-------------------

**Contents:**

A :doc:`fix phonon <fix_phonon>` command that calculates dynamical
matrices, which can then be used to compute phonon dispersion
relations, directly from molecular dynamics simulations.
And a :doc:`dynamical_matrix <dynamical_matrix>` as well as a
:doc:`third_order <third_order>` command to compute the dynamical matrix
and third order tensor from finite differences.

**Install:**

The PHONON package requires that also the `KSPACE <PKG-KSPACE>`_
package is installed.


**Authors:** Ling-Ti Kong (Shanghai Jiao Tong University) for "fix phonon"
and Charlie Sievers (UC Davis) for "dynamical_matrix" and "third_order"

**Supporting info:**

* src/PHONON: filenames -> commands
* src/PHONON/README
* :doc:`fix phonon <fix_phonon>`
* :doc:`dynamical_matrix <dynamical_matrix>`
* :doc:`third_order <third_order>`
* examples/PACKAGES/phonon

----------

.. _PKG-PLUGIN:

PLUGIN package
--------------

**Contents:**

A :doc:`plugin <plugin>` command that can load and unload several
kind of styles in LAMMPS from shared object files at runtime without
having to recompile and relink LAMMPS.

When the environment variable ``LAMMPS_PLUGIN_PATH`` is set, then LAMMPS
will search the directory (or directories) listed in this path for files
with names that end in ``plugin.so`` (e.g. ``helloplugin.so``) and will
try to load the contained plugins automatically at start-up.

**Authors:** Axel Kohlmeyer (Temple U)

.. versionadded:: 8Apr2021

**Supporting info:**

* src/PLUGIN: filenames -> commands
* :doc:`plugin command <plugin>`
* :doc:`Information on writing plugins <Developer_plugins>`
* examples/plugin

----------

.. _PKG-PLUMED:

PLUMED package
-------------------

**Contents:**

The fix plumed command allows you to use the PLUMED free energy plugin
for molecular dynamics to analyze and bias your LAMMPS trajectory on
the fly.  The PLUMED library is called from within the LAMMPS input
script by using the :doc:`fix plumed <fix_plumed>` command.

**Authors:** The `PLUMED library <https://www.plumed.org>`_ is written
and maintained by Massimilliano Bonomi, Giovanni Bussi, Carlo Camiloni,
and Gareth Tribello.

**Install:**

This package has :ref:`specific installation instructions <plumed>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/PLUMED/README
* lib/plumed/README
* :doc:`fix plumed <fix_plumed>`
* examples/PACKAGES/plumed

----------

.. _PKG-POEMS:

POEMS package
-------------

**Contents:**

A fix that wraps the Parallelizable Open source Efficient Multibody
Software (POEMS) library, which is able to simulate the dynamics of
articulated body systems.  These are systems with multiple rigid
bodies (collections of particles) whose motion is coupled by
connections at hinge points.

**Author:** Rudra Mukherjee (JPL) while at RPI.

**Install:**

This package has :ref:`specific installation instructions <poems>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/POEMS: filenames -> commands
* src/POEMS/README
* lib/poems/README
* :doc:`fix poems <fix_poems>`
* examples/rigid

----------

.. _PKG-PTM:

PTM package
----------------

**Contents:**

A :doc:`compute ptm/atom <compute_ptm_atom>` command that calculates
local structure characterization using the Polyhedral Template
Matching methodology.

**Author:** Peter Mahler Larsen (MIT).

**Supporting info:**

* src/PTM: filenames not starting with ptm\_ -> commands
* src/PTM: filenames starting with ptm\_ -> supporting code
* src/PTM/LICENSE
* :doc:`compute ptm/atom <compute_ptm_atom>`

----------

.. _PKG-PYTHON:

PYTHON package
--------------

**Contents:**

A :doc:`python <python>` command which allow you to execute Python code
from a LAMMPS input script.  The code can be in a separate file or
embedded in the input script itself.  See the :doc:`Python call
<Python_call>` page for an overview of using Python from LAMMPS in this
manner and all the :doc:`Python <Python_head>` manual pages for other
ways to use LAMMPS and Python together.

.. note::

   Building with the PYTHON package assumes you have a Python development
   environment (headers and libraries) available on your system, which needs
   to be either Python version 2.7 or Python 3.5 and later.

**Install:**

This package has :ref:`specific installation instructions <python>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/PYTHON: filenames -> commands
* :doc:`Python call <Python_head>`
* lib/python/README
* examples/python

----------

.. _PKG-QEQ:

QEQ package
-----------

**Contents:**

Several fixes for performing charge equilibration (QEq) via different
algorithms.  These can be used with pair styles that perform QEq as
part of their formulation.

**Supporting info:**

* src/QEQ: filenames -> commands
* :doc:`fix qeq/\* <fix_qeq>`
* examples/qeq
* examples/streitz

----------

.. _PKG-QMMM:

QMMM package
-----------------

**Contents:**

A :doc:`fix qmmm <fix_qmmm>` command which allows LAMMPS to be used as
the MM code in a QM/MM simulation.  This is currently only available
in combination with the `Quantum ESPRESSO <espresso_>`_ package.

.. _espresso: https://www.quantum-espresso.org

To use this package you must have Quantum ESPRESSO (QE) available on
your system and include its coupling library in the compilation and
then compile LAMMPS as a library.  For QM/MM calculations you then
build a custom binary with MPI support, that sets up 3 partitions with
MPI sub-communicators (for inter- and intra-partition communication)
and then calls the corresponding library interfaces on each partition
(2x LAMMPS and 1x QE).

The current implementation supports an ONIOM style mechanical coupling
and a multi-pole based electrostatic coupling to the Quantum ESPRESSO
plane wave DFT package.  The QM/MM interface has been written in a
manner that coupling to other QM codes should be possible without
changes to LAMMPS itself.

**Authors:** Axel Kohlmeyer (Temple U). Mariella Ippolito and Carlo Cavazzoni (CINECA, Italy)

**Install:**

This package has :ref:`specific installation instructions <qmmm>`
on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/QMMM: filenames -> commands
* src/QMMM/README
* lib/qmmm/README
* :doc:`fix phonon <fix_phonon>`
* lib/qmmm/example-ec/README
* lib/qmmm/example-mc/README

----------

.. _PKG-QTB:

QTB package
----------------

**Contents:**

Two fixes which provide a self-consistent quantum treatment of
vibrational modes in a classical molecular dynamics simulation.  By
coupling the MD simulation to a colored thermostat, it introduces zero
point energy into the system, altering the energy power spectrum and
the heat capacity to account for their quantum nature. This is useful
when modeling systems at temperatures lower than their classical
limits or when temperatures ramp across the classical limits in a
simulation.

**Author:** Yuan Shen (Stanford U).

**Supporting info:**

* src/QTB: filenames -> commands
* src/QTB/README
* :doc:`fix qtb <fix_qtb>`
* :doc:`fix qbmsst <fix_qbmsst>`
* examples/PACKAGES/qtb

----------

.. _PKG-REACTION:

REACTION package
---------------------

**Contents:**

This package implements the REACTER protocol, which allows for complex
bond topology changes (reactions) during a running MD simulation when
using classical force fields. Topology changes are defined in pre- and
post-reaction molecule templates and can include creation and deletion
of bonds, angles, dihedrals, impropers, atom types, bond types, angle
types, dihedral types, improper types, and/or atomic charges. Other
options currently available include reaction constraints (e.g., angle
and Arrhenius constraints), deletion of reaction byproducts or other
small molecules, creation of new atoms or molecules bonded to existing
atoms, and using LAMMPS variables for input parameters.

**Author:** Jacob R. Gissinger (NASA Langley Research Center).

**Supporting info:**

* src/REACTION: filenames -> commands
* src/REACTION/README
* :doc:`fix bond/react <fix_bond_react>`
* examples/PACKAGES/reaction
* `2017 LAMMPS Workshop <https://www.lammps.org/workshops/Aug17/pdf/gissinger.pdf>`_
* `2019 LAMMPS Workshop <https://www.lammps.org/workshops/Aug19/talk_gissinger.pdf>`_
* `2021 LAMMPS Workshop <https://www.lammps.org/workshops/Aug21/talk/jacob-gissinger/>`_
* `REACTER website (reacter.org) <https://www.reacter.org/>`_

----------

.. _PKG-REAXFF:

REAXFF package
------------------

**Contents:**

A pair style which implements the ReaxFF potential in C/C++.  ReaxFF
is a universal reactive force field.  See the src/REAXFF/README file
for more info on differences between the two packages.  Also two fixes
for monitoring molecules as bonds are created and destroyed.

**Author:** Hasan Metin Aktulga (MSU) while at Purdue University.

**Supporting info:**

* src/REAXFF: filenames -> commands
* src/REAXFF/README
* :doc:`pair_style reaxff <pair_reaxff>`
* :doc:`fix reaxff/bonds <fix_reaxff_bonds>`
* :doc:`fix reaxff/species <fix_reaxff_species>`
* examples/reaxff

----------

.. _PKG-REPLICA:

REPLICA package
---------------

**Contents:**

A collection of multi-replica methods which can be used when running
multiple LAMMPS simulations (replicas).  See the :doc:`Howto replica <Howto_replica>` page for an overview of how to run
multi-replica simulations in LAMMPS.  Methods in the package include
nudged elastic band (NEB), parallel replica dynamics (PRD),
temperature accelerated dynamics (TAD), parallel tempering, and a
verlet/split algorithm for performing long-range Coulombics on one set
of processors, and the remainder of the force field calculation on
another set.

**Supporting info:**

* src/REPLICA: filenames -> commands
* :doc:`Howto replica <Howto_replica>`
* :doc:`neb <neb>`
* :doc:`prd <prd>`
* :doc:`tad <tad>`
* :doc:`temper <temper>`,
* :doc:`temper/npt <temper_npt>`,
* :doc:`temper/grem <temper_grem>`,
* :doc:`run_style verlet/split <run_style>`
* examples/neb
* examples/prd
* examples/tad
* examples/PACKAGES/grem

----------

.. _PKG-RIGID:

RIGID package
-------------

**Contents:**

Fixes which enforce rigid constraints on collections of atoms or
particles.  This includes SHAKE and RATTLE, as well as various
rigid-body integrators for a few large bodies or many small bodies.
Also several computes which calculate properties of rigid bodies.

**Supporting info:**

* src/RIGID: filenames -> commands
* :doc:`compute erotate/rigid <compute_erotate_rigid>`
* :doc:`fix shake <fix_shake>`
* :doc:`fix rattle <fix_shake>`
* :doc:`fix rigid/\* <fix_rigid>`
* examples/ASPHERE
* examples/rigid
* bench/in.rhodo
* https://www.lammps.org/movies.html#box
* https://www.lammps.org/movies.html#star

----------

.. _PKG-SCAFACOS:

SCAFACOS package
---------------------

**Contents:**

A KSpace style which wraps the `ScaFaCoS Coulomb solver library <http://www.scafacos.de>`_ to compute long-range Coulombic
interactions.

To use this package you must have the ScaFaCoS library available on
your system.

**Author:** Rene Halver (JSC) wrote the scafacos LAMMPS command.

ScaFaCoS itself was developed by a consortium of German research
facilities with a BMBF (German Ministry of Science and Education)
funded project in 2009-2012. Participants of the consortium were the
Universities of Bonn, Chemnitz, Stuttgart, and Wuppertal as well as
the Forschungszentrum Juelich.

**Install:**

This package has :ref:`specific installation instructions <scafacos>` on the :doc:`Build extras <Build_extras>` page.
The SCAFACOS package requires that LAMMPS is build in :ref:`MPI parallel mode <serial>`.

**Supporting info:**

* src/SCAFACOS: filenames -> commands
* src/SCAFACOS/README
* :doc:`kspace_style scafacos <kspace_style>`
* :doc:`kspace_modify <kspace_modify>`
* examples/PACKAGES/scafacos

----------

.. _PKG-SHOCK:

SHOCK package
-------------

**Contents:**

Fixes for running impact simulations where a shock-wave passes through
a material.

**Supporting info:**

* src/SHOCK: filenames -> commands
* :doc:`fix append/atoms <fix_append_atoms>`
* :doc:`fix msst <fix_msst>`
* :doc:`fix nphug <fix_nphug>`
* :doc:`fix wall/piston <fix_wall_piston>`
* examples/hugoniostat
* examples/msst

----------

.. _PKG-SMTBQ:

SMTBQ package
------------------

**Contents:**

Pair styles which implement Second Moment Tight Binding models.
One with QEq charge equilibration (SMTBQ) for the description of
ionocovalent bonds in oxides, and two more as plain SMATB models.

**Authors:** SMTBQ: Nicolas Salles, Emile Maras, Olivier Politano, and Robert
Tetot (LAAS-CNRS, France);
SMATB: Daniele Rapetti (Politecnico di Torino)

**Supporting info:**

* src/SMTBQ: filenames -> commands
* src/SMTBQ/README
* :doc:`pair_style smtbq <pair_smtbq>`
* :doc:`pair_style smatb <pair_smatb>`, :doc:`pair_style smatb/single <pair_smatb>`
* examples/PACKAGES/smtbq

----------

.. _PKG-SPH:

SPH package
----------------

**Contents:**

An atom style, fixes, computes, and several pair styles which
implements smoothed particle hydrodynamics (SPH) for liquids.  See the
related :ref:`MACHDYN package <PKG-MACHDYN>` package for smooth Mach dynamics
(SMD) for solids.

This package contains ideal gas, Lennard-Jones equation of states,
Tait, and full support for complete (i.e. internal-energy dependent)
equations of state.  It allows for plain or Monaghans XSPH integration
of the equations of motion.  It has options for density continuity or
density summation to propagate the density field.  It has
:doc:`set <set>` command options to set the internal energy and density
of particles from the input script and allows the same quantities to
be output with thermodynamic output or to dump files via the :doc:`compute property/atom <compute_property_atom>` command.

**Author:** Georg Ganzenmuller (Fraunhofer-Institute for High-Speed
Dynamics, Ernst Mach Institute, Germany).

**Supporting info:**

* src/SPH: filenames -> commands
* src/SPH/README
* `doc/PDF/SPH_LAMMPS_userguide.pdf <PDF/SPH_LAMMPS_userguide.pdf>`_
* examples/PACKAGES/sph
* https://www.lammps.org/movies.html#sph

----------

.. _PKG-SPIN:

SPIN package
------------

**Contents:**

Model atomic magnetic spins classically, coupled to atoms moving in
the usual manner via MD.  Various pair, fix, and compute styles.

**Author:** Julien Tranchida (Sandia).

**Supporting info:**

* src/SPIN: filenames -> commands
* :doc:`Howto spins <Howto_spins>`
* :doc:`pair_style spin/dipole/cut <pair_spin_dipole>`
* :doc:`pair_style spin/dipole/long <pair_spin_dipole>`
* :doc:`pair_style spin/dmi <pair_spin_dmi>`
* :doc:`pair_style spin/exchange <pair_spin_exchange>`
* :doc:`pair_style spin/exchange/biquadratic <pair_spin_exchange>`
* :doc:`pair_style spin/magelec <pair_spin_magelec>`
* :doc:`pair_style spin/neel <pair_spin_neel>`
* :doc:`fix nve/spin <fix_nve_spin>`
* :doc:`fix langevin/spin <fix_langevin_spin>`
* :doc:`fix precession/spin <fix_precession_spin>`
* :doc:`compute spin <compute_spin>`
* :doc:`neb/spin <neb_spin>`
* examples/SPIN

----------

.. _PKG-SRD:

SRD package
-----------

**Contents:**

A pair of fixes which implement the Stochastic Rotation Dynamics (SRD)
method for coarse-graining of a solvent, typically around large
colloidal particles.

**Supporting info:**

* src/SRD: filenames -> commands
* :doc:`fix srd <fix_srd>`
* :doc:`fix wall/srd <fix_wall_srd>`
* examples/srd
* examples/ASPHERE
* https://www.lammps.org/movies.html#tri
* https://www.lammps.org/movies.html#line
* https://www.lammps.org/movies.html#poly

----------

.. _PKG-TALLY:

TALLY package
------------------

**Contents:**

Several compute styles that can be called when pairwise interactions
are calculated to tally information (forces, heat flux, energy,
stress, etc) about individual interactions.

**Author:** Axel Kohlmeyer (Temple U).

**Supporting info:**

* src/TALLY: filenames -> commands
* src/TALLY/README
* :doc:`compute \*/tally <compute_tally>`
* examples/PACKAGES/tally

----------

.. _PKG-UEF:

UEF package
----------------

**Contents:**

A fix style for the integration of the equations of motion under
extensional flow with proper boundary conditions, as well as several
supporting compute styles and an output option.

**Author:** David Nicholson (MIT).

**Supporting info:**

* src/UEF: filenames -> commands
* src/UEF/README
* :doc:`fix nvt/uef <fix_nh_uef>`
* :doc:`fix npt/uef <fix_nh_uef>`
* :doc:`compute pressure/uef <compute_pressure_uef>`
* :doc:`compute temp/uef <compute_temp_uef>`
* :doc:`dump cfg/uef <dump_cfg_uef>`
* examples/uef

----------

.. _PKG-VORONOI:

VORONOI package
---------------

**Contents:**

A compute command which calculates the Voronoi tesselation of a
collection of atoms by wrapping the `Voro++ library <voro-home_>`_.  This
can be used to calculate the local volume or each atoms or its near
neighbors.

.. _voro-home: https://math.lbl.gov/voro++

To use this package you must have the Voro++ library available on your
system.

**Author:** Daniel Schwen (INL) while at LANL.  The open-source Voro++
library was written by Chris Rycroft (Harvard U) while at UC Berkeley
and LBNL.

**Install:**

This package has :ref:`specific installation instructions <voronoi>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/VORONOI: filenames -> commands
* src/VORONOI/README
* lib/voronoi/README
* :doc:`compute voronoi/atom <compute_voronoi_atom>`
* examples/voronoi

----------

.. _PKG-VTK:

VTK package
----------------

**Contents:**

A :doc:`dump vtk <dump_vtk>` command which outputs snapshot info in the
`VTK format <vtk_>`_, enabling visualization by `Paraview <paraview_>`_ or
other visualization packages.

.. _vtk: https://www.vtk.org

.. _paraview: https://www.paraview.org

To use this package you must have VTK library available on your
system.

**Authors:** Richard Berger (JKU) and Daniel Queteschiner (DCS Computing).

**Install:**

This package has :ref:`specific installation instructions <vtk>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/VTK: filenames -> commands
* src/VTK/README
* lib/vtk/README
* :doc:`dump vtk <dump_vtk>`

----------

.. _PKG-YAFF:

YAFF package
-----------------

**Contents:**

Some potentials that are also implemented in the Yet Another Force Field (`YAFF <yaff_>`_) code.
The expressions and their use are discussed in the following papers

* Vanduyfhuys et al., J. Comput. Chem., 36 (13), 1015-1027 (2015) `link <vanduyfhuys2015_>`_
* Vanduyfhuys et al., J. Comput. Chem., 39 (16), 999-1011 (2018) `link <vanduyfhuys2018_>`_

which discuss the `QuickFF <quickff_>`_ methodology.

.. _vanduyfhuys2015: https://doi.org/10.1002/jcc.23877
.. _vanduyfhuys2018: https://doi.org/10.1002/jcc.25173
.. _quickff: https://molmod.github.io/QuickFF
.. _yaff: https://github.com/molmod/yaff

**Author:** Steven Vandenbrande.

.. versionadded:: 1Feb2019

**Supporting info:**

* src/YAFF/README
* :doc:`angle_style cross <angle_cross>`
* :doc:`angle_style mm3 <angle_mm3>`
* :doc:`bond_style mm3 <bond_mm3>`
* :doc:`improper_style distharm <improper_distharm>`
* :doc:`improper_style sqdistharm <improper_sqdistharm>`
* :doc:`pair_style mm3/switch3/coulgauss/long <pair_lj_switch3_coulgauss_long>`
* :doc:`pair_style lj/switch3/coulgauss/long <pair_lj_switch3_coulgauss_long>`
* examples/PACKAGES/yaff
