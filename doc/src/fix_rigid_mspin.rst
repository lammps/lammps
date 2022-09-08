.. index:: fix rigid/nvt/mspin

fix rigid/nvt/mspin command
===========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rigid/nvt/mspin bodystyle keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* rigid/nvt/mspin = style name of this fix
* bodystyle = molecule or group
* one or more keyword/value pairs may be appended
* keyword = *temp* or *bfield* or *dpcut* or *alpha* or *beta*

  .. parsed-literal::

       *temp* values = Tstart Tstop Tdamp
         Tstart,Tstop = desired temperature at start/stop of run (temperature units)
         Tdamp = temperature damping parameter (time units)
       *bfield* values = Bx By Bz uniform
         Bx,By,Bz = external magnetic field vector components (Tesla)
         uniform = optional, if magnetic external field is uniform
       *dpcut* value = cutoff
         cutoff = a dipolar interaction cutoff (distance units)
       *alpha* value = Alpha
         Alpha = scaling factor for dipolar interactions
       *beta* value = Beta
         Beta = scaling factor for magnetic dipole moment

Examples
""""""""

.. code-block:: LAMMPS

    fix 1 all rigid/nvt/mspin molecule temp 300 300 100 dpcut 64
    fix 1 all rigid/nvt/mspin molecule temp 300 300 100 bfield 5.0 0.0 0.0
    fix 1 all rigid/nvt/mspin molecule temp 300 300 100 bfield 0.0 0.0 0.3 uniform dpcut 64 beta 8.0

Description
"""""""""""
LAMMPS package for atomistic molecular dynamics simulation of
magnetic nanoparticles. This fix calculates magnetic dipole moment
of a nanoparticle core at each timestep and update force and
energy of the particle given an external magnetic field or
due to dipolar interactions between multiple nanoparticles.

The package assumes a net magnetic dipole moment vector embedded in
the nanoparticle core. To determine the magnitude and direction of
the dipole moment :math:`\vec{\mu}`, two individual atoms of the
nanoparticle core are set with non-zero "magnetic charge" values,
:math:`|q_m|`. A detailed description of the method can be found in
(:ref:`Mahmood2022 <Mahmood2022>`).

Prior to this fix, a new atom property using the :doc:`fix property/atom <fix_property_atom>`
command should be defined to allocate memory for the :math:`q_m` values.

.. parsed-literal::

  fix   qm    all   property/atom   d_qm

The magnetic charge value can be computed using

.. math::

  q_m = \frac{\mu}{d}

where :math:`\mu` is the dipole moment of magnetic nanoparticle in real units,
and :math:`d` is the distance between the two chosen atoms from the
nanoparticle core. If the chosen atoms are located at the particle surface,
:math:`d` will be same as the particle diameter.

This package currently supports only the LAMMPS :doc:`real units <units>`.
Coversion of the dipole moment from SI (i.e., :math:`\text{A m}^{-2}`) to *real* units
can be done using the following equation.

.. math::

  \mu_{\text{real}} = \frac{\mu_{\text{SI}}}{1.6022 \times 10^{-21}}

In current implementation, each NP can carry only 1 magnetic dipole moment,
specified by a pair of atoms carrying non-zero :math:`q_m` values.
The :math:`q_m` values of the selected atoms can be set using the
:doc:`set command <set>` and their atomic indices. The magnitude
of the two atoms should be same while the signs should be opposite.
For example,

.. parsed-literal::

  set     atom     4    d_qm  0.0353
  set     atom     5    d_qm -0.0353

Depending on the interatomic distance between the two atoms and
the specified :math:`q_m` values, the package calculates the embedded
dipole moment of a particle from the :math:`-q_m` atom towards the
:math:`+q_m` atom during simulation and evaluates the dynamics
of the particle utilizing the :ref:`RIGID <PKG-RIGID>` package.
The dynamics of the solvent and ligands attached to the surface of
the nanoparticle are evaluated using :doc:`fix NVT <fix_nh>` simulation.

.. note::

  The force field parameters of the selected atoms should be kept unchanged,
  since the atoms are used only to identify the position and direction
  of the dipole moment. The magnetic charge values do not affect the
  interatomic interactions.

The keyword *bfield* is used to specify an external magnetic field
vector acting on the particle, with the field vector components towards the
x, y, z directions should follow the keyword.

The optional *uniform* keyword can be specified after the external field arguments.
The default field type is assumed to be non-uniform.

The *dpcut* keyword can be used to set a cutoff value for magnetic
dipolar interaction between multiple magnetic particles. If this keyword
is not specified, magnetic dipolar interaction is turned off.

Two scaling factors *alpha* and *beta* can be used to scale the dipolar
interaction and the magnitude of the dipole moment respectively. The
default values of the scaling factors are set to 1.0.

Restrictions
""""""""""""

This fix is part of the :ref:`MSPIN <PKG-MSPIN>` package. It is only enabled if
LAMMPS was built with that package.

The :ref:`MSPIN <PKG-MSPIN>` package requires
LAMMPS be built with the :ref:`RIGID <PKG-RIGID>` package.
See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rigid <fix_rigid>`,
:doc:`compute mspin/distance <compute_mspin_distance>`,
:doc:`compute mspin/energy <compute_mspin_energy>`

----------

.. _Mahmood2022:

**(Mahmood)** A.U. Mahmood and Y.G. Yingling,
*All-Atom Simulation Method for Zeeman Alignment and Dipolar Assembly
of Magnetic Nanoparticles*, `Journal of Chemical Theory and Computation 18 [5],
3122-3135 (2022) <https://doi.org/10.1021/acs.jctc.1c01253>`_
