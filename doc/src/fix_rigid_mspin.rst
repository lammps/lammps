.. index:: fix rigid/nvt/mspin

fix rigid/nvt/mspin command
===========================

Examples
""""""""

.. code-block:: LAMMPS

    fix 1 all rigid/nvt/mspin molecule temp 300 300 10 dpcut 64
    fix 1 all rigid/nvt/mspin molecule temp 300 300 10 &
            bfield 5.0 0.0 0.0
    fix 1 all rigid/nvt/mspin molecule temp 300 300 10 &
            bfield 0.0 0.0 0.3 uniform dpcut 64 beta 8.0

Description
"""""""""""
LAMMPS packages for atomistic molecular dynamics simulations of
magnetic nanoparticles. This fix calculates magnetic dipole moment
of a rigid nanoparticle core at each timestep and update force and
energy of the particle given an external magnetic field or
due to magnetic dipolar interactions between multiple particles.

To determine the magnitude and direction of the dipole moment vector,
:math:`\mathbf{m}` two individual atoms of a nanoparticle are set with
non-zero and equal "magnetic charge" value, :math:`\abs{q_m}`.

The force field parameters of the atoms should be kept unchanged.
A detailed description of this method can be found in (:ref:`Mahmood2022 <Mahmood2022>`).

The keyword *bfield* is used to specify an external magnetic field
vector acting on the particle, with the vector components towards the
x, y, z directions should follow the keyword.

The *uniform* keyword can be specified after the external field arguments.
The default field type is assumed to be non-uniform.

The *dpcut* keyword can be used to set a cutoff value for magnetic
dipolar interaction between multiple magnetic particles. If this keyword
is not specified, magnetic dipolar interaction is turned off.

Two scaling factors *alpha* and *beta* can be used to adjust the dipolar
interaction and the magnitude of the dipole moment respectively. The
default values of the scaling factors are set to 1.0.

Restrictions
""""""""""""

This fix is part of the MSPIN package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` page for more info.

This fix requires LAMMPS be built with the RIGID package.
See the :doc:`Build settings <Build_settings>` page for details.

----------

.. _Mahmood2022:

**(Mahmood)** A.U. Mahmood and Y.G. Yingling,
*All-Atom Simulation Method for Zeeman Alignment and Dipolar Assembly
of Magnetic Nanoparticles*, `Journal of Chemical Theory and Computation 18 [5], 3122&ndash;3135 (2022) <https://doi.org/10.1021/acs.jctc.1c01253>`_
