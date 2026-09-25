.. index:: fix tumble

fix tumble command
==================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID tumble mode rate seed keyword

* ID, group-ID are documented in :doc:`fix <fix>` command
* tumble = style name of this fix command
* mode = *dipole*
* rate = average number of tumbles per unit time (inverse time units)
* seed = random number seed to use for the tumbles (positive integer)
* zero or one keyword may be appended
* keyword = *planar_rotation*

  .. parsed-literal::

       *planar_rotation* value = none
         confine new orientations to the xy plane in 3d simulations

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all tumble dipole 0.5 12345
   fix 1 active tumble dipole 2.0 4321 planar_rotation

Description
"""""""""""

.. versionadded:: TBD

Reorient the particles in the group at random times, as in the
run-and-tumble model of the motion of swimming bacteria such as *E.
coli* :ref:`(Berg1972) <Berg1972>`, :ref:`(Tailleur2008) <Tailleur2008>`.
Combined with a self-propulsion force along the orientation of the
particles from :doc:`fix propel/self <fix_propel_self>`, this results in
"run-and-tumble particles" which move in straight runs at constant
speed, interrupted by tumbles that select a new direction.  See the
:doc:`Howto active matter <Howto_active>` page for an overview of active
matter models in LAMMPS.

The tumbles of each particle form a Poisson process with the given
*rate*: on every timestep a particle in the group tumbles with the
probability :math:`1 - \exp(-\lambda\,\Delta t)`, where :math:`\lambda`
is the *rate* and :math:`\Delta t` the size of the timestep.  Hence the
average time between tumbles (the "run time") is :math:`1/\lambda`,
independent of the size of the timestep.  When a particle tumbles, its
orientation vector is replaced by a new direction that is uniformly
distributed on the unit circle (for 2d simulations) or on the unit
sphere (for 3d simulations); the length of the orientation vector is
preserved.  Between tumbles the orientation is not changed by this fix.
Other fixes may still rotate it, for instance :doc:`fix brownian/sphere
<fix_brownian>` (rotational diffusion) or :doc:`fix nve/sphere
<fix_nve_sphere>` with the *update dipole* keyword (torques).

Since every tumble completely decorrelates the orientation, the
orientation autocorrelation function of non-interacting run-and-tumble
particles decays as :math:`\langle \mathbf{e}(t) \cdot \mathbf{e}(0)
\rangle = \exp(-\lambda t)`.  For comparison, for active Brownian
particles with the rotational diffusion coefficient :math:`D_r` it
decays as :math:`\exp(-(d-1) D_r t)` in :math:`d` dimensions, so that
run-and-tumble particles with :math:`\lambda = (d-1) D_r` have the same
long-time diffusive behavior :ref:`(Cates2013) <Cates2013>`.

For mode *dipole*, the orientation of a particle is its dipole vector
:math:`\mathbf{\mu}_i` as defined by :doc:`atom_style dipole
<atom_style>`.  All atoms in the group must have a dipole moment of
non-zero length, which can be set with the :doc:`set <set>` command
(e.g. with the *dipole/random* keyword).  A tumble keeps the length of
the dipole vector and only changes its direction.

The reorientation is applied after the positions have been updated by
the time integrator on a timestep and before the forces are computed, so
that a self-propulsion force computed from the orientation uses the new
direction on the same timestep.

The optional *planar_rotation* keyword confines the new orientations to
the xy plane in 3d simulations, i.e. the z component of the orientation
is set to zero.  This corresponds to the *planar_rotation* keyword of
:doc:`fix brownian/sphere <fix_brownian>` and should be used together
with it.  In 2d simulations the new orientation is always in the xy
plane and this keyword is not allowed.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of its random number generator to :doc:`binary
restart files <restart>`, so that a simulation continued from a restart
reproduces the same stochastic trajectory (provided the number of MPI
processes is unchanged).  No global or per-atom quantities are stored by
this fix for access by various :doc:`output commands <Howto_output>`.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

When using :doc:`run_style respa <run_style>`, the tumbles are performed
once per outermost timestep.

Restrictions
""""""""""""

This fix is part of the BROWNIAN package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Mode *dipole* requires that atoms store a dipole moment as defined by
the :doc:`atom_style dipole <atom_style>` command, which is part of the
DIPOLE package.

Related commands
""""""""""""""""

:doc:`fix propel/self <fix_propel_self>`, :doc:`fix brownian/sphere <fix_brownian>`,
:doc:`fix propel/ou <fix_propel_ou>`, :doc:`fix align/neighbor <fix_align_neighbor>`,
:doc:`set <set>`

Default
"""""""

none

----------

.. _Berg1972:

**(Berg1972)** H. C. Berg and D. A. Brown, Chemotaxis in Escherichia coli analysed by three-dimensional tracking, Nature 239, 500 (1972).

.. _Tailleur2008:

**(Tailleur2008)** J. Tailleur and M. E. Cates, Statistical Mechanics of Interacting Run-and-Tumble Bacteria, Phys. Rev. Lett. 100, 218103 (2008).

.. _Cates2013:

**(Cates2013)** M. E. Cates and J. Tailleur, When are active Brownian particles and run-and-tumble particles equivalent? Consequences for motility-induced phase separation, EPL 101, 20010 (2013).
