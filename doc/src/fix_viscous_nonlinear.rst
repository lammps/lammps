.. index:: fix viscous/nonlinear

fix viscous/nonlinear command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID viscous/nonlinear rho_fluid mu_fluid keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* viscous/nonlinear = style name of this fix command
* rho_fluid = mass density of the surrounding fluid (mass/volume units)
* mu_fluid = dynamic viscosity of the surrounding fluid (pressure\*time units)
* zero or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *velocity*
       *velocity* values = Vx Vy Vz
         Vx,Vy,Vz = components of the (uniform) fluid velocity (velocity units)

Examples
""""""""

.. code-block:: LAMMPS

   fix drag all viscous/nonlinear 1.2 1.8e-5
   fix drag flow viscous/nonlinear 1.2 1.8e-5 velocity 0.0 0.0 0.4

Description
"""""""""""

.. versionadded:: 4Jul2026

Add a nonlinear (Reynolds-number dependent) drag force to each
finite-size spherical particle in the group, modeling the interaction
with a uniform background fluid (e.g. an upward gas stream).  Unlike
:doc:`fix viscous <fix_viscous>`, which applies a drag force strictly
proportional to the particle velocity (Stokes drag), this fix uses the
standard drag-coefficient relation with the Schiller-Naumann
correlation, which is accurate over a much wider range of particle
Reynolds numbers.

The drag force on particle *i* is

.. math::

   \vec{F}_i = -\frac{1}{2}\, C_d\, \rho_f\, \pi r_i^2\, |\vec{v}_{rel}|\, \vec{v}_{rel}

where :math:`r_i` is the particle radius, :math:`\rho_f` is the fluid
mass density, :math:`\vec{v}_{rel} = \vec{v}_i - \vec{v}_f` is the
particle velocity relative to the fluid, and the drag coefficient
:math:`C_d` follows the Schiller-Naumann correlation

.. math::

   C_d = \frac{24}{Re}\left(1 + 0.15\, Re^{0.687}\right), \qquad
   Re = \frac{\rho_f\, |\vec{v}_{rel}|\, (2 r_i)}{\mu_f}

with :math:`Re` the particle Reynolds number based on the diameter
:math:`2 r_i`, the fluid density :math:`\rho_f`, and the dynamic viscosity
of the fluid :math:`\mu_f`.  In the low-Reynolds-number limit
(:math:`Re \rightarrow 0`) the correlation reduces to
:math:`C_d = 24/Re` and the force becomes the Stokes drag
:math:`\vec{F}_i = -6 \pi \mu_f r_i \vec{v}_{rel}`.

By default the fluid is at rest.  The optional *velocity* keyword sets a
uniform fluid velocity :math:`\vec{v}_f`, so the drag is computed
from the particle velocity relative to the moving fluid.  This fix only
applies a drag force; buoyancy and gravity (if desired) must be added
separately, e.g. with :doc:`fix gravity <fix_gravity>`.

Restrictions
""""""""""""

This fix is part of the GRANULAR package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This fix requires that atoms store a radius as defined by the
:doc:`atom_style sphere <atom_style>` command.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.  No global or per-atom quantities are stored by
this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix.  This allows one to set at which level of the :doc:`r-RESPA
<run_style>` integrator the fix is modifying forces.  Default is the
outermost level.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

Related commands
""""""""""""""""

:doc:`fix viscous <fix_viscous>`,
:doc:`fix viscous/sphere <fix_viscous_sphere>`,
:doc:`fix gravity <fix_gravity>`

Default
"""""""

The fluid velocity is (0,0,0).
