.. index:: fix drude/transform/direct

fix drude/transform/direct command
==================================

fix drude/transform/inverse command
===================================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID style keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *drude/transform/direct* or *drude/transform/inverse*

Examples
""""""""


.. parsed-literal::

   fix 3 all drude/transform/direct
   fix 1 all drude/transform/inverse

Description
"""""""""""

Transform the coordinates of Drude oscillators from real to reduced
and back for thermalizing the Drude oscillators as described in
:ref:`(Lamoureux) <Lamoureux1>` using a Nose-Hoover thermostat.  This fix is
designed to be used with the :doc:`thermalized Drude oscillator model <Howto_drude>`.  Polarizable models in LAMMPS are described
on the :doc:`Howto polarizable <Howto_polarizable>` doc page.

Drude oscillators are a pair of atoms representing a single
polarizable atom.  Ideally, the mass of Drude particles would vanish
and their positions would be determined self-consistently by iterative
minimization of the energy, the cores' positions being fixed.  It is
however more efficient and it yields comparable results, if the Drude
oscillators (the motion of the Drude particle relative to the core)
are thermalized at a low temperature.  In that case, the Drude
particles need a small mass.

The thermostats act on the reduced degrees of freedom, which are
defined by the following equations.  Note that in these equations
upper case denotes atomic or center of mass values and lower case
denotes Drude particle or dipole values. Primes denote the transformed
(reduced) values, while bare letters denote the original values.

Masses:

.. math::

    M' = M + m


.. math::

    m' = \frac {M\, m } {M'}

Positions:

.. math::

    X' = \frac {M\, X + m\, x} {M'}


.. math::

    x' = x - X

Velocities:

.. math::

    V' = \frac {M\, V + m\, v} {M'}


.. math::

    v' = v - V

Forces:

.. math::

    F' = F + f


.. math::

    f' = \frac { M\, f - m\, F} {M'}

This transform conserves the total kinetic energy

.. math::

    \frac 1 2 \, (M\, V^2\ + m\, v^2)
   = \frac 1 2 \, (M'\, V'^2\ + m'\, v'^2)

and the virial defined with absolute positions

.. math::

    X\, F + x\, f = X'\, F' + x'\, f'


----------


This fix requires each atom know whether it is a Drude particle or
not.  You must therefore use the :doc:`fix drude <fix_drude>` command to
specify the Drude status of each atom type.

.. note::

   only the Drude core atoms need to be in the group specified for
   this fix. A Drude electron will be transformed together with its core
   even if it is not itself in the group.  It is safe to include Drude
   electrons or non-polarizable atoms in the group. The non-polarizable
   atoms will simply not be transformed.


----------


This fix does NOT perform time integration. It only transform masses,
coordinates, velocities and forces. Thus you must use separate time
integration fixes, like :doc:`fix nve <fix_nve>` or :doc:`fix npt <fix_nh>` to actually update the velocities and positions of
atoms.  In order to thermalize the reduced degrees of freedom at
different temperatures, two Nose-Hoover thermostats must be defined,
acting on two distinct groups.

.. note::

   The *fix drude/transform/direct* command must appear before any
   Nose-Hoover thermostatting fixes.  The *fix drude/transform/inverse*
   command must appear after any Nose-Hoover thermostatting fixes.

Example:


.. parsed-literal::

   fix fDIRECT all drude/transform/direct
   fix fNVT gCORES nvt temp 300.0 300.0 100
   fix fNVT gDRUDES nvt temp 1.0 1.0 100
   fix fINVERSE all drude/transform/inverse
   compute TDRUDE all temp/drude
   thermo_style custom step cpu etotal ke pe ebond ecoul elong press vol temp c_TDRUDE[1] c_TDRUDE[2]

In this example, *gCORES* is the group of the atom cores and *gDRUDES*
is the group of the Drude particles (electrons). The centers of mass
of the Drude oscillators will be thermostatted at 300.0 and the
internal degrees of freedom will be thermostatted at 1.0.  The
temperatures of cores and Drude particles, in center-of-mass and
relative coordinates, are calculated using :doc:`compute temp/drude <compute_temp_drude>`

In addition, if you want to use a barostat to simulate a system at
constant pressure, only one of the Nose-Hoover fixes must be *npt*\ ,
the other one should be *nvt*\ . You must add a *compute temp/com* and a
*fix\_modify* command so that the temperature of the *npt* fix be just
that of its group (the Drude cores) but the pressure be the overall
pressure *thermo\_press*.

Example:


.. parsed-literal::

   compute cTEMP_CORE gCORES temp/com
   fix fDIRECT all drude/transform/direct
   fix fNPT gCORES npt temp 298.0 298.0 100 iso 1.0 1.0 500
   fix_modify fNPT temp cTEMP_CORE press thermo_press
   fix fNVT gDRUDES nvt temp 5.0 5.0 100
   fix fINVERSE all drude/transform/inverse

In this example, *gCORES* is the group of the atom cores and *gDRUDES*
is the group of the Drude particles. The centers of mass of the Drude
oscillators will be thermostatted at 298.0 and the internal degrees of
freedom will be thermostatted at 5.0. The whole system will be
barostatted at 1.0.

In order to avoid the flying ice cube problem (irreversible transfer
of linear momentum to the center of mass of the system), you may need
to add a *fix momentum* command:


.. parsed-literal::

   fix fMOMENTUM all momentum 100 linear 1 1 1


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix drude <fix_drude>`,
:doc:`fix langevin/drude <fix_langevin_drude>`,
:doc:`compute temp/drude <compute_temp_drude>`,
:doc:`pair_style thole <pair_thole>`

**Default:** none


----------


.. _Lamoureux1:



**(Lamoureux)** Lamoureux and Roux, J Chem Phys, 119, 3025-3039 (2003).
