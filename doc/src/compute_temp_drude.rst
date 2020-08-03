.. index:: compute temp/drude

compute temp/drude command
==========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID temp/drude

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/drude = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute TDRUDE all temp/drude

Description
"""""""""""

Define a computation that calculates the temperatures of core-Drude
pairs. This compute is designed to be used with the :doc:`thermalized Drude oscillator model <Howto_drude>`.  Polarizable models in LAMMPS
are described on the :doc:`Howto polarizable <Howto_polarizable>` doc
page.

Drude oscillators consist of a core particle and a Drude particle
connected by a harmonic bond, and the relative motion of these Drude
oscillators is usually maintained cold by a specific thermostat that
acts on the relative motion of the core-Drude particle
pairs. Therefore, because LAMMPS considers Drude particles as normal
atoms in its default temperature compute (:doc:`compute temp <compute_temp>` command), the reduced temperature of the
core-Drude particle pairs is not calculated correctly.

By contrast, this compute calculates the temperature of the cores
using center-of-mass velocities of the core-Drude pairs, and the
reduced temperature of the Drude particles using the relative
velocities of the Drude particles with respect to their cores.
Non-polarizable atoms are considered as cores.  Their velocities
contribute to the temperature of the cores.

**Output info:**

This compute calculates a global scalar (the temperature) and a global
vector of length 6, which can be accessed by indices 1-6, whose components
are

1. temperature of the centers of mass (temperature units)
2. temperature of the dipoles (temperature units)
3. number of degrees of freedom of the centers of mass
4. number of degrees of freedom of the dipoles
5. kinetic energy of the centers of mass (energy units)
6. kinetic energy of the dipoles (energy units)

These values can be used by any command that uses global scalar or
vector values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

Both the scalar value and the first two values of the vector
calculated by this compute are "intensive".  The other 4 vector values
are "extensive".

Restrictions
""""""""""""

The number of degrees of freedom contributing to the temperature is
assumed to be constant for the duration of the run unless the
*fix_modify* command sets the option *dynamic yes*\ .

Related commands
""""""""""""""""

:doc:`fix drude <fix_drude>`, :doc:`fix langevin/drude <fix_langevin_drude>`, :doc:`fix drude/transform <fix_drude_transform>`, :doc:`pair_style thole <pair_thole>`, :doc:`compute temp <compute_temp>`

**Default:** none
