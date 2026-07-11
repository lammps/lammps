Calculate temperature
=====================

Temperature is computed as kinetic energy divided by some number of
degrees of freedom (and the Boltzmann constant).  Since kinetic energy
is a function of particle velocity, there is often a need to
distinguish between a particle's advection velocity (due to some
aggregate motion of particles) and its thermal velocity.  The sum of
the two is the particle's total velocity, but the latter is often what
is wanted to compute a temperature.

LAMMPS has several options for computing temperatures, any of which
can be used in :doc:`thermostatting <Howto_thermostat>` and
:doc:`barostatting <Howto_barostat>`.  These :doc:`compute commands <compute>` calculate temperature:

* :doc:`compute temp <compute_temp>`
* :doc:`compute temp/sphere <compute_temp_sphere>`
* :doc:`compute temp/asphere <compute_temp_asphere>`
* :doc:`compute temp/com <compute_temp_com>`
* :doc:`compute temp/deform <compute_temp_deform>`
* :doc:`compute temp/partial <compute_temp_partial>`
* :doc:`compute temp/profile <compute_temp_profile>`
* :doc:`compute temp/ramp <compute_temp_ramp>`
* :doc:`compute temp/region <compute_temp_region>`

All but the first 3 calculate velocity biases directly (e.g. advection
velocities) that are removed when computing the thermal temperature.
:doc:`Compute temp/sphere <compute_temp_sphere>` and :doc:`compute temp/asphere <compute_temp_asphere>` compute kinetic energy for
finite-size particles that includes rotational degrees of freedom.
They both allow for velocity biases indirectly, via an optional extra
argument which is another temperature compute that subtracts a
velocity bias.  This allows the translational velocity of spherical or
aspherical particles to be adjusted in prescribed ways.
