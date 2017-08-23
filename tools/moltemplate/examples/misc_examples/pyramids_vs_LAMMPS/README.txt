   Description:

This is a simulation of pyramid-shaped objects resting on an immobile surface
(resembling graphene).  Each pyramid is built from spherical particles stacked
like cannon-balls (or fruit).  Ordinarily, the stack does not move
because the particles at the ground layer are immobilized.  However,
given an initial (small) perturbation the pyramids collapse in an avalanche.

(In this example, the perturbation is due to shock because we (intentionally)
 did not minimize the system before starting the simulation.  This shock
 causes an avalanche to begin approximately 5000 timesteps later.)

The particles roll down the pyramid and bounce off the "ground". The bouncing
is due to a repulsive external force which is added artificially.
(See the "run.in" file.)  The simulation looks weird without something
to bounce off of.  So I added a graphene surface at the bottom as scenery.
(It does not exert any force on the atoms.)

(Random comment: This could be a fun example to illustrate the Boltzmann
 distribution.  Because there is no damping, in a small region, I'm guessing
 the particle heights should eventually approach the Boltzmann distribution
 for some temperature consistent with the initial potential energy of the
 system.)
