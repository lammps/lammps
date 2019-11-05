Polarizable models
==================

In polarizable force fields the charge distributions in molecules and
materials respond to their electrostatic environments. Polarizable
systems can be simulated in LAMMPS using three methods:

* the fluctuating charge method, implemented in the :doc:`QEQ <fix_qeq>`
  package,
* the adiabatic core-shell method, implemented in the
  :doc:`CORESHELL <Howto_coreshell>` package,
* the thermalized Drude dipole method, implemented in the
  :doc:`USER-DRUDE <Howto_drude>` package.

The fluctuating charge method calculates instantaneous charges on
interacting atoms based on the electronegativity equalization
principle. It is implemented in the :doc:`fix qeq <fix_qeq>` which is
available in several variants. It is a relatively efficient technique
since no additional particles are introduced. This method allows for
charge transfer between molecules or atom groups. However, because the
charges are located at the interaction sites, off-plane components of
polarization cannot be represented in planar molecules or atom groups.

The two other methods share the same basic idea: polarizable atoms are
split into one core atom and one satellite particle (called shell or
Drude particle) attached to it by a harmonic spring.  Both atoms bear
a charge and they represent collectively an induced electric dipole.
These techniques are computationally more expensive than the QEq
method because of additional particles and bonds. These two
charge-on-spring methods differ in certain features, with the
core-shell model being normally used for ionic/crystalline materials,
whereas the so-called Drude model is normally used for molecular
systems and fluid states.

The core-shell model is applicable to crystalline materials where the
high symmetry around each site leads to stable trajectories of the
core-shell pairs. However, bonded atoms in molecules can be so close
that a core would interact too strongly or even capture the Drude
particle of a neighbor. The Drude dipole model is relatively more
complex in order to remedy this and other issues. Specifically, the
Drude model includes specific thermostatting of the core-Drude pairs
and short-range damping of the induced dipoles.

The three polarization methods can be implemented through a
self-consistent calculation of charges or induced dipoles at each
timestep. In the fluctuating charge scheme this is done by the matrix
inversion method in :doc:`fix qeq/point <fix_qeq>`, but for core-shell
or Drude-dipoles the relaxed-dipoles technique would require an slow
iterative procedure. These self-consistent solutions yield accurate
trajectories since the additional degrees of freedom representing
polarization are massless.  An alternative is to attribute a mass to
the additional degrees of freedom and perform time integration using
an extended Lagrangian technique. For the fluctuating charge scheme
this is done by :doc:`fix qeq/dynamic <fix_qeq>`, and for the
charge-on-spring models by the methods outlined in the next two
sections. The assignment of masses to the additional degrees of
freedom can lead to unphysical trajectories if care is not exerted in
choosing the parameters of the polarizable models and the simulation
conditions.

In the core-shell model the vibration of the shells is kept faster
than the ionic vibrations to mimic the fast response of the
polarizable electrons.  But in molecular systems thermalizing the
core-Drude pairs at temperatures comparable to the rest of the
simulation leads to several problems (kinetic energy transfer, too
short a timestep, etc.) In order to avoid these problems the relative
motion of the Drude particles with respect to their cores is kept
"cold" so the vibration of the core-Drude pairs is very slow,
approaching the self-consistent regime.  In both models the
temperature is regulated using the velocities of the center of mass of
core+shell (or Drude) pairs, but in the Drude model the actual
relative core-Drude particle motion is thermostatted separately as
well.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
