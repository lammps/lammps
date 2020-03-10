.. index:: compute stress/mop

compute stress/mop command
==========================

compute stress/mop/profile command
==================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID style dir args keywords ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* style = *stress/mop* or *stress/mop/profile*
* dir = *x* or *y* or *z* is the direction normal to the plane
* args = argument specific to the compute style
* keywords = *kin* or *conf* or *total* (one of more can be specified)

.. parsed-literal::

     *stress/mop* args = pos
       pos = *lower* or *center* or *upper* or coordinate value (distance units) is the position of the plane
     *stress/mop/profile* args = origin delta
       origin = *lower* or *center* or *upper* or coordinate value (distance units) is the position of the first plane
       delta = value (distance units) is the distance between planes

   compute 1 all stress/mop x lower total
   compute 1 liquid stress/mop z 0.0 kin conf
   fix 1 all ave/time 10 1000 10000 c_1[\*] file mop.time
   fix 1 all ave/time 10 1000 10000 c_1[2] file mop.time

   compute 1 all stress/mop/profile x lower 0.1 total
   compute 1 liquid stress/mop/profile z 0.0 0.25 kin conf
   fix 1 all ave/time 500 20 10000 c_1[\*] ave running overwrite file mopp.time mode vector

Description
"""""""""""

Compute *stress/mop* and compute *stress/mop/profile* define computations that
calculate components of the local stress tensor using the method of
planes :ref:`(Todd) <mop-todd>`.  Specifically in compute *stress/mop* calculates 3
components are computed in directions *dir*\ ,\ *x*\ ; *dir*\ ,\ *y*\ ; and
*dir*\ ,\ *z*\ ; where *dir* is the direction normal to the plane, while
in compute *stress/mop/profile* the profile of the stress is computed.

Contrary to methods based on histograms of atomic stress (i.e. using
:doc:`compute stress/atom <compute_stress_atom>`), the method of planes is
compatible with mechanical balance in heterogeneous systems and at
interfaces :ref:`(Todd) <mop-todd>`.

The stress tensor is the sum of a kinetic term and a configurational
term, which are given respectively by Eq. (21) and Eq. (16) in
:ref:`(Todd) <mop-todd>`. For the kinetic part, the algorithm considers that
atoms have crossed the plane if their positions at times t-dt and t are
one on either side of the plane, and uses the velocity at time t-dt/2
given by the velocity-Verlet algorithm.

Between one and three keywords can be used to indicate which
contributions to the stress must be computed: kinetic stress (kin),
configurational stress (conf), and/or total stress (total).

NOTE 1: The configurational stress is computed considering all pairs of atoms where at least one atom belongs to group group-ID.

NOTE 2: The local stress does not include any Lennard-Jones tail
corrections to the pressure added by the :doc:`pair_modify tail yes <pair_modify>` command, since those are contributions to the global system pressure.

**Output info:**

Compute *stress/mop* calculates a global vector (indices starting at 1), with 3
values for each declared keyword (in the order the keywords have been
declared). For each keyword, the stress tensor components are ordered as
follows: stress\_dir,x, stress\_dir,y, and stress\_dir,z.

Compute *stress/mop/profile* instead calculates a global array, with 1 column
giving the position of the planes where the stress tensor was computed,
and with 3 columns of values for each declared keyword (in the order the
keywords have been declared). For each keyword, the profiles of stress
tensor components are ordered as follows: stress\_dir,x; stress\_dir,y;
and stress\_dir,z.

The values are in pressure :doc:`units <units>`.

The values produced by this compute can be accessed by various :doc:`output commands <Howto_output>`. For instance, the results can be written to a file using the :doc:`fix ave/time <fix_ave_time>` command. Please see the example in the examples/USER/mop folder.

Restrictions
""""""""""""

These styles are part of the USER-MISC package. They are only enabled if
LAMMPS is built with that package. See the :doc:`Build package <Build_package>`
doc page on for more info.

The method is only implemented for 3d orthogonal simulation boxes whose
size does not change in time, and axis-aligned planes.

The method only works with two-body pair interactions, because it
requires the class method pair->single() to be implemented. In
particular, it does not work with more than two-body pair interactions,
intra-molecular interactions, and long range (kspace) interactions.

Related commands
""""""""""""""""

:doc:`compute stress/atom <compute_stress_atom>`

**Default:** none

----------

.. _mop-todd:

**(Todd)** B. D. Todd, Denis J. Evans, and Peter J. Daivis: "Pressure tensor for inhomogeneous fluids",
Phys. Rev. E 52, 1627 (1995).
