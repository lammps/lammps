.. index:: fix orient/fcc

fix orient/fcc command
======================

fix orient/bcc command
======================


.. parsed-literal::

   fix ID group-ID orient/fcc nstats dir alat dE cutlo cuthi file0 file1
   fix ID group-ID orient/bcc nstats dir alat dE cutlo cuthi file0 file1

* ID, group-ID are documented in :doc:`fix <fix>` command
* nstats = print stats every this many steps, 0 = never
* dir = 0/1 for which crystal is used as reference
* alat = fcc/bcc cubic lattice constant (distance units)
* dE = energy added to each atom (energy units)
* cutlo,cuthi = values between 0.0 and 1.0, cutlo < cuthi
* file0,file1 = files that specify orientation of each grain

Examples
""""""""


.. parsed-literal::

   fix gb all orient/fcc 0 1 4.032008 0.001 0.25 0.75 xi.vec chi.vec
   fix gb all orient/bcc 0 1 2.882 0.001 0.25 0.75 ngb.left ngb.right

Description
"""""""""""

The fix applies an orientation-dependent force to atoms near a planar
grain boundary which can be used to induce grain boundary migration
(in the direction perpendicular to the grain boundary plane).  The
motivation and explanation of this force and its application are
described in :ref:`(Janssens) <Janssens>`. The adaptation to bcc crystals
is described in :ref:`(Wicaksono1) <Wicaksono1>`. The computed force is only
applied to atoms in the fix group.

The basic idea is that atoms in one grain (on one side of the
boundary) have a potential energy dE added to them.  Atoms in the
other grain have 0.0 potential energy added.  Atoms near the boundary
(whose neighbor environment is intermediate between the two grain
orientations) have an energy between 0.0 and dE added.  This creates
an effective driving force to reduce the potential energy of atoms
near the boundary by pushing them towards one of the grain
orientations.  For dir = 1 and dE > 0, the boundary will thus move so
that the grain described by file0 grows and the grain described by
file1 shrinks.  Thus this fix is designed for simulations of two-grain
systems, either with one grain boundary and free surfaces parallel to
the boundary, or a system with periodic boundary conditions and two
equal and opposite grain boundaries.  In either case, the entire
system can displace during the simulation, and such motion should be
accounted for in measuring the grain boundary velocity.

The potential energy added to atom I is given by these formulas

.. math::

  \xi_{i} = \sum_{j=1}^{12} \left| \br_{j} - \br_{j}^{\rm I} \right|

>>>image was here
  \xi_{\rm IJ} = \sum_{j=1}^{12} \left| \br_{j}^{\rm J} - \br_{j}^{\rm I} \right|


which are fully explained in :ref:`(Janssens) <Janssens>`.  For fcc crystals
this order parameter Xi for atom I in equation (1) is a sum over the
12 nearest neighbors of atom I. For bcc crystals it is the
corresponding sum of the 8 nearest neighbors. Rj is the vector from
atom I to its neighbor J, and RIj is a vector in the reference
(perfect) crystal.  That is, if dir = 0/1, then RIj is a vector to an
atom coord from file 0/1.  Equation (2) gives the expected value of
the order parameter XiIJ in the other grain.  Hi and lo cutoffs are
defined in equations (3) and (4), using the input parameters *cutlo*
and *cuthi* as thresholds to avoid adding grain boundary energy when
the deviation in the order parameter from 0 or 1 is small (e.g. due to
thermal fluctuations in a perfect crystal).  The added potential
energy Ui for atom I is given in equation (6) where it is interpolated
between 0 and dE using the two threshold Xi values and the Wi value of
equation (5).

The derivative of this energy expression gives the force on each atom
which thus depends on the orientation of its neighbors relative to the
2 grain orientations.  Only atoms near the grain boundary feel a net
force which tends to drive them to one of the two grain orientations.

In equation (1), the reference vector used for each neighbor is the
reference vector closest to the actual neighbor position.  This means
it is possible two different neighbors will use the same reference
vector.  In such cases, the atom in question is far from a perfect
orientation and will likely receive the full dE addition, so the
effect of duplicate reference vector usage is small.

The *dir* parameter determines which grain wants to grow at the
expense of the other.  A value of 0 means the first grain will shrink;
a value of 1 means it will grow.  This assumes that *dE* is positive.
The reverse will be true if *dE* is negative.

The *alat* parameter is the cubic lattice constant for the fcc or bcc
material and is only used to compute a cutoff distance of 1.57 \* alat
/ sqrt(2) for finding the 12 or 8 nearest neighbors of each atom
(which should be valid for an fcc or bcc crystal).  A longer/shorter
cutoff can be imposed by adjusting *alat*\ .  If a particular atom has
less than 12 or 8 neighbors within the cutoff, the order parameter of
equation (1) is effectively multiplied by 12 or 8 divided by the
actual number of neighbors within the cutoff.

The *dE* parameter is the maximum amount of additional energy added to
each atom in the grain which wants to shrink.

The *cutlo* and *cuthi* parameters are used to reduce the force added
to bulk atoms in each grain far away from the boundary.  An atom in
the bulk surrounded by neighbors at the ideal grain orientation would
compute an order parameter of 0 or 1 and have no force added.
However, thermal vibrations in the solid will cause the order
parameters to be greater than 0 or less than 1.  The cutoff parameters
mask this effect, allowing forces to only be added to atoms with
order-parameters between the cutoff values.

*File0* and *file1* are filenames for the two grains which each
contain 6 vectors (6 lines with 3 values per line) which specify the
grain orientations.  Each vector is a displacement from a central atom
(0,0,0) to a nearest neighbor atom in an fcc lattice at the proper
orientation.  The vector lengths should all be identical since an fcc
lattice has a coordination number of 12.  Only 6 are listed due to
symmetry, so the list must include one from each pair of
equal-and-opposite neighbors.  A pair of orientation files for a
Sigma=5 tilt boundary are shown below. A tutorial that can help for
writing the orientation files is given in :ref:`(Wicaksono2) <Wicaksono2>`

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix\_modify <fix_modify>` *energy* option is supported by this
fix to add the potential energy of atom interactions with the grain
boundary driving force to the system's potential energy as part of
:doc:`thermodynamic output <thermo_style>`.

The :doc:`fix\_modify <fix_modify>` *respa* option is supported by these
fixes. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator a fix is adding its forces. Default is the outermost level.

This fix calculates a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the potential
energy change due to this fix.  The scalar value calculated by this
fix is "extensive".

This fix also calculates a per-atom array which can be accessed by
various :doc:`output commands <Howto_output>`.  The array stores the
order parameter Xi and normalized order parameter (0 to 1) for each
atom.  The per-atom values can be accessed on any timestep.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the MISC package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This fix should only be used with fcc or bcc lattices.

Related commands
""""""""""""""""

:doc:`fix\_modify <fix_modify>`

**Default:** none


----------


.. _Janssens:



**(Janssens)** Janssens, Olmsted, Holm, Foiles, Plimpton, Derlet, Nature
Materials, 5, 124-127 (2006).

.. _Wicaksono1:



**(Wicaksono1)** Wicaksono, Sinclair, Militzer, Computational Materials
Science, 117, 397-405 (2016).

.. _Wicaksono2:



**(Wicaksono2)** Wicaksono, figshare,
https://dx.doi.org/10.6084/m9.figshare.1488628.v1 (2015).


----------


For illustration purposes, here are example files that specify a
Sigma=5 <100> tilt boundary.  This is for a lattice constant of 3.5706
Angs.

file0:


.. parsed-literal::

        0.798410432046075    1.785300000000000    1.596820864092150
       -0.798410432046075    1.785300000000000   -1.596820864092150
        2.395231296138225    0.000000000000000    0.798410432046075
        0.798410432046075    0.000000000000000   -2.395231296138225
        1.596820864092150    1.785300000000000   -0.798410432046075
        1.596820864092150   -1.785300000000000   -0.798410432046075

file1:


.. parsed-literal::

       -0.798410432046075    1.785300000000000    1.596820864092150
        0.798410432046075    1.785300000000000   -1.596820864092150
        0.798410432046075    0.000000000000000    2.395231296138225
        2.395231296138225    0.000000000000000   -0.798410432046075
        1.596820864092150    1.785300000000000    0.798410432046075
        1.596820864092150   -1.785300000000000    0.798410432046075


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
