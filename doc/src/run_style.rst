.. index:: run_style

run_style command
=================

Syntax
""""""

.. code-block:: LAMMPS

   run_style style args

* style = *verlet* or *verlet/split* or *respa* or *respa/omp*

  .. parsed-literal::

       *verlet* args = none
       *verlet/split* args = none
       *respa* args = N n1 n2 ... keyword values ...
         N = # of levels of rRESPA
         n1, n2, ... = loop factors between rRESPA levels (N-1 values)
         zero or more keyword/value pairings may be appended to the loop factors
         keyword = *bond* or *angle* or *dihedral* or *improper* or
                   *pair* or *inner* or *middle* or *outer* or *hybrid* or *kspace*
           *bond* value = M
             M = which level (1-N) to compute bond forces in
           *angle* value = M
             M = which level (1-N) to compute angle forces in
           *dihedral* value = M
             M = which level (1-N) to compute dihedral forces in
           *improper* value = M
             M = which level (1-N) to compute improper forces in
           *pair* value = M
             M = which level (1-N) to compute pair forces in
           *inner* values = M cut1 cut2
             M = which level (1-N) to compute pair inner forces in
             cut1 = inner cutoff between pair inner and
                    pair middle or outer  (distance units)
             cut2 = outer cutoff between pair inner and
                    pair middle or outer  (distance units)
           *middle* values = M cut1 cut2
             M = which level (1-N) to compute pair middle forces in
             cut1 = inner cutoff between pair middle and pair outer (distance units)
             cut2 = outer cutoff between pair middle and pair outer (distance units)
           *outer* value = M
             M = which level (1-N) to compute pair outer forces in
           *hybrid* values = M1 M2 ... (as many values as there are hybrid sub-styles
             M1 = which level (1-N) to compute the first pair_style hybrid sub-style in
             M2 = which level (1-N) to compute the second pair_style hybrid sub-style in
             M3,etc
           *kspace* value = M
             M = which level (1-N) to compute kspace forces in

Examples
""""""""

.. code-block:: LAMMPS

   run_style verlet
   run_style respa 4 2 2 2 bond 1 dihedral 2 pair 3 kspace 4
   run_style respa 4 2 2 2 bond 1 dihedral 2 inner 3 5.0 6.0 outer 4 kspace 4
   run_style respa 3 4 2 bond 1 hybrid 2 2 1 kspace 3

Description
"""""""""""

Choose the style of time integrator used for molecular dynamics
simulations performed by LAMMPS.

The *verlet* style is a standard velocity-Verlet integrator.

----------

The *verlet/split* style is also a velocity-Verlet integrator, but it
splits the force calculation within each timestep over 2 partitions of
processors.  See the :doc:`-partition command-line switch <Run_options>`
for info on how to run LAMMPS with multiple partitions.

Specifically, this style performs all computation except the
:doc:`kspace_style <kspace_style>` portion of the force field on the first
partition.  This include the :doc:`pair style <pair_style>`, :doc:`bond style <bond_style>`, :doc:`neighbor list building <neighbor>`,
:doc:`fixes <fix>` including time integration, and output.  The
:doc:`kspace_style <kspace_style>` portion of the calculation is
performed on the second partition.

This is most useful for the PPPM kspace_style when its performance on
a large number of processors degrades due to the cost of communication
in its 3d FFTs.  In this scenario, splitting your P total processors
into 2 subsets of processors, P1 in the first partition and P2 in the
second partition, can enable your simulation to run faster.  This is
because the long-range forces in PPPM can be calculated at the same
time as pair-wise and bonded forces are being calculated, and the FFTs
can actually speed up when running on fewer processors.

To use this style, you must define 2 partitions where P1 is a multiple
of P2.  Typically having P1 be 3x larger than P2 is a good choice.
The 3d processor layouts in each partition must overlay in the
following sense.  If P1 is a Px1 by Py1 by Pz1 grid, and P2 = Px2 by
Py2 by Pz2, then Px1 must be an integer multiple of Px2, and similarly
for Py1 a multiple of Py2, and Pz1 a multiple of Pz2.

Typically the best way to do this is to let the first partition choose
its onn optimal layout, then require the second partition's layout to
match the integer multiple constraint.  See the
:doc:`processors <processors>` command with its *part* keyword for a way
to control this, e.g.

.. code-block:: LAMMPS

   processors * * * part 1 2 multiple

You can also use the :doc:`partition <partition>` command to explicitly
specify the processor layout on each partition.  E.g. for 2 partitions
of 60 and 15 processors each:

.. code-block:: LAMMPS

   partition yes 1 processors 3 4 5
   partition yes 2 processors 3 1 5

When you run in 2-partition mode with the *verlet/split* style, the
thermodynamic data for the entire simulation will be output to the log
and screen file of the first partition, which are log.lammps.0 and
screen.0 by default; see the :doc:`-plog and -pscreen command-line switches <Run_options>` to change this.  The log and screen file
for the second partition will not contain thermodynamic output beyond the
first timestep of the run.

See the :doc:`Speed packages <Speed_packages>` doc page for performance
details of the speed-up offered by the *verlet/split* style.  One
important performance consideration is the assignment of logical
processors in the 2 partitions to the physical cores of a parallel
machine.  The :doc:`processors <processors>` command has options to
support this, and strategies are discussed in :doc:`Section 5 <Speed>` of the manual.

----------

The *respa* style implements the rRESPA multi-timescale integrator
:ref:`(Tuckerman) <Tuckerman3>` with N hierarchical levels, where level 1 is
the innermost loop (shortest timestep) and level N is the outermost
loop (largest timestep).  The loop factor arguments specify what the
looping factor is between levels.  N1 specifies the number of
iterations of level 1 for a single iteration of level 2, N2 is the
iterations of level 2 per iteration of level 3, etc.  N-1 looping
parameters must be specified.

Thus with a 4-level respa setting of "2 2 2" for the 3 loop factors,
you could choose to have bond interactions computed 8x per large
timestep, angle interactions computed 4x, pair interactions computed
2x, and long-range interactions once per large timestep.

The :doc:`timestep <timestep>` command sets the large timestep for the
outermost rRESPA level.  Thus if the 3 loop factors are "2 2 2" for
4-level rRESPA, and the outer timestep is set to 4.0 fmsec, then the
inner timestep would be 8x smaller or 0.5 fmsec.  All other LAMMPS
commands that specify number of timesteps (e.g. :doc:`thermo <thermo>`
for thermo output every N steps, :doc:`neigh_modify delay/every <neigh_modify>` parameters, :doc:`dump <dump>` every N
steps, etc) refer to the outermost timesteps.

The rRESPA keywords enable you to specify at what level of the
hierarchy various forces will be computed.  If not specified, the
defaults are that bond forces are computed at level 1 (innermost
loop), angle forces are computed where bond forces are, dihedral
forces are computed where angle forces are, improper forces are
computed where dihedral forces are, pair forces are computed at the
outermost level, and kspace forces are computed where pair forces are.
The inner, middle, outer forces have no defaults.

For fixes that support it, the rRESPA level at which a given fix is
active, can be selected through the :doc:`fix_modify <fix_modify>` command.

The *inner* and *middle* keywords take additional arguments for
cutoffs that are used by the pairwise force computations.  If the 2
cutoffs for *inner* are 5.0 and 6.0, this means that all pairs up to
6.0 apart are computed by the inner force.  Those between 5.0 and 6.0
have their force go ramped to 0.0 so the overlap with the next regime
(middle or outer) is smooth.  The next regime (middle or outer) will
compute forces for all pairs from 5.0 outward, with those from 5.0 to
6.0 having their value ramped in an inverse manner.

Note that you can use *inner* and *outer* without using *middle* to
split the pairwise computations into two portions instead of three.
Unless you are using a very long pairwise cutoff, a 2-way split is
often faster than a 3-way split, since it avoids too much duplicate
computation of pairwise interactions near the intermediate cutoffs.

Also note that only a few pair potentials support the use of the
*inner* and *middle* and *outer* keywords.  If not, only the *pair*
keyword can be used with that pair style, meaning all pairwise forces
are computed at the same rRESPA level.  See the doc pages for
individual pair styles for details.

Another option for using pair potentials with rRESPA is with the
*hybrid* keyword, which requires the use of the :doc:`pair_style hybrid or hybrid/overlay <pair_hybrid>` command.  In this scenario, different
sub-styles of the hybrid pair style are evaluated at different rRESPA
levels.  This can be useful, for example, to set different timesteps
for hybrid coarse-grained/all-atom models.  The *hybrid* keyword
requires as many level assignments as there are hybrid sub-styles,
which assigns each sub-style to a rRESPA level, following their order
of definition in the pair_style command. Since the *hybrid* keyword
operates on pair style computations, it is mutually exclusive with
either the *pair* or the *inner*\ /\ *middle*\ /\ *outer* keywords.

When using rRESPA (or for any MD simulation) care must be taken to
choose a timestep size(s) that insures the Hamiltonian for the chosen
ensemble is conserved.  For the constant NVE ensemble, total energy
must be conserved.  Unfortunately, it is difficult to know *a priori*
how well energy will be conserved, and a fairly long test simulation
(~10 ps) is usually necessary in order to verify that no long-term
drift in energy occurs with the trial set of parameters.

With that caveat, a few rules-of-thumb may be useful in selecting
*respa* settings.  The following applies mostly to biomolecular
simulations using the CHARMM or a similar all-atom force field, but
the concepts are adaptable to other problems.  Without SHAKE, bonds
involving hydrogen atoms exhibit high-frequency vibrations and require
a timestep on the order of 0.5 fmsec in order to conserve energy.  The
relatively inexpensive force computations for the bonds, angles,
impropers, and dihedrals can be computed on this innermost 0.5 fmsec
step.  The outermost timestep cannot be greater than 4.0 fmsec without
risking energy drift.  Smooth switching of forces between the levels
of the rRESPA hierarchy is also necessary to avoid drift, and a 1-2
angstrom "healing distance" (the distance between the outer and inner
cutoffs) works reasonably well.  We thus recommend the following
settings for use of the *respa* style without SHAKE in biomolecular
simulations:

.. code-block:: LAMMPS

   timestep  4.0
   run_style respa 4 2 2 2 inner 2 4.5 6.0 middle 3 8.0 10.0 outer 4

With these settings, users can expect good energy conservation and
roughly a 2.5 fold speedup over the *verlet* style with a 0.5 fmsec
timestep.

If SHAKE is used with the *respa* style, time reversibility is lost,
but substantially longer time steps can be achieved.  For biomolecular
simulations using the CHARMM or similar all-atom force field, bonds
involving hydrogen atoms exhibit high frequency vibrations and require
a time step on the order of 0.5 fmsec in order to conserve energy.
These high frequency modes also limit the outer time step sizes since
the modes are coupled.  It is therefore desirable to use SHAKE with
respa in order to freeze out these high frequency motions and increase
the size of the time steps in the respa hierarchy.  The following
settings can be used for biomolecular simulations with SHAKE and
rRESPA:

.. code-block:: LAMMPS

   fix             2 all shake 0.000001 500 0 m 1.0 a 1
   timestep        4.0
   run_style       respa 2 2 inner 1 4.0 5.0 outer 2

With these settings, users can expect good energy conservation and
roughly a 1.5 fold speedup over the *verlet* style with SHAKE and a
2.0 fmsec timestep.

For non-biomolecular simulations, the *respa* style can be
advantageous if there is a clear separation of time scales - fast and
slow modes in the simulation.  For example, a system of slowly-moving
charged polymer chains could be setup as follows:

.. code-block:: LAMMPS

   timestep 4.0
   run_style respa 2 8

This is two-level rRESPA with an 8x difference between the short and
long timesteps.  The bonds, angles, dihedrals will be computed every
0.5 fs (assuming real units), while the pair and kspace interactions
will be computed once every 4 fs.  These are the default settings for
each kind of interaction, so no additional keywords are necessary.

Even a LJ system can benefit from rRESPA if the interactions are
divided by the inner, middle and outer keywords.  A 2-fold or more
speedup can be obtained while maintaining good energy conservation.
In real units, for a pure LJ fluid at liquid density, with a sigma of
3.0 angstroms, and epsilon of 0.1 Kcal/mol, the following settings
seem to work well:

.. code-block:: LAMMPS

   timestep  36.0
   run_style respa 3 3 4 inner 1 3.0 4.0 middle 2 6.0 7.0 outer 3

----------

The *respa/omp* style is a variant of *respa* adapted for use with
pair, bond, angle, dihedral, improper, or kspace styles with an *omp*
suffix. It is functionally equivalent to *respa* but performs
additional operations required for managing *omp* styles.  For more on
*omp* styles see the :doc:`Speed omp <Speed_omp>` doc page.  Accelerated
styles take the same arguments and should produce the same results,
except for round-off and precision issues.

You can specify *respa/omp* explicitly in your input script, or you
can use the :doc:`-suffix command-line switch <Run_options>` when you
invoke LAMMPS, or you can use the :doc:`suffix <suffix>` command in your
input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.

----------

Restrictions
""""""""""""

The *verlet/split* style can only be used if LAMMPS was built with the
REPLICA package. Correspondingly the *respa/omp* style is available
only if the USER-OMP package was included. See the :doc:`Build package <Build_package>` doc page for more info.

Whenever using rRESPA, the user should experiment with trade-offs in
speed and accuracy for their system, and verify that they are
conserving energy to adequate precision.

Related commands
""""""""""""""""

:doc:`timestep <timestep>`, :doc:`run <run>`

Default
"""""""

.. code-block:: LAMMPS

   run_style verlet

For run_style respa, the default assignment of interactions
to rRESPA levels is as follows:

* bond forces = level 1 (innermost loop)
* angle forces = same level as bond forces
* dihedral forces = same level as angle forces
* improper forces = same level as dihedral forces
* pair forces = level N (outermost level)
* kspace forces = same level as pair forces
* inner, middle, outer forces = no default

----------

.. _Tuckerman3:

**(Tuckerman)** Tuckerman, Berne and Martyna, J Chem Phys, 97, p 1990
(1992).
