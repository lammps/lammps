.. index:: fix plumed

fix plumed command
==================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID plumed keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* plumed = style name of this fix command
* keyword = *plumedfile* or *outfile* or *path_integral* or *pimd_fix*

  .. parsed-literal::

       *plumedfile* arg = name of PLUMED input file to use (default: NULL)
       *outfile* arg = name of file on which to write the PLUMED log (default: NULL)
       *path_integral* arg = *off*, *centroid*, *bead_mean*, or *bead_density* (default: off)
       *pimd_fix* arg = ID of the coupled fix pimd/langevin

Examples
""""""""

.. code-block:: LAMMPS

   fix pl all plumed plumedfile plumed.dat outfile p.log
   fix pl all plumed plumedfile plumed.dat outfile p.log path_integral centroid pimd_fix fpimd
   fix pl all plumed plumedfile plumed.dat outfile p.log path_integral bead_mean pimd_fix fpimd
   fix pl all plumed plumedfile plumed.dat outfile p.log path_integral bead_density pimd_fix fpimd

Description
"""""""""""

This fix instructs LAMMPS to call the `PLUMED <plumedhome_>`_ library, which
allows one to perform various forms of trajectory analysis on the fly
and to also use methods such as umbrella sampling and metadynamics to
enhance the sampling of phase space.

The documentation included here only describes the fix plumed command
itself.  This command is LAMMPS specific, whereas most of the
functionality implemented in PLUMED will work with a range of MD codes,
and when PLUMED is used as a stand alone code for analysis.  The full
`documentation for PLUMED <plumeddocs_>`_ is available online and included
in the PLUMED source code.  The PLUMED library development is hosted at
`https://github.com/plumed/plumed2 <https://github.com/plumed/plumed2>`_
A detailed discussion of the code can be found in :ref:`(Tribello) <Tribello>`.

There is an example input for using this package with LAMMPS in the
examples/PACKAGES/plumed directory.

----------

The command to make LAMMPS call PLUMED during a run requires two keyword
value pairs pointing to the PLUMED input file and an output file for the
PLUMED log. The user must specify these arguments every time PLUMED is
to be used.  Furthermore, the fix plumed command should appear in the
LAMMPS input file **after** relevant input parameters (e.g. the timestep)
have been set.

The *group-ID* entry is ignored. LAMMPS will always pass all the atoms
to PLUMED and there can only be one instance of the plumed fix at a
time. The way the plumed fix is implemented ensures that the minimum
amount of information required is communicated.  Furthermore, PLUMED
supports multiple, completely independent collective variables, multiple
independent biases and multiple independent forms of analysis.  There is
thus really no restriction in functionality by only allowing only one
plumed fix in the LAMMPS input.

The *plumedfile* keyword allows the user to specify the name of the
PLUMED input file.  Instructions as to what should be included in a
plumed input file can be found in the `documentation for PLUMED
<plumeddocs_>`_

The *outfile* keyword allows the user to specify the name of a file in
which to output the PLUMED log.  This log file normally just repeats the
information that is contained in the input file to confirm it was
correctly read and parsed.  The names of the files in which the results
are stored from the various analysis options performed by PLUMED will
be specified by the user in the PLUMED input file.

.. versionadded:: TBD

.. versionchanged:: TBD
   The *centroid*, *bead_mean*, and *bead_density* modes support normal-mode
   PIMD with NVT, NPH, or NPT and multiple MPI ranks per bead.

For these path-integral modes, :doc:`fix pimd/langevin <fix_pimd>` uses
the ring-polymer inverse temperature :math:`\beta/P` with unscaled physical
forces.  For a fixed physical bias :math:`U_B` reported by this fix, the
stationary configurational density is

.. math::

   \rho_B(X)\propto\exp[-(\beta/P)(H_{\mathrm{ring}}(X)+P U_B(X))].

Thus the dynamical bias force is :math:`-P\nabla_b U_B`, while the
reported scalar remains :math:`U_B`.  The PIMD energy estimator accounts
for this distinction whether or not ``fix_modify energy yes`` includes
the scalar in the potential energy.

The *path_integral centroid* setting couples PLUMED to the Cartesian coordinate
centroid provided by the :doc:`fix pimd/langevin <fix_pimd>` command selected
with *pimd_fix*.  The PIMD fix must be defined before fix plumed.  One PLUMED
state is created on partition zero, so there is one bias history rather than an
independent history for each bead.  With *method pimd* and *ensemble nvt*, a
centroid bias force :math:`\mathbf{F}_c` adds :math:`\mathbf{F}_c` to every
one of the :math:`P` Cartesian beads.  With *method nmpimd* and *ensemble nvt*,
*nph*, or *npt*, the centroid coordinate and force obey

.. math::

   \mathbf{q}_0=\frac{1}{\sqrt{P}}\sum_b\mathbf{R}_b,
   \qquad
   \mathbf{F}_{q_0}=\sqrt{P}\mathbf{F}_c,

and the non-centroid modes receive no bias force.  The once-owned bias virial
is included in the current-step centroid-virial pressure.  NPH and NPT use
that pressure for the normal-mode barostat; NVT does not initialize or update
a barostat or the simulation cell.

This mode biases a collective variable evaluated from the coordinate centroid,
which is generally different from averaging the collective variable over the
beads.  The scalar bias energy and bias virial are nonzero only on partition
zero.  The virial uses the dynamical bias normalization, while the scalar
reports the physical bias before the PIMD energy correction.  Bead-resolved
trajectories are still required to reconstruct a bead-defined quantum free
energy.

The *path_integral bead_mean* setting creates one PLUMED instance on every
bead partition and enables PLUMED's multiple-replica communication.  Use the
PLUMED ``ENSEMBLE`` action to define the arithmetic bead mean of a
collective variable and apply biases only to that mean.  For example:

.. code-block:: text

   d: DISTANCE ATOMS=1,2
   mean: ENSEMBLE ARG=d
   bias: RESTRAINT ARG=mean.d AT=0.5 KAPPA=10

If :math:`S=P^{-1}\sum_b s(\mathbf{R}_b)`, PLUMED propagates the bias force
with the chain-rule factor :math:`1/P`.  LAMMPS multiplies only the bias
force increment and its virial by :math:`P` to match the dynamical bias
potential, so bead :math:`b` receives
:math:`-(\partial U_B/\partial S)\nabla_b s`.  Physical forces are unchanged.
Applying a bias directly to
``d`` instead of ``mean.d`` creates independent per-bead biases and is not a
bead-mean calculation.

All PLUMED instances evaluate the same bead-mean bias in lockstep.  PLUMED
adds the partition suffix to its output and restart files.  The LAMMPS scalar
bias energy is reported only on partition zero so that it is counted once;
the chain-rule force and virial contributions remain local to every bead.
Do not use a multiple-walker option to combine the PIMD beads: they are parts
of one ring polymer, not statistically independent walkers.

With *method nmpimd*, the same Cartesian bead coordinates are passed to
PLUMED, and the complete Cartesian force after the PLUMED contribution is
transformed to normal modes.  A nonlinear collective variable can therefore
generate nonzero forces on internal modes; replacing this transformation by a
centroid-only force would be incorrect.  NVT keeps the cell fixed.  NPH and
NPT include the current-step bead-bias virial in the centroid pressure used by
the BZP barostat.

For every *path_integral* mode used with *method nmpimd*, define ``fix
plumed`` after all other fixes that have a post-force callback.  This keeps
Cartesian force contributions ahead of the normal-mode force transformation.
LAMMPS stops with an error if a post-force fix is defined after ``fix
plumed``; fixes without a post-force callback may still follow it.

The same mode can construct the instantaneous path spread without another
LAMMPS communication routines.  For a bead-local scalar :math:`s_b`, define

.. math::

   \sigma_s^2=\frac{1}{P}\sum_b s_b^2-
   \left(\frac{1}{P}\sum_b s_b\right)^2.

For example:

.. code-block:: text

   s: DISTANCE ATOMS=1,2
   s2: CUSTOM ARG=s FUNC=x*x PERIODIC=NO
   mean: ENSEMBLE ARG=s
   mean2: ENSEMBLE ARG=s2
   spread2: CUSTOM ARG=mean.s,mean2.s2 FUNC=y-x*x PERIODIC=NO
   bias: RESTRAINT ARG=spread2 AT=0.0 KAPPA=10

This evaluates :math:`B(\sigma_s^2)` with PLUMED's existing chain-rule
derivatives and preserves one partition-zero scalar bias owner.  For a distance
CV with consistently unwrapped coordinates, a deterministic coupled graph can
also form :math:`s_c=s(\overline{\mathbf R})` from the replica-averaged
Cartesian displacement:

.. code-block:: text

   dvec: DISTANCE ATOMS=1,2 COMPONENTS NOPBC
   s: CUSTOM ARG=dvec.x,dvec.y,dvec.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
   s2: CUSTOM ARG=s FUNC=x*x PERIODIC=NO
   mean: ENSEMBLE ARG=s
   mean2: ENSEMBLE ARG=s2
   spread2: CUSTOM ARG=mean.s,mean2.s2 FUNC=y-x*x PERIODIC=NO
   centroid: ENSEMBLE ARG=dvec.x,dvec.y,dvec.z
   sc: CUSTOM ARG=centroid.dvec.x,centroid.dvec.y,centroid.dvec.z \
       FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
   coupled: CUSTOM ARG=sc,spread2 FUNC=x*y PERIODIC=NO
   bias: BIASVALUE ARG=coupled

This example evaluates the test bias
:math:`B(s_c,\sigma_s^2)=s_c\sigma_s^2` and propagates both centroid and
spread derivatives through the existing ``ENSEMBLE`` graph.  It does not
replace :math:`s_c` with the generally different bead mean
:math:`P^{-1}\sum_b s_b`.  General CVs must be evaluated from the Cartesian
centroid with their own correct periodic-coordinate convention.  A literal
:math:`\sigma_s` coordinate has a singular derivative at zero spread and
requires an explicit regularization; this example does not enable automatic
switching between bias modes.

The *path_integral bead_density* setting implements a symmetric bias of the
instantaneous bead density,

.. math::

   U_B(X,t)=\frac{1}{P}\sum_{b=1}^{P} B(s(\mathbf{R}_b),t).

Each bead evaluates the same PLUMED bias function at its local collective
variable.  LAMMPS retains the local bias force and virial without an
additional :math:`1/P` factor, and
reports the mean of the :math:`P` local bias energies on partition zero.  A
fixed bias therefore needs no replica-averaging action:

.. code-block:: text

   d: DISTANCE ATOMS=1,2
   bias: RESTRAINT ARG=d AT=0.5 KAPPA=10

For *method nmpimd*, the local bias force and virial retain this dynamical
normalization; the physical Cartesian force is left unchanged before the full
force is transformed to normal modes.  NVT keeps the cell fixed, while NPH and
NPT pass the refreshed current-step virial to the BZP pressure path.

For a history-dependent bias, every PLUMED instance must share one field.
For example, ``METAD`` can use ``WALKERS_MPI`` as the field-communication
mechanism.  Because all :math:`P` beads deposit at the same physical time, the
hill height on each bead must be the intended shared height divided by
:math:`P`.  The beads remain correlated coordinates of one ring polymer; the
multiple-walker machinery does not make them independent statistical samples.
Using separate histories, or using an unscaled hill height on every bead,
does not implement this mode's shared bead-density Hamiltonian.
LAMMPS does not inspect arbitrary PLUMED inputs to verify field sharing or
deposition normalization.  The native regressions cover fixed biases, a
matched five-step centroid/bead-density linear-bias dynamics limit,
``METAD WALKERS_MPI`` with single- and multi-rank bead partitions, and a
four-bead ``OPES_METAD WALKERS_MPI`` restart.  They validate one shared HILLS
stream with :math:`1/P` metadynamics hill heights, shared OPES KERNELS and STATE files,
partition-zero bias ownership, zero-local-atom ranks, and PLUMED file-restart
continuity.  This is an interface contract, not production admission for OPES
reweighting, binary-restart dynamics, performance, or sampling efficiency.

The *bead_density* and *bead_mean* modes are different.  The former evaluates
:math:`P^{-1}\sum_b B(s_b)`, while the latter evaluates
:math:`B(P^{-1}\sum_b s_b)` through ``ENSEMBLE``.  They are generally unequal
for nonlinear collective variables or nonlinear biases.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When performing a restart of a calculation that involves PLUMED you
must include a RESTART command in the PLUMED input file as detailed in
the `PLUMED documentation <plumeddocs_>`_.  When the restart command
is found in the PLUMED input PLUMED will append to the files that were
generated in the run that was performed previously.  No part of the
PLUMED restart data is included in the LAMMPS restart files.
Furthermore, any history dependent bias potentials that were
accumulated in previous calculations will be read in when the RESTART
command is included in the PLUMED input.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by
this fix to add the energy change from the biasing force added by
PLUMED to the global potential energy of the system as part of
:doc:`thermodynamic output <thermo_style>`.  The default setting for
this fix is :doc:`fix_modify energy yes <fix_modify>`.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by
this fix to add the contribution from the biasing force to the global
pressure of the system via the :doc:`compute pressure
<compute_pressure>` command.  This can be accessed by
:doc:`thermodynamic output <thermo_style>`.  The default setting for
this fix is :doc:`fix_modify virial yes <fix_modify>`.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the PLUMED
energy mentioned above.  The scalar value calculated by this fix is
"extensive".

Note that other quantities of interest can be output by commands that
are native to PLUMED.

Fixed conditional path functions
--------------------------------

A fixed complete-path function combining a Cartesian-centroid CV and a
bead-averaged score can use the same *bead_mean* adapter when both inputs
are constructed correctly in the PLUMED graph.  For example, average
Cartesian components first with ``ENSEMBLE`` and then apply a nonlinear
function to form the centroid CV.  Average the bead-local score separately.
Do not replace a nonlinear centroid CV with the mean of that nonlinear CV.
The component construction requires consistent coordinate images and does
not provide an automatic Cartesian-centroid interface for arbitrary CVs.

For positive score :math:`a(X)`, a frozen positive normalizer :math:`m(c)`,
and :math:`0\le\lambda<1`, a conditional correction may have the form

.. math::

   U_B(X)=B_c(c)-k_B T\log[(1-\lambda)+\lambda a(X)/m(c)].

Both the score and normalizer derivatives must remain in the graph.
The optional PLUMED ``CONDITIONAL_PATH`` function evaluates the log mixture;
ordinary ``CUSTOM`` and ``BIASVALUE`` actions can apply its energy correction.
This does not introduce another *path_integral* mode.  The adapter retains
its existing physical-energy and dynamical-force normalization; do not
multiply this function by another bead-count factor.  Keep all fields fixed
and verify input/model identities before a restart.  Equilibrium reweighting
uses one total-bias weight per complete path, not one independent weight per
bead.  The existence of this force graph does not establish sampling gains.

Frozen probability-ratio mixtures
---------------------------------

A common frozen field :math:`v(s)` can instead define the total path bias

.. math::

   U_A(X)=-k_B T\log\left[\frac{1}{P}\sum_b\exp(-v(s_b)/k_B T)\right].

The optional PLUMED ``PATH_LOGMEANEXP`` action provides a numerically stable
replica reduction and its bead-local softmax derivative. Use *bead_mean*
for this complete-path graph. Its physical force coefficient already
contains the normalization; another :math:`1/P` factor is incorrect.
``EXPECTED_REPLICAS`` must equal the number of bead partitions, not the
number of spatial MPI ranks. The native normalization tests include pure,
centroid-only and mixed frozen graphs using ordinary PLUMED functions.

For a frozen active OPES action, apply :math:`U_A-v(s_b)` as a correction
so that the original local bias force and energy are canceled. Applying
both the full :math:`U_A` and the local OPES bias would double count the
field. Verify immutable state identity and native update suppression;
current shared-density OPES deposition weights do not implement adaptive
learning under this new Hamiltonian. ``PROBABILITY_MIX`` can combine a
non-bias centroid scalar and :math:`U_A` using a frozen global normalizer.
The initial path-mixture qualification is fixed-volume NVT; this does not
admit pressure-coupled or arbitrary molecular-centroid variants.

Restrictions
""""""""""""

This fix is part of the PLUMED package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

There can only be one fix plumed command active at a time.

The *centroid*, *bead_mean*, and *bead_density* modes require
:doc:`fix pimd/langevin <fix_pimd>` with either *method pimd* and *ensemble
nvt*, or *method nmpimd* and *ensemble nvt*, *nph*, or *npt*.  Normal-mode
pressure ensembles use the BZP barostat supported by fix pimd/langevin;
normal-mode NVE and Cartesian-PIMD pressure coupling are not supported.  All
three modes require a fixed atom count, consecutive atom IDs, and an atom map,
and can distribute each bead over multiple MPI ranks.  The *bead_mean* and
*bead_density* modes additionally require multiple LAMMPS partitions.  A
*bead_mean* input must explicitly form a complete-path bias with
``ENSEMBLE`` or another differentiated replica reduction such as
``PATH_LOGMEANEXP``.  A
history-dependent *bead_density* input must use one shared bias field and scale
each bead's deposition by :math:`1/P`.  None of the path-integral modes
supports energy-dependent PLUMED actions, minimization, or r-RESPA.
The default *path_integral off* setting remains incompatible with path-integral
fixes.

Related commands
""""""""""""""""

:doc:`fix smd <fix_smd>`
:doc:`fix colvars <fix_colvars>`

Default
"""""""

The default options are plumedfile = NULL, outfile = NULL, and path_integral = off.

----------

.. _Tribello:

**(Tribello)** G.A. Tribello, M. Bonomi, D. Branduardi, C. Camilloni and G. Bussi, Comp. Phys. Comm 185, 604 (2014)

.. _plumeddocs: https://www.plumed.org/doc.html

.. _plumedhome: https://www.plumed.org/
