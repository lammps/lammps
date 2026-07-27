Tests for granular (DEM) models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionadded:: 4Jul2026

The unittest/granular folder contains a YAML-driven test suite for
discrete element method (DEM) / granular contact models, built in the same
spirit as the force-style tests above but specialized for time-resolved
trajectories of small, analytically tractable granular systems. These tests
are only enabled if the :ref:GRANULAR package<PKG-GRANULAR> is enabled.

There are 8 test programs, test_dem_01 through test_dem_08, covering
particle-impact-level benchmarks: two-sphere and sphere-wall collisions,
oblique and spinning-sphere impacts, rolling and slipping contact, and
cohesive pull-off. These follow the software-agnostic DEM benchmark of
:ref:Mohajeri et al. <dem_Mohajeri2024> (rolling resistance and cohesion)
and the particle-impact benchmark of :ref:Chung and Ooi <dem_Chung2011>
(normal, oblique, and spinning-sphere collisions). The test programs are:

.. list-table::
:header-rows: 1

Program
Scenario
Analytic model
test_dem_01
two-sphere head-on normal collision
collision_restitution
test_dem_02
two-sphere elastic Hertzian normal impact (peak force)
hertz_normal_impact
test_dem_03
sphere-wall elastic Hertzian normal impact (peak force)
hertz_normal_impact
test_dem_04
oblique impact on a wall (gross-sliding friction)
oblique_impact
test_dem_05
sphere sliding then rolling without slipping
slip_cessation
test_dem_06
spinning sphere impact (rebound with friction)
spin_impact
test_dem_07
rolling-resistance decay
rolling_decay
test_dem_08
cohesive DMT pull-off force
pulloff_dmt

Every test program shares the same driver logic, implemented in
unittest/granular/test_dem_common.cpp and compiled into the
granular_tests support library; each test_dem_0N.cpp only contains
the two GoogleTest fixtures (newton_on and newton_off). As with the
force-style tests, the reference systems are defined by YAML files in the
unittest/granular/tests folder and registered as CTest cases by their
file name (dem0N-*.yaml becomes test DEM0N:*); adding or removing a
YAML file requires re-running CMake. A given driver may cover several
variants of one scenario -- across contact models
(hooke, hooke/history, hertz, hertz/material, mindlin,
mindlin/rescale), dimensionality, or unit systems -- since it simply
runs every dem0N-*.yaml file that matches its number.

Unlike the force-style tests, the entire system is built from the YAML
file rather than from a fixed input template. A YAML file provides an
optional variables block (emitted as :doc:index variables <variable>
so they can be substituted as ${name} anywhere in the command strings),
pre_commands that create the geometry, pair_style / pair_coeff
that select the contact model, and post_commands that add the
integrator and any walls. The trajectory is then advanced in a sequence of
run_segments and, after each segment, the per-atom positions,
velocities, torques, and angular velocities are compared against the
recorded reference. A minimal example (dem03-hertz-wall-3d-si.yaml)
looks like:

.. code-block:: yaml

lammps_version: 4 Jul 2026
tags: granular, generated
date_generated: Tue Jul 21 21:40:51 2026
epsilon: 1e-10
prerequisites: |
atom sphere
pair granular
pre_commands: |
units si
dimension 3
boundary f f f
atom_style sphere
region box block -0.05 0.05 -0.05 0.05 -0.05 0.1 units box
create_box 1 box
create_atoms 1 single 0.0 0.0 ${z0} units box
set group all diameter ${diam} density ${dens}
velocity all set 0.0 0.0 -${vin}
comm_modify vel yes
neighbor 0.001 bin
neigh_modify delay 0
timestep ${dt}
post_commands: |
fix integr all nve/sphere
fix zwall all wall/gran granular hertz ${kn} ${en} tangential linear_nohistory 0.0 0.0 damping coeff_restitution zplane 0.0 NULL
natoms: 1
variables: |
diam 0.005
dens 7000.0
kn 7.11111e+10
en 1.0
vin 3.9
vrela 3.9
mred_factor 1.0
radius 0.0025
dt 1.0e-8
z0 0.002550
pair_style: granular
pair_coeff: |
* * hertz ${kn} ${en} tangential linear_nohistory 0.0 0.0 damping coeff_restitution
run_segments: |-
1300 838 900
analytic_enable: yes
analytic_model: hertz_normal_impact
analytic_tol: 0.01
analytic_segment: 1

run_pos / run_vel / run_torque / run_omega blocks follow

The following table describes the available keys:

.. list-table::
:header-rows: 1

Key:
Description:
epsilon
relative precision required for the recorded (regression) reference data
prerequisites
list of style kind / style name pairs required to run the test
variables
name/value pairs exposed as ${name} index variables for substitution
pre_commands
commands that build the geometry (units, box, atoms, set, timestep)
pair_style / pair_coeff
the particle-particle (or particle-wall) contact model
post_commands
fixes added after the geometry (integrator, walls)
run_segments
whitespace-separated list of run lengths; state is captured after each
run_pos, run_vel
reference positions and velocities, as segment tag x y z rows
run_torque, run_omega
reference torque and angular velocity (when applicable)
analytic_enable
yes to also assert a closed-form (analytic) model
analytic_model
which analytic model to evaluate (see below)
analytic_tol
relative tolerance for the analytic assertion (looser than epsilon)
analytic_segment
run segment at which the analytic model is checked (-1 means the last)

The per-atom reference blocks use a segment tag x y z row format, so a
single block holds the data for all run segments and the row order does not
matter. Because granular/atomic systems do not build an atom map by
default, the reference generator iterates over local atoms by tag rather
than calling Atom::map().

Each test runs as a pure regression check (the recorded data is reproduced
to within epsilon) under both the newton on and newton off
fixtures, which are expected to give identical results. In addition, every
test opts in to an analytic check that compares a derived quantity against
a closed-form solution implemented in
unittest/granular/test_analytic_models.cpp. The analytic tolerance is
deliberately loose, because the soft-sphere DEM result only approaches the
idealized (hard-sphere or instantaneous-contact) solution. The models
currently implemented are:

.. list-table::
:header-rows: 1

Model:
Checks:
collision_restitution
two-sphere momentum conservation and restitution :math:e = -(v_1'-v_2')/(v_1-v_2)
hertz_normal_impact
Hertzian peak energy balance :math:\tfrac{1}{2}\mu_{red} V_{rela}^2 = \tfrac{2}{5} P_{max}\alpha_{max}
oblique_impact
gross-sliding rebound :math:v_x' = v_x - \mu(1+e)v_z, :math:\omega_y = \tfrac{5}{2}\mu(1+e)v_z/r
slip_cessation
rolling-without-slipping limit :math:u = 5 u_0/7, :math:\omega = u/r
spin_impact
gross-sliding rebound of a spinning sphere: :math:v_x' = \mu(1+e)v_n, :math:\omega_y' = \omega_0 - \tfrac{5}{2}\mu(1+e)v_n/r
rolling_decay
linear spin-down under rolling resistance: :math:\omega = \omega_0 - \tfrac{5 \mu_r g}{2 r} t
pulloff_dmt
DMT pull-off force at contact :math:|F| = 4 \pi \gamma R_{\mathrm{eff}}

Adding a new reference (YAML) file
""""""""""""""""""""""""""""""""""

Copy an existing dem0N-*.yaml
for a similar scenario, adjust the variables, pre_commands,
pair_style/pair_coeff and post_commands for the new model, and
give it a new name matching the dem0N-*.yaml pattern of the test program
it belongs to. Leave out the reference data blocks initially, then
(re)generate them in place with:

.. code-block:: bash

TEST_ARGS=-u ctest -R DEM0N:myvariant

or by running the driver directly (test_dem_0N dem0N-myvariant.yaml -u).
Do not write the generated file to a sibling dem0N-*.yaml name (for
example with the -g newfile.yaml option pointing into the tests
folder), because the CONFIGURE_DEPENDS glob would then register it as an
extra, stale test. After adding the file, re-run CMake so the new test is
registered, then verify it with ctest -V -R DEM0N:myvariant (the -s
option of the driver reports per-quantity error statistics, which helps when
choosing epsilon and the analytic tolerance).

Adding a new test program
"""""""""""""""""""""""""

Create test_dem_0N.cpp as a thin copy of an existing one (only the
GoogleTest suite name changes), add an
add_executable/register_dem_tests pair to
unittest/granular/CMakeLists.txt, and add the corresponding
dem0N-*.yaml reference files. If the new scenario needs a
closed-form check, add a named model to test_analytic_models.cpp
that reads its parameters from the variables block (and reads
masses, radii, etc. from the live LAMMPS instance to avoid depending on
derived quantities) and assert it with EXPECT_LE on the relative
error.

References
""""""""""

.. _dem_Mohajeri2024:

(Mohajeri et al., 2024) M. J. Mohajeri, C. Coetzee, and D. L. Schott,
A software-agnostic benchmark for DEM simulation of cohesive and
non-cohesive materials, Powder Technology, 447, 120136 (2024),
https://doi.org/10.1016/j.powtec.2024.120136

.. _dem_Chung2011:

(Chung and Ooi, 2011) Y. C. Chung and J. Y. Ooi, Benchmark tests for
verifying discrete element modelling codes at particle impact level,
Granular Matter, 13, 643-656 (2011),
https://doi.org/10.1007/s10035-011-0277-0
