Tutorial for Thermalized Drude oscillators in LAMMPS
====================================================

This tutorial explains how to use Drude oscillators in LAMMPS to
simulate polarizable systems using the USER-DRUDE package. As an
illustration, the input files for a simulation of 250 phenol molecules
are documented. First of all, LAMMPS has to be compiled with the
USER-DRUDE package activated. Then, the data file and input scripts
have to be modified to include the Drude dipoles and how to handle
them.

----------

**Overview of Drude induced dipoles**

Polarizable atoms acquire an induced electric dipole moment under the
action of an external electric field, for example the electric field
created by the surrounding particles.  Drude oscillators represent
these dipoles by two fixed charges: the core (DC) and the Drude
particle (DP) bound by a harmonic potential. The Drude particle can be
thought of as the electron cloud whose center can be displaced from
the position of the corresponding nucleus.

The sum of the masses of a core-Drude pair should be the mass of the
initial (unsplit) atom, :math:`m_C + m_D = m`.  The sum of their charges
should be the charge of the initial (unsplit) atom, :math:`q_C + q_D = q`.
A harmonic potential between the core and Drude partners should be
present, with force constant :math:`k_D` and an equilibrium distance of
zero. The (half-)stiffness of the :doc:`harmonic bond <bond_harmonic>`
:math:`K_D = k_D/2` and the Drude charge :math:`q_D` are related to the atom
polarizability :math:`\alpha` by

.. math::

   K_D = \frac 1 2\, \frac {q_D^2} \alpha

Ideally, the mass of the Drude particle should be small, and the
stiffness of the harmonic bond should be large, so that the Drude
particle remains close ot the core. The values of Drude mass, Drude
charge, and force constant can be chosen following different
strategies, as in the following examples of polarizable force
fields:

* :ref:`Lamoureux and Roux <Lamoureux2>` suggest adopting a global half-stiffness, :math:`K_D` = 500 kcal/(mol Ang :math:`{}^2`) - which corresponds to a force constant :math:`k_D` = 4184 kJ/(mol Ang :math:`{}^2`) - for all types of core-Drude bond, a global mass :math:`m_D` = 0.4 g/mol (or u) for all types of Drude particles, and to calculate the Drude charges for individual atom types from the atom polarizabilities using equation (1). This choice is followed in the polarizable CHARMM force field.
* Alternately :ref:`Schroeder and Steinhauser <Schroeder>` suggest adopting a global charge :math:`q_D` = -1.0e and a global mass :math:`m_D` = 0.1 g/mol (or u) for all Drude particles, and to calculate the force constant for each type of core-Drude bond from equation (1). The timesteps used by these authors are between 0.5 and 2 fs, with the degrees of freedom of the Drude oscillators kept cold at 1 K.
* In both these force fields hydrogen atoms are treated as non-polarizable.

The motion of of the Drude particles can be calculated by minimizing
the energy of the induced dipoles at each timestep, by an iterative,
self-consistent procedure. The Drude particles can be massless and
therefore do not contribute to the kinetic energy. However, the
relaxed method is computational slow. An extended-lagrangian method
can be used to calculate the positions of the Drude particles, but
this requires them to have mass. It is important in this case to
decouple the degrees of freedom associated with the Drude oscillators
from those of the normal atoms. Thermalizing the Drude dipoles at
temperatures comparable to the rest of the simulation leads to several
problems (kinetic energy transfer, very short timestep, etc.), which
can be remedied by the "cold Drude" technique (:ref:`Lamoureux and Roux <Lamoureux2>`).

Two closely related models are used to represent polarization through
"charges on a spring": the core-shell model and the Drude
model. Although the basic idea is the same, the core-shell model is
normally used for ionic/crystalline materials, whereas the Drude model
is normally used for molecular systems and fluid states. In ionic
crystals the symmetry around each ion and the distance between them
are such that the core-shell model is sufficiently stable. But to be
applicable to molecular/covalent systems the Drude model includes two
important features:

#. The possibility to thermostat the additional degrees of freedom associated with the induced dipoles at very low temperature, in terms of the reduced coordinates of the Drude particles with respect to their cores. This makes the trajectory close to that of relaxed induced dipoles.
#. The Drude dipoles on covalently bonded atoms interact too strongly due to the short distances, so an atom may capture the Drude particle (shell) of a neighbor, or the induced dipoles within the same molecule may align too much.  To avoid this, damping at short of the interactions between the point charges composing the induced dipole can be done by :ref:`Thole <Thole2>` functions.

----------

**Preparation of the data file**

The data file is similar to a standard LAMMPS data file for
*atom\_style full*.  The DPs and the *harmonic bonds* connecting them
to their DC should appear in the data file as normal atoms and bonds.

You can use the *polarizer* tool (Python script distributed with the
USER-DRUDE package) to convert a non-polarizable data file (here
*data.102494.lmp*\ ) to a polarizable data file (\ *data-p.lmp*\ )

.. parsed-literal::

   polarizer -q -f phenol.dff data.102494.lmp data-p.lmp

This will automatically insert the new atoms and bonds.
The masses and charges of DCs and DPs are computed
from *phenol.dff*\ , as well as the DC-DP bond constants.  The file
*phenol.dff* contains the polarizabilities of the atom types
and the mass of the Drude particles, for instance:

.. parsed-literal::

   # units: kJ/mol, A, deg
   # kforce is in the form k/2 r_D\^2
   # type  m_D/u   q_D/e    k_D   alpha/A3  thole
   OH      0.4    -1.0    4184.0   0.63     0.67
   CA      0.4    -1.0    4184.0   1.36     2.51
   CAI     0.4    -1.0    4184.0   1.09     2.51

The hydrogen atoms are absent from this file, so they will be treated
as non-polarizable atoms.  In the non-polarizable data file
*data.102494.lmp*\ , atom names corresponding to the atom type numbers
have to be specified as comments at the end of lines of the *Masses*
section.  You probably need to edit it to add these names. It should
look like

.. parsed-literal::

   Masses

   1 12.011 # CAI
   2 12.011 # CA
   3 15.999 # OH
   4 1.008  # HA
   5 1.008  # HO

----------

**Basic input file**

The atom style should be set to (or derive from) *full*\ , so that you
can define atomic charges and molecular bonds, angles, dihedrals...

The *polarizer* tool also outputs certain lines related to the input
script (the use of these lines will be explained below).  In order for
LAMMPS to recognize that you are using Drude oscillators, you should
use the fix *drude*\ . The command is

.. code-block:: LAMMPS

   fix DRUDE all drude C C C N N D D D

The N, C, D following the *drude* keyword have the following meaning:
There is one tag for each atom type. This tag is C for DCs, D for DPs
and N for non-polarizable atoms.  Here the atom types 1 to 3 (C and O
atoms) are DC, atom types 4 and 5 (H atoms) are non-polarizable and
the atom types 6 to 8 are the newly created DPs.

By recognizing the fix *drude*\ , LAMMPS will find and store matching
DC-DP pairs and will treat DP as equivalent to their DC in the
*special bonds* relations.  It may be necessary to extend the space
for storing such special relations.  In this case extra space should
be reserved by using the *extra/special/per/atom* keyword of either
the :doc:`read_data <read_data>` or :doc:`create_box <create_box>`
command.  With our phenol, there is 1 more special neighbor for which
space is required.  Otherwise LAMMPS crashes and gives the required
value.

.. code-block:: LAMMPS

   read_data data-p.lmp extra/special/per/atom 1

Let us assume we want to run a simple NVT simulation at 300 K.  Note
that Drude oscillators need to be thermalized at a low temperature in
order to approximate a self-consistent field (SCF), therefore it is not
possible to simulate an NVE ensemble with this package.  Since dipoles
are approximated by a charged DC-DP pair, the *pair\_style* must
include Coulomb interactions, for instance *lj/cut/coul/long* with
*kspace\_style pppm*. For example, with a cutoff of 10. and a precision
1.e-4:

.. code-block:: LAMMPS

   pair_style lj/cut/coul/long 10.0
   kspace_style pppm 1.0e-4

As compared to the non-polarizable input file, *pair\_coeff* lines need
to be added for the DPs.  Since the DPs have no Lennard-Jones
interactions, their :math:`\epsilon` is 0. so the only *pair\_coeff* line
that needs to be added is

.. code-block:: LAMMPS

   pair_coeff * 6* 0.0 0.0 # All-DPs

Now for the thermalization, the simplest choice is to use the :doc:`fix langevin/drude <fix_langevin_drude>`.

.. code-block:: LAMMPS

   fix LANG all langevin/drude 300. 100 12435 1. 20 13977

This applies a Langevin thermostat at temperature 300. to the centers
of mass of the DC-DP pairs, with relaxation time 100 and with random
seed 12345.  This fix applies also a Langevin thermostat at temperature
1. to the relative motion of the DPs around their DCs, with relaxation
time 20 and random seed 13977.  Only the DCs and non-polarizable
atoms need to be in this fix's group.  LAMMPS will thermostat the DPs
together with their DC.  For this, ghost atoms need to know their
velocities. Thus you need to add the following command:

.. code-block:: LAMMPS

   comm_modify vel yes

In order to avoid that the center of mass of the whole system
drifts due to the random forces of the Langevin thermostat on DCs, you
can add the *zero yes* option at the end of the fix line.

If the fix *shake* is used to constrain the C-H bonds, it should be
invoked after the fix *langevin/drude* for more accuracy.

.. code-block:: LAMMPS

   fix SHAKE ATOMS shake 0.0001 20 0 t 4 5

.. note::

   The group of the fix *shake* must not include the DPs.  If the
   group *ATOMS* is defined by non-DPs atom types, you could use

Since the fix *langevin/drude* does not perform time integration (just
modification of forces but no position/velocity updates), the fix
*nve* should be used in conjunction.

.. code-block:: LAMMPS

   fix NVE all nve

Finally, do not forget to update the atom type elements if you use
them in a *dump\_modify ... element ...* command, by adding the element
type of the DPs. Here for instance

.. code-block:: LAMMPS

   dump DUMP all custom 10 dump.lammpstrj id mol type element x y z ix iy iz
   dump_modify DUMP element C C O H H D D D

The input file should now be ready for use!

You will notice that the global temperature *thermo\_temp* computed by
LAMMPS is not 300. K as wanted.  This is because LAMMPS treats DPs as
standard atoms in his default compute.  If you want to output the
temperatures of the DC-DP pair centers of mass and of the DPs relative
to their DCs, you should use the :doc:`compute temp\_drude <compute_temp_drude>`

.. code-block:: LAMMPS

   compute TDRUDE all temp/drude

And then output the correct temperatures of the Drude oscillators
using *thermo\_style custom* with respectively *c\_TDRUDE[1]* and
*c\_TDRUDE[2]*. These should be close to 300.0 and 1.0 on average.

.. code-block:: LAMMPS

   thermo_style custom step temp c_TDRUDE[1] c_TDRUDE[2]

----------

**Thole screening**

Dipolar interactions represented by point charges on springs may not
be stable, for example if the atomic polarizability is too high for
instance, a DP can escape from its DC and be captured by another DC,
which makes the force and energy diverge and the simulation
crash. Even without reaching this extreme case, the correlation
between nearby dipoles on the same molecule may be exaggerated.  Often,
special bond relations prevent bonded neighboring atoms to see the
charge of each other's DP, so that the problem does not always appear.
It is possible to use screened dipole-dipole interactions by using the
:doc:`*pair\_style thole* <pair_thole>`.  This is implemented as a
correction to the Coulomb pair\_styles, which dampens at short distance
the interactions between the charges representing the induced dipoles.
It is to be used as *hybrid/overlay* with any standard *coul* pair
style.  In our example, we would use

.. code-block:: LAMMPS

   pair_style hybrid/overlay lj/cut/coul/long 10.0 thole 2.6 10.0

This tells LAMMPS that we are using two pair\_styles.  The first one is
as above (\ *lj/cut/coul/long 10.0*\ ).  The second one is a *thole*
pair\_style with default screening factor 2.6 (:ref:`Noskov <Noskov2>`) and
cutoff 10.0.

Since *hybrid/overlay* does not support mixing rules, the interaction
coefficients of all the pairs of atom types with i < j should be
explicitly defined.  The output of the *polarizer* script can be used
to complete the *pair\_coeff* section of the input file.  In our
example, this will look like:

.. code-block:: LAMMPS

   pair_coeff    1    1 lj/cut/coul/long    0.0700   3.550
   pair_coeff    1    2 lj/cut/coul/long    0.0700   3.550
   pair_coeff    1    3 lj/cut/coul/long    0.1091   3.310
   pair_coeff    1    4 lj/cut/coul/long    0.0458   2.985
   pair_coeff    2    2 lj/cut/coul/long    0.0700   3.550
   pair_coeff    2    3 lj/cut/coul/long    0.1091   3.310
   pair_coeff    2    4 lj/cut/coul/long    0.0458   2.985
   pair_coeff    3    3 lj/cut/coul/long    0.1700   3.070
   pair_coeff    3    4 lj/cut/coul/long    0.0714   2.745
   pair_coeff    4    4 lj/cut/coul/long    0.0300   2.420
   pair_coeff    *    5 lj/cut/coul/long    0.0000   0.000
   pair_coeff    *   6* lj/cut/coul/long    0.0000   0.000
   pair_coeff    1    1 thole   1.090   2.510
   pair_coeff    1    2 thole   1.218   2.510
   pair_coeff    1    3 thole   0.829   1.590
   pair_coeff    1    6 thole   1.090   2.510
   pair_coeff    1    7 thole   1.218   2.510
   pair_coeff    1    8 thole   0.829   1.590
   pair_coeff    2    2 thole   1.360   2.510
   pair_coeff    2    3 thole   0.926   1.590
   pair_coeff    2    6 thole   1.218   2.510
   pair_coeff    2    7 thole   1.360   2.510
   pair_coeff    2    8 thole   0.926   1.590
   pair_coeff    3    3 thole   0.630   0.670
   pair_coeff    3    6 thole   0.829   1.590
   pair_coeff    3    7 thole   0.926   1.590
   pair_coeff    3    8 thole   0.630   0.670
   pair_coeff    6    6 thole   1.090   2.510
   pair_coeff    6    7 thole   1.218   2.510
   pair_coeff    6    8 thole   0.829   1.590
   pair_coeff    7    7 thole   1.360   2.510
   pair_coeff    7    8 thole   0.926   1.590
   pair_coeff    8    8 thole   0.630   0.670

For the *thole* pair style the coefficients are

#. the atom polarizability in units of cubic length
#. the screening factor of the Thole function (optional, default value
   specified by the pair\_style command)
#. the cutoff (optional, default value defined by the pair\_style command)

The special neighbors have charge-charge and charge-dipole
interactions screened by the *coul* factors of the *special\_bonds*
command (0.0, 0.0, and 0.5 in the example above).  Without using the
pair\_style *thole*\ , dipole-dipole interactions are screened by the
same factor.  By using the pair\_style *thole*\ , dipole-dipole
interactions are screened by Thole's function, whatever their special
relationship (except within each DC-DP pair of course).  Consider for
example 1-2 neighbors: using the pair\_style *thole*\ , their dipoles
will see each other (despite the *coul* factor being 0.) and the
interactions between these dipoles will be damped by Thole's function.

----------

**Thermostats and barostats**

Using a Nose-Hoover barostat with the *langevin/drude* thermostat is
straightforward using fix *nph* instead of *nve*\ .  For example:

.. code-block:: LAMMPS

   fix NPH all nph iso 1. 1. 500

It is also possible to use a Nose-Hoover instead of a Langevin
thermostat.  This requires to use :doc:`\ *fix drude/transform*\ <fix_drude_transform>` just before and after the
time integration fixes.  The *fix drude/transform/direct* converts the
atomic masses, positions, velocities and forces into a reduced
representation, where the DCs transform into the centers of mass of
the DC-DP pairs and the DPs transform into their relative position
with respect to their DC. The *fix drude/transform/inverse* performs
the reverse transformation.  For a NVT simulation, with the DCs and
atoms at 300 K and the DPs at 1 K relative to their DC one would use

.. code-block:: LAMMPS

   fix DIRECT all drude/transform/direct
   fix NVT1 ATOMS nvt temp 300. 300. 100
   fix NVT2 DRUDES nvt temp 1. 1. 20
   fix INVERSE all drude/transform/inverse

For our phenol example, the groups would be defined as

.. code-block:: LAMMPS

   group ATOMS  type 1 2 3 4 5 # DCs and non-polarizable atoms
   group CORES  type 1 2 3     # DCs
   group DRUDES type 6 7 8     # DPs

Note that with the fixes *drude/transform*\ , it is not required to
specify *comm\_modify vel yes* because the fixes do it anyway (several
times and for the forces also).  To avoid the flying ice cube artifact
:ref:`(Lamoureux) <Lamoureux2>`, where the atoms progressively freeze and the
center of mass of the whole system drifts faster and faster, the *fix
momentum* can be used. For instance:

.. code-block:: LAMMPS

   fix MOMENTUM all momentum 100 linear 1 1 1

It is a bit more tricky to run a NPT simulation with Nose-Hoover
barostat and thermostat.  First, the volume should be integrated only
once. So the fix for DCs and atoms should be *npt* while the fix for
DPs should be *nvt* (or vice versa).  Second, the *fix npt* computes a
global pressure and thus a global temperature whatever the fix group.
We do want the pressure to correspond to the whole system, but we want
the temperature to correspond to the fix group only.  We must then use
the *fix\_modify* command for this.  In the end, the block of
instructions for thermostatting and barostatting will look like

.. code-block:: LAMMPS

   compute TATOMS ATOMS temp
   fix DIRECT all drude/transform/direct
   fix NPT ATOMS npt temp 300. 300. 100 iso 1. 1. 500
   fix_modify NPT temp TATOMS press thermo_press
   fix NVT DRUDES nvt temp 1. 1. 20
   fix INVERSE all drude/transform/inverse

----------

**Rigid bodies**

You may want to simulate molecules as rigid bodies (but polarizable).
Common cases are water models such as :ref:`SWM4-NDP <SWM4-NDP>`, which is a
kind of polarizable TIP4P water.  The rigid bodies and the DPs should
be integrated separately, even with the Langevin thermostat.  Let us
review the different thermostats and ensemble combinations.

NVT ensemble using Langevin thermostat:

.. code-block:: LAMMPS

   comm_modify vel yes
   fix LANG all langevin/drude 300. 100 12435 1. 20 13977
   fix RIGID ATOMS rigid/nve/small molecule
   fix NVE DRUDES nve

NVT ensemble using Nose-Hoover thermostat:

.. code-block:: LAMMPS

   fix DIRECT all drude/transform/direct
   fix RIGID ATOMS rigid/nvt/small molecule temp 300. 300. 100
   fix NVT DRUDES nvt temp 1. 1. 20
   fix INVERSE all drude/transform/inverse

NPT ensemble with Langevin thermostat:

.. code-block:: LAMMPS

   comm_modify vel yes
   fix LANG all langevin/drude 300. 100 12435 1. 20 13977
   fix RIGID ATOMS rigid/nph/small molecule iso 1. 1. 500
   fix NVE DRUDES nve

NPT ensemble using Nose-Hoover thermostat:

.. code-block:: LAMMPS

   compute TATOM ATOMS temp
   fix DIRECT all drude/transform/direct
   fix RIGID ATOMS rigid/npt/small molecule temp 300. 300. 100 iso 1. 1. 500
   fix_modify RIGID temp TATOM press thermo_press
   fix NVT DRUDES nvt temp 1. 1. 20
   fix INVERSE all drude/transform/inverse

----------

.. _Lamoureux2:

**(Lamoureux)** Lamoureux and Roux, J Chem Phys, 119, 3025-3039 (2003)

.. _Schroeder:

**(Schroeder)**  Schroeder and Steinhauser, J Chem Phys, 133,
154511 (2010).

.. _Jiang2:

**(Jiang)** Jiang, Hardy, Phillips, MacKerell, Schulten, and Roux,
 J Phys Chem Lett, 2, 87-92 (2011).

.. _Thole2:

**(Thole)** Chem Phys, 59, 341 (1981).

.. _Noskov2:

**(Noskov)** Noskov, Lamoureux and Roux, J Phys Chem B, 109, 6705 (2005).

.. _SWM4-NDP:

**(SWM4-NDP)** Lamoureux, Harder, Vorobyov, Roux, MacKerell, Chem Phys
Let, 418, 245-249 (2006)
