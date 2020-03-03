Adiabatic core/shell model
==========================

The adiabatic core-shell model by :ref:`Mitchell and Fincham <MitchellFincham>` is a simple method for adding polarizability
to a system.  In order to mimic the electron shell of an ion, a
satellite particle is attached to it. This way the ions are split into
a core and a shell where the latter is meant to react to the
electrostatic environment inducing polarizability.  See the :doc:`Howto polarizable <Howto_polarizable>` doc page for a discussion of all
the polarizable models available in LAMMPS.

Technically, shells are attached to the cores by a spring force f =
k\*r where k is a parameterized spring constant and r is the distance
between the core and the shell. The charges of the core and the shell
add up to the ion charge, thus q(ion) = q(core) + q(shell). This
setup introduces the ion polarizability (alpha) given by
alpha = q(shell)\^2 / k. In a
similar fashion the mass of the ion is distributed on the core and the
shell with the core having the larger mass.

To run this model in LAMMPS, :doc:`atom_style <atom_style>` *full* can
be used since atom charge and bonds are needed.  Each kind of
core/shell pair requires two atom types and a bond type.  The core and
shell of a core/shell pair should be bonded to each other with a
harmonic bond that provides the spring force. For example, a data file
for NaCl, as found in examples/coreshell, has this format:


.. parsed-literal::

   432   atoms  # core and shell atoms
   216   bonds  # number of core/shell springs

   4     atom types  # 2 cores and 2 shells for Na and Cl
   2     bond types

   0.0 24.09597 xlo xhi
   0.0 24.09597 ylo yhi
   0.0 24.09597 zlo zhi

   Masses       # core/shell mass ratio = 0.1

   1 20.690784  # Na core
   2 31.90500   # Cl core
   3 2.298976   # Na shell
   4 3.54500    # Cl shell

   Atoms

   1    1    2   1.5005    0.00000000   0.00000000   0.00000000 # core of core/shell pair 1
   2    1    4  -2.5005    0.00000000   0.00000000   0.00000000 # shell of core/shell pair 1
   3    2    1   1.5056    4.01599500   4.01599500   4.01599500 # core of core/shell pair 2
   4    2    3  -0.5056    4.01599500   4.01599500   4.01599500 # shell of core/shell pair 2
   (...)

   Bonds   # Bond topology for spring forces

   1     2     1     2   # spring for core/shell pair 1
   2     2     3     4   # spring for core/shell pair 2
   (...)

Non-Coulombic (e.g. Lennard-Jones) pairwise interactions are only
defined between the shells.  Coulombic interactions are defined
between all cores and shells.  If desired, additional bonds can be
specified between cores.

The :doc:`special_bonds <special_bonds>` command should be used to
turn-off the Coulombic interaction within core/shell pairs, since that
interaction is set by the bond spring.  This is done using the
:doc:`special_bonds <special_bonds>` command with a 1-2 weight = 0.0,
which is the default value.  It needs to be considered whether one has
to adjust the :doc:`special_bonds <special_bonds>` weighting according
to the molecular topology since the interactions of the shells are
bypassed over an extra bond.

Note that this core/shell implementation does not require all ions to
be polarized.  One can mix core/shell pairs and ions without a
satellite particle if desired.

Since the core/shell model permits distances of r = 0.0 between the
core and shell, a pair style with a "cs" suffix needs to be used to
implement a valid long-range Coulombic correction.  Several such pair
styles are provided in the CORESHELL package.  See :doc:`this doc page <pair_cs>` for details.  All of the core/shell enabled pair
styles require the use of a long-range Coulombic solver, as specified
by the :doc:`kspace_style <kspace_style>` command.  Either the PPPM or
Ewald solvers can be used.

For the NaCL example problem, these pair style and bond style settings
are used:


.. code-block:: LAMMPS

   pair_style      born/coul/long/cs 20.0 20.0
   pair_coeff      * *      0.0 1.000   0.00  0.00   0.00
   pair_coeff      3 3    487.0 0.23768 0.00  1.05   0.50 #Na-Na
   pair_coeff      3 4 145134.0 0.23768 0.00  6.99   8.70 #Na-Cl
   pair_coeff      4 4 405774.0 0.23768 0.00 72.40 145.40 #Cl-Cl

   bond_style      harmonic
   bond_coeff      1 63.014 0.0
   bond_coeff      2 25.724 0.0

When running dynamics with the adiabatic core/shell model, the
following issues should be considered.  The relative motion of
the core and shell particles corresponds to the polarization,
hereby an instantaneous relaxation of the shells is approximated
and a fast core/shell spring frequency ensures a nearly constant
internal kinetic energy during the simulation.
Thermostats can alter this polarization behavior, by scaling the
internal kinetic energy, meaning the shell will not react freely to
its electrostatic environment.
Therefore it is typically desirable to decouple the relative motion of
the core/shell pair, which is an imaginary degree of freedom, from the
real physical system.  To do that, the :doc:`compute temp/cs <compute_temp_cs>` command can be used, in conjunction with
any of the thermostat fixes, such as :doc:`fix nvt <fix_nh>` or :doc:`fix langevin <fix_langevin>`.  This compute uses the center-of-mass velocity
of the core/shell pairs to calculate a temperature, and insures that
velocity is what is rescaled for thermostatting purposes.  This
compute also works for a system with both core/shell pairs and
non-polarized ions (ions without an attached satellite particle).  The
:doc:`compute temp/cs <compute_temp_cs>` command requires input of two
groups, one for the core atoms, another for the shell atoms.
Non-polarized ions which might also be included in the treated system
should not be included into either of these groups, they are taken
into account by the *group-ID* (2nd argument) of the compute.  The
groups can be defined using the :doc:`group *type*\ <group>` command.
Note that to perform thermostatting using this definition of
temperature, the :doc:`fix modify temp <fix_modify>` command should be
used to assign the compute to the thermostat fix.  Likewise the
:doc:`thermo_modify temp <thermo_modify>` command can be used to make
this temperature be output for the overall system.

For the NaCl example, this can be done as follows:


.. code-block:: LAMMPS

   group cores type 1 2
   group shells type 3 4
   compute CSequ all temp/cs cores shells
   fix thermoberendsen all temp/berendsen 1427 1427 0.4    # thermostat for the true physical system
   fix thermostatequ all nve                               # integrator as needed for the berendsen thermostat
   fix_modify thermoberendsen temp CSequ
   thermo_modify temp CSequ                                # output of center-of-mass derived temperature

The pressure for the core/shell system is computed via the regular
LAMMPS convention by :ref:`treating the cores and shells as individual particles <MitchellFincham2>`. For the thermo output of the pressure
as well as for the application of a barostat, it is necessary to
use an additional :doc:`pressure <compute_pressure>` compute based on
the default :doc:`temperature <compute_temp>` and specifying it as a
second argument in :doc:`fix modify <fix_modify>` and
:doc:`thermo_modify <thermo_modify>` resulting in:


.. code-block:: LAMMPS

   (...)
   compute CSequ all temp/cs cores shells
   compute thermo_press_lmp all pressure thermo_temp       # pressure for individual particles
   thermo_modify temp CSequ press thermo_press_lmp         # modify thermo to regular pressure
   fix press_bar all npt temp 300 300 0.04 iso 0 0 0.4
   fix_modify press_bar temp CSequ press thermo_press_lmp  # pressure modification for correct kinetic scalar

If :doc:`compute temp/cs <compute_temp_cs>` is used, the decoupled
relative motion of the core and the shell should in theory be
stable.  However numerical fluctuation can introduce a small
momentum to the system, which is noticeable over long trajectories.
Therefore it is recommendable to use the :doc:`fix momentum <fix_momentum>` command in combination with :doc:`compute temp/cs <compute_temp_cs>` when equilibrating the system to
prevent any drift.

When initializing the velocities of a system with core/shell pairs, it
is also desirable to not introduce energy into the relative motion of
the core/shell particles, but only assign a center-of-mass velocity to
the pairs.  This can be done by using the *bias* keyword of the
:doc:`velocity create <velocity>` command and assigning the :doc:`compute temp/cs <compute_temp_cs>` command to the *temp* keyword of the
:doc:`velocity <velocity>` command, e.g.


.. code-block:: LAMMPS

   velocity all create 1427 134 bias yes temp CSequ
   velocity all scale 1427 temp CSequ

To maintain the correct polarizability of the core/shell pairs, the
kinetic energy of the internal motion shall remain nearly constant.
Therefore the choice of spring force and mass ratio need to ensure
much faster relative motion of the 2 atoms within the core/shell pair
than their center-of-mass velocity. This allows the shells to
effectively react instantaneously to the electrostatic environment and
limits energy transfer to or from the core/shell oscillators.
This fast movement also dictates the timestep that can be used.

The primary literature of the adiabatic core/shell model suggests that
the fast relative motion of the core/shell pairs only allows negligible
energy transfer to the environment.
The mentioned energy transfer will typically lead to a small drift
in total energy over time.  This internal energy can be monitored
using the :doc:`compute chunk/atom <compute_chunk_atom>` and :doc:`compute temp/chunk <compute_temp_chunk>` commands.  The internal kinetic
energies of each core/shell pair can then be summed using the sum()
special function of the :doc:`variable <variable>` command.  Or they can
be time/averaged and output using the :doc:`fix ave/time <fix_ave_time>`
command.  To use these commands, each core/shell pair must be defined
as a "chunk".  If each core/shell pair is defined as its own molecule,
the molecule ID can be used to define the chunks.  If cores are bonded
to each other to form larger molecules, the chunks can be identified
by the :doc:`fix property/atom <fix_property_atom>` via assigning a
core/shell ID to each atom using a special field in the data file read
by the :doc:`read_data <read_data>` command.  This field can then be
accessed by the :doc:`compute property/atom <compute_property_atom>`
command, to use as input to the :doc:`compute chunk/atom <compute_chunk_atom>` command to define the core/shell
pairs as chunks.

For example if core/shell pairs are the only molecules:


.. code-block:: LAMMPS

   read_data NaCl_CS_x0.1_prop.data
   compute prop all property/atom molecule
   compute cs_chunk all chunk/atom c_prop
   compute cstherm all temp/chunk cs_chunk temp internal com yes cdof 3.0     # note the chosen degrees of freedom for the core/shell pairs
   fix ave_chunk all ave/time 10 1 10 c_cstherm file chunk.dump mode vector

For example if core/shell pairs and other molecules are present:


.. code-block:: LAMMPS

   fix csinfo all property/atom i_CSID                       # property/atom command
   read_data NaCl_CS_x0.1_prop.data fix csinfo NULL CS-Info  # atom property added in the data-file
   compute prop all property/atom i_CSID
   (...)

The additional section in the date file would be formatted like this:


.. parsed-literal::

   CS-Info         # header of additional section

   1   1           # column 1 = atom ID, column 2 = core/shell ID
   2   1
   3   2
   4   2
   5   3
   6   3
   7   4
   8   4
   (...)


----------


.. _MitchellFincham:



**(Mitchell and Fincham)** Mitchell, Fincham, J Phys Condensed Matter,
5, 1031-1038 (1993).

.. _MitchellFincham2:



**(Fincham)** Fincham, Mackrodt and Mitchell, J Phys Condensed Matter,
6, 393-404 (1994).
