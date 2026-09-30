.. index:: pair_style runner

pair_style runner command
=========================

Syntax
------

.. code-block:: LAMMPS

   pair_style runner keyword value ...

* zero or more keyword/value pairs may be appended
* keyword = *dir* or *cflength* or *cfenergy* or *f_comm* or *q_comm* or *committee_size* or *total_charge* or *use_prev_q* or *check_extrap* or *max_extrap* or *show_ew* or *sum_ew_freq* or *reset_ew_freq*
* value depends on the preceding keyword:

.. list-table::
   :widths: 20 20 60
   :header-rows: 1

   * - Keyword
     - Value
     - Description
   * - *dir*
     - path
     - Path to RuNNer configuration files
   * - *cflength*
     - factor
     - Length unit conversion factor
   * - *cfenergy*
     - factor
     - Energy unit conversion factor
   * - *committee_size*
     - N
     - Number of committee members
   * - *f_comm*
     - *yes* or *no*
     - Write individual committee forces to per-atom array
   * - *q_comm*
     - *yes* or *no*
     - Write individual committee charges to per-atom array
   * - *total_charge*
     - Q
     - Total system charge
   * - *use_prev_q*
     - *yes* or *no*
     - Use charges from previous step as QEq guess
   * - *check_extrap*
     - *yes* or *no*
     - Enable extrapolation warning (EW) monitoring
   * - *max_extrap*
     - N
     - Stop simulation if EW count exceeds N
   * - *show_ew*
     - *yes* or *no*
     - Write EWs to the log file
   * - *sum_ew_freq*
     - N
     - Write EW summary every N steps
   * - *reset_ew_freq*
     - N
     - Reset EW counters every N steps

Examples
--------

.. code-block:: LAMMPS

   pair_style runner dir "./potential_files"
   pair_coeff * * 1 8

   fix 1 all property/atom d2_f_comm 24 ghost yes
   pair_style runner dir "./potential_files" cflength 1.8897261328 &
      cfenergy 0.0367493254 committee_size 8 f_comm yes
   pair_coeff * * 1 3 8 25

   fix 1 all property/atom d2_q_comm 4 ghost yes
   pair_style runner dir "./potential_files" committee_size 4 q_comm yes total_charge 0.0
   pair_coeff * * 8 12 79

Description
-----------

.. versionadded:: 4Jul2026

This pair style provides an interface to the `RuNNer 2
<https://gitlab.com/runner-suite/runner2>`_ (Ruhr University Neural
Network Energy Representation) library. It implements High-Dimensional
Neural Network Potentials (HDNNPs) as introduced in (:ref:`Behler and
Parrinello 2007 <Behler_2007>`).  HDNNPs are machine learning potentials
that represent the total energy of a system as a sum of
environment-dependent atomic contributions.

The pair style supports several "generations" of HDNNPs as categorized
in (:ref:`Behler 2021 <Behler_2021>`):

* **Second-generation (2G):** Short-range many-body potentials where the
  total energy is the sum of atomic energies predicted from local
  chemical environments (:ref:`Behler and Parrinello 2007
  <Behler_2007>`).
* **Third-generation (3G):** Extends 2G by adding explicit long-range
  electrostatic interactions based on environment-dependent atomic
  partial charges (:ref:`Artrith, Morawietz and Behler 2011
  <Artrith_Morawietz_Behler_2011>`).
* **Fourth-generation (4G):** Includes global charge equilibration (QEq)
  based on environment-dependent electronegativities. These charges are
  used to calculate long-range electrostatics and serve as a global
  descriptor for the atomic energies, allowing for the description of
  nonlocal charge transfer (:ref:`Ko et al 2021
  <Ko_Finkler_Goedecker_Behler_2021>`; :ref:`Gubler et al 2024
  <Gubler_Schaefer_Behler_Goedecker>`).

Additionally, all generations can be augmented with:

* **Hirshfeld-based dispersion:** Long-range van der Waals interactions
  based on the Tkatchenko-Scheffler dispersion model (:ref:`Tkatchenko
  and Scheffler 2009 <Tkatchenko_Scheffler_2009>`).
* **Repulsive potentials:** Screened nuclear repulsion at short
  interatomic distances based on the Ziegler-Biersack-Littmark (ZBL)
  model.

Only a single :doc:`pair_coeff <pair_coeff>` command with two asterisk
wildcards is used with this pair style. Its additional arguments define
the mapping of LAMMPS atom types to RuNNer atomic numbers.

.. code-block:: LAMMPS

   pair_coeff * * 1 8

The example above maps LAMMPS atom types 1 and 2 to atomic numbers 1
("H") and 8 ("O") in RuNNer.

----

General
^^^^^^^

For more information on the options of RuNNer, `visit the official
RuNNer documentation. <https://runner-suite.gitlab.io/runner2/>`_.

Use the *dir* keyword to specify the directory containing the RuNNer
configuration files. The directory must contain:

* ``input.nn``: The HDNNP generation, architecture and feature map
  specifications.
* ``scaling_*.data``: Feature map scaling data for the model.
* ``weights_*.???.data``: Neural network parameters (weights and biases) for each element and model.

The RuNNer library is unit-agnostic.  Use *cflength* and *cfenergy* to
scale LAMMPS coordinates and energies to the units in which the
potential was trained (typically Bohr and Hartree).  If the HDNNP was
trained in Bohr/Hartree and the LAMMPS simulation uses *metal* units
(Angstroms, eV), then *cflength* and *cfenergy* must be the
multiplicative factors required to convert LAMMPS units to the
respective quantities in native HDNNP units:

.. math::

   R_{\text{native}} = R_{\text{LAMMPS}} \times \text{cflength}

.. math::

   E_{\text{native}} = E_{\text{LAMMPS}} \times \text{cfenergy}



Example for *metal* units to a Bohr/Hartree potential:

.. code-block:: LAMMPS

   cflength 1.8897261328   # Angstrom to Bohr
   cfenergy 0.0367493254   # eV to Hartree

Since machine learning potentials are most reliable within their
training data range, the *runner* pair style can monitor whether the
features representing local atomic environments extrapolate beyond the
training range.

* Set *check_extrap* to *yes* to enable monitoring.
* The keyword *show_ew* enables writing individual extrapolation
  warnings (EWs) to the log file.
* Use *sum_ew_freq* to write a summary of EW counts at specific
  intervals instead of logging every occurrence.
* The *max_extrap* threshold allows termination of a simulation if the
  total EW count exceeds this value. Setting *max_extrap* to a negative
  number disables the termination threshold.
* Use *reset_ew_freq* to reset the EW counters at specific intervals.

Committees
^^^^^^^^^^

The pair style supports **Committees**, where multiple HDNNPs sharing
atomic descriptors are evaluated simultaneously.  The forces, energies,
and virials used to propagate the simulation are the average of all
committee members.  This is useful for Query-by-Committee-based Active
Learning approaches and uncertainty estimation for production
simulations.

In the case of a *committee_size* greater than 1, the *dir* keyword must
point to a directory that contains the global configuration files
(``input.nn`` and ``scaling_*.data``) and a set of subdirectories named
**1** to **N** (where N is the *committee_size*).  Each subdirectory must
contain the ``weights_*.???.data`` files for that specific committee
member.

.. code-block:: text

  potential_files/
  |-- 1/
  |   |-- weights_short.001.data
  |   `-- weights_short.008.data
  |-- 2/
  |   |-- weights_short.001.data
  |   `-- weights_short.008.data
  |-- input.nn
  `-- scaling.data

**Accessing Member Energies**

The individual potential energies of each committee member can be
accessed using the :doc:`compute pair <compute_pair>` command:

.. code-block:: LAMMPS

   compute e_comm all pair runner

The energies are stored in a global vector *e_comm* of length
*committee_size*. They can be accessed in subsequent commands (like
:doc:`thermo_style <thermo_style>`) as:

* ``c_e_comm[1]``: Potential energy of member 1
* ``c_e_comm[2]``: Potential energy of member 2
* ``c_e_comm[N]``: Potential energy of member N

**Accessing Member Forces and Charges**

To extract individual member forces or charges (e.g., to compute
per-atom variance), you must define a custom per-atom array using the
:doc:`fix property/atom <fix_property_atom>` command **before** the
*pair_style* command.

For forces, set *f_comm* to *yes*.  The array **must** be a
floating-point array (type ``d2``) named *f_comm* with three times
*committee_size* columns.  Ghost atom communication must be enabled.

.. code-block:: LAMMPS

   # Example: committee_size 8 requires 8 * 3 = 24 columns
   fix 1 all property/atom d2_f_comm 24 ghost yes
   pair_style runner committee_size 8 f_comm yes

Forces are stored sequentially by member and dimension (*fx1*, *fy1*,
*fz1*, *fx2*, ...) and can be accessed by subsequent commands (like
:doc:`dump <dump>`) as:

* ``d2_f_comm[1]``, ``d2_f_comm[2]``, ``d2_f_comm[3]``: Atomic forces (x,y,z) for member 1
* ``d2_f_comm[4]``, ``d2_f_comm[5]``, ``d2_f_comm[6]``: Atomic forces (x,y,z) for member 2

For 3G and 4G potentials, set *q_comm* to *yes* to extract individual
member atomic charges. The array **must** be a floating-point array
(type ``d2``) named *q_comm* with *committee_size* columns.

.. code-block:: LAMMPS

   # Example: committee_size 2
   fix 2 all property/atom d2_q_comm 2 ghost yes
   pair_style runner committee_size 2 q_comm yes

The committee atomic charges can be accessed by subsequent commands
(like :doc:`dump <dump>`) as:

* ``d2_q_comm[1]``: Atomic charges for member 1
* ``d2_q_comm[2]``: Atomic charges for member 2
* ``d2_q_comm[N]``: Atomic charges for member N


3G / 4G only
^^^^^^^^^^^^

For 3G and 4G HDNNPs, the total charge of the system can be specified
with the *total_charge* keyword.  For periodic systems, the system must
be charge neutral.

For 4G HDNNPs only, the *use_prev_q* keyword controls the initial guess
for the iterative Charge Equilibration (QEq) solver.  Setting this to
*yes* uses the predicted charges from the previous time step as the
starting point for the current time step, which can reduce the number of
iterations required for convergence.

.. admonition:: Periodicity limitations
   :class: note

   3G and 4G HDNNPs require either full 3D periodicity (``boundary p p p``) or no periodicity (``boundary f f f``). Partial periodicity (e.g., ``boundary p p f``) is currently not supported for long-range electrostatics.

.. admonition:: MPI scaling
   :class: note

   3G and 4G HDNNPs require a global collection of the atomic structure
   on a single process to perform long-range electrostatics and global
   QEq calculations.  This creates an MPI bottleneck that leads to
   inferior strong scaling as the number of MPI tasks increases.  Since
   3G and 4G HDNNPs are heavily optimized for OpenMP in RuNNer, it is
   highly recommended to use a small number of MPI tasks and a large
   number of OpenMP threads per task to achieve the best performance.

----

Mixing, shift, table, tail correction, restart, rRESPA info
-----------------------------------------------------------

This style does not support mixing.  The :doc:`pair_coeff <pair_coeff>`
command should only be invoked with asterisk wildcards to define the
mapping for all atom types.

This style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This style does not write information to :doc:`binary restart files <restart>`.
You must re-specify the *pair_style* and *pair_coeff* commands in any input
script that reads a restart file.

This style can only be used via the *pair* keyword of the :doc:`run_style respa <run_style>`
command.  It does not support the *inner*, *middle*, or *outer* keywords.

Restrictions
------------

This pair style is part of the ML-RUNNER package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` doc page for more info.

Currently, only one instance of ``pair_style runner`` can be initialized
per simulation.  The style does not support the use of :doc:`pair_style
hybrid <pair_hybrid>` where multiple ``runner`` instances are defined.

Related commands
----------------

:doc:`pair_coeff <pair_coeff>`, :doc:`fix property/atom <fix_property_atom>`,
:doc:`compute pair <compute_pair>`

Default
-------

The default options are:

* *dir* = "./"
* *cflength* = 1.0
* *cfenergy* = 1.0
* *committee_size* = 1
* *f_comm* = no
* *q_comm* = no
* *total_charge* = 0.0
* *use_prev_q* = no
* *check_extrap* = no
* *max_extrap* = 100
* *show_ew* = no
* *sum_ew_freq* = 0
* *reset_ew_freq* = 0

----

References
----------

.. _Behler_2007:

**(Behler and Parrinello 2007)** Behler, J.; Parrinello, M., Phys. Rev. Lett. 2007, 98 (14), 146401.

.. _Tkatchenko_Scheffler_2009:

**(Tkatchenko and Scheffler 2009)** Tkatchenko, A.; Scheffler, M., Phys. Rev. Lett. 2009, 102 (7), 073005.

.. _Behler_2011:

**(Behler 2011)** Behler, J., J. Chem. Phys. 2011, 134 (7), 074106.

.. _Artrith_Morawietz_Behler_2011:

**(Artrith, Morawietz and Behler 2011)** Artrith, N.; Morawietz, T.; Behler, J., Phys. Rev. B 2011, 83 (15), 153101.

.. _Ko_Finkler_Goedecker_Behler_2021:

**(Ko et al 2021)** Ko, T. W.; Finkler, J. A.; Goedecker, S.; Behler, J., Nat. Commun. 2021, 12, 398.

.. _Behler_2021:

**(Behler 2021)** Behler, J., Chem. Rev. 2021, 121 (16), 10037-10072.

.. _Gubler_Schaefer_Behler_Goedecker:

**(Gubler et al 2024)** Gubler, M.; Finkler, J. A.; Schaefer, M. R.; Behler, J.; Goedecker S., J. Chem. Theory Comput. 2024, 20 (16), 7264-7271.
