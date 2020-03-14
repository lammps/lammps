.. index:: kim_init, kim_interactions, kim_query, kim_param

kim_init command
=================

kim_interactions command
=========================

kim_query command
==================

kim_param command
==================

Syntax
""""""

.. code-block:: LAMMPS

   kim_init model user_units unitarg
   kim_interactions typeargs
   kim_query variable formatarg query_function queryargs
   kim_param get param_name index_range variables formatarg
   kim_param set param_name index_range values

.. _formatarg_options:

* model = name of the KIM interatomic model (the KIM ID for models archived in OpenKIM)
* user_units = the LAMMPS :doc:`units <units>` style assumed in the LAMMPS input script
* unitarg = *unit_conversion_mode* (optional)
* typeargs = atom type to species mapping (one entry per atom type) or *fixed_types* for models with a preset fixed mapping
* variable(s) = single name or list of names of (string style) LAMMPS variable(s) where a query result or parameter get result is stored. Variables that do not exist will be created by the command.
* formatarg = *list, split, or explicit* (optional):

  .. parsed-literal::

     *list* = returns a single string with a list of space separated values
            (e.g. "1.0 2.0 3.0"), which is placed in a LAMMPS variable as
            defined by the *variable* argument. [default for *kim_query*]
     *split* = returns the values separately in new variables with names based
            on the prefix specified in *variable* and a number appended to
            indicate which element in the list of values is in the variable.
     *explicit* = returns the values separately in one more more variable names
            provided as arguments that precede *formatarg*\ . [default for *kim_param*]

* query_function = name of the OpenKIM web API query function to be used
* queryargs = a series of *keyword=value* pairs that represent the web query; supported keywords depend on the query function
* param_name = name of a KIM portable model parameter
* index_range = KIM portable model parameter index range (an integer for a single element, or pair of integers separated by a colon for a range of elements)
* values = new value(s) to replace the current value(s) of a KIM portable model parameter

Examples
""""""""

.. code-block:: LAMMPS

   kim_init SW_StillingerWeber_1985_Si__MO_405512056662_005 metal
   kim_interactions Si
   kim_init Sim_LAMMPS_ReaxFF_StrachanVanDuinChakraborty_2003_CHNO__SM_107643900657_000 real
   kim_init Sim_LAMMPS_ReaxFF_StrachanVanDuinChakraborty_2003_CHNO__SM_107643900657_000 metal unit_conversion_mode
   kim_interactions C H O
   kim_init Sim_LAMMPS_IFF_PCFF_HeinzMishraLinEmami_2015Ver1v5_FccmetalsMineralsSolventsPolymers__SM_039297821658_000 real
   kim_interactions fixed_types
   kim_query a0 get_lattice_constant_cubic crystal=["fcc"] species=["Al"] units=["angstrom"]
   kim_param get gamma 1 varGamma
   kim_param set gamma 1 3.0

Description
"""""""""""

The set of *kim_commands* provide a high-level wrapper around the
`Open Knowledgebase of Interatomic Models (OpenKIM) <https://openkim.org>`_
repository of interatomic models (IMs) (potentials and force fields),
so that they can be used by LAMMPS scripts.  These commands do not implement
any computations directly, but rather generate LAMMPS input commands based
on the information retrieved from the OpenKIM repository to initialize and
activate OpenKIM IMs and query their predictions for use in the LAMMPS script.
All LAMMPS input commands generated and executed by *kim_commands* are
echoed to the LAMMPS log file.

Benefits of Using OpenKIM IMs
-----------------------------

Employing OpenKIM IMs provides LAMMPS users with multiple benefits:

Reliability
^^^^^^^^^^^

* All content archived in OpenKIM is reviewed by the `KIM Editor <https://openkim.org/governance/>`_ for quality.
* IMs in OpenKIM are archived with full provenance control. Each is associated with a maintainer responsible for the integrity of the content. All changes are tracked and recorded.
* IMs in OpenKIM are exhaustively tested using `KIM Tests <https://openkim.org/doc/evaluation/kim-tests/>`_ that compute a host of material properties, and `KIM Verification Checks <https://openkim.org/doc/evaluation/kim-verification-checks/>`_ that provide the user with information on various aspects of the IM behavior and coding correctness. This information is displayed on the IM's page accessible through the  `OpenKIM browse interface <https://openkim.org/browse>`_.

Reproducibility
^^^^^^^^^^^^^^^

* Each IM in OpenKIM is issued a unique identifier (`KIM ID <https://openkim.org/doc/schema/kim-ids/>`_), which includes a version number (last three digits).  Any changes that can result in different numerical values lead to a version increment in the KIM ID. This makes it possible to reproduce simulations since the specific version of a specific IM used can be retrieved using its KIM ID.
* OpenKIM is a member organization of `DataCite <https://datacite.org/>`_ and issues digital object identifiers (DOIs) to all IMs archived in OpenKIM. This makes it possible to cite the IM code used in a simulation in a publications to give credit to the developers and further facilitate reproducibility.

Convenience
^^^^^^^^^^^

* IMs in OpenKIM are distributed in binary form along with LAMMPS and can be used in a LAMMPS input script simply by providing their KIM ID in the *kim_init* command documented on this page.
* The *kim_query* web query tool provides the ability to use the predictions of IMs for supported material properties (computed via `KIM Tests <https://openkim.org/doc/evaluation/kim-tests/>`_) as part of a LAMMPS input script setup and analysis.
* Support is provided for unit conversion between the :doc:`unit style <units>` used in the LAMMPS input script and the units required by the OpenKIM IM. This makes it possible to use a single input script with IMs using different units without change and minimizes the likelihood of errors due to incompatible units.

.. _IM_types:

Types of IMs in OpenKIM
-----------------------

There are two types of IMs archived in OpenKIM:

.. _PM_type:

1. The first type is called a *KIM Portable Model* (PM). A KIM PM is an independent computer implementation of an IM written in one of the languages supported by KIM (C, C++, Fortran) that conforms to the KIM Application Programming Interface (`KIM API <https://openkim.org/kim-api/>`_) Portable Model Interface (PMI) standard. A KIM PM will work seamlessly with any simulation code that supports the KIM API/PMI standard (including LAMMPS; see `complete list of supported codes <https://openkim.org/projects-using-kim/>`_).
2. The second type is called a *KIM Simulator Model* (SM). A KIM SM is an IM that is implemented natively within a simulation code (\ *simulator*\ ) that supports the KIM API Simulator Model Interface (SMI); in this case LAMMPS. A separate SM package is archived in OpenKIM for each parameterization of the IM, which includes all of the necessary parameter files, LAMMPS commands, and metadata (supported species, units, etc.) needed to run the IM in LAMMPS.

With these two IM types, OpenKIM can archive and test almost all IMs that
can be used by LAMMPS. (It is easy to contribute new IMs to OpenKIM, see
the `upload instructions <https://openkim.org/doc/repository/adding-content/>`_.)

OpenKIM IMs are uniquely identified by a
`KIM ID <https://openkim.org/doc/schema/kim-ids/>`_.
The extended KIM ID consists of
a human-readable prefix identifying the type of IM, authors, publication year,
and supported species, separated by two underscores from the KIM ID itself,
which begins with an IM code
(\ *MO* for a KIM Portable Model, and *SM* for a KIM Simulator Model)
followed by a unique 12-digit code and a 3-digit version identifier.
By convention SM prefixes begin with *Sim_* to readily identify them.

.. parsed-literal::

   SW_StillingerWeber_1985_Si__MO_405512056662_005
   Sim_LAMMPS_ReaxFF_StrachanVanDuinChakraborty_2003_CHNO__SM_107643900657_000

Each OpenKIM IM has a dedicated "Model Page" on `OpenKIM <https://openkim.org>`_
providing all the information on the IM including a title, description,
authorship and citation information, test and verification check results,
visualizations of results, a wiki with documentation and user comments, and
access to raw files, and other information.
The URL for the Model Page is constructed from the
`extended KIM ID <https://openkim.org/doc/schema/kim-ids/>`_ of the IM:

.. parsed-literal::

   https://openkim.org/id/extended_KIM_ID

For example, for the Stillinger--Weber potential
listed above the Model Page is located at:

.. parsed-literal::

   `https://openkim.org/id/SW_StillingerWeber_1985_Si__MO_405512056662_005 <https://openkim.org/id/SW_StillingerWeber_1985_Si__MO_405512056662_005>`_

See the `current list of KIM PMs and SMs archived in OpenKIM <https://openkim.org/browse/models/by-species>`_.
This list is sorted by species and can be filtered to display only
IMs for certain species combinations.

See `Obtaining KIM Models <http://openkim.org/doc/usage/obtaining-models>`_ to
learn how to install a pre-build binary of the OpenKIM Repository of Models.

.. note::

   It is also possible to locally install IMs not archived in OpenKIM,
   in which case their names do not have to conform to the KIM ID format.

Using OpenKIM IMs with LAMMPS
-----------------------------

Two commands are employed when using OpenKIM IMs, one to select the
IM and perform necessary initialization (*kim_init*), and the second
to set up the IM for use by executing any necessary LAMMPS commands
(*kim_interactions*). Both are required.

See the *examples/kim* directory for example input scripts that use KIM PMs
and KIM SMs.

OpenKIM IM Initialization (*kim_init*)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *kim_init* mode command must be issued **before**
the simulation box is created (normally at the top of the file).
This command sets the OpenKIM IM that will be used and may issue
additional commands changing LAMMPS default settings that are required
for using the selected IM (such as :doc:`units <units>` or
:doc:`atom_style <atom_style>`). If needed, those settings can be overridden,
however, typically a script containing a *kim_init* command
would not include *units* and *atom_style* commands.

The required arguments of *kim_init* are the *model* name of the
IM to be used in the simulation (for an IM archived in OpenKIM this is
its `extended KIM ID <https://openkim.org/doc/schema/kim-ids/>`_, and
the *user_units*, which are the LAMMPS :doc:`units style <units>` used
in the input script.  (Any dimensioned numerical values in the input
script and values read in from files are expected to be in the
*user_units* system.)

The selected IM can be either a :ref:`KIM PM or a KIM SM <IM_types>`.
For a KIM SM, the *kim_init* command verifies that the SM is designed
to work with LAMMPS (and not another simulation code).
In addition, the LAMMPS version used for defining
the SM and the LAMMPS version being currently run are
printed to help diagnose any incompatible changes to input script or
command syntax between the two LAMMPS versions.

Based on the selected model *kim_init* may modify the
:doc:`atom_style <atom_style>`.
Some SMs have requirements for this setting. If this is the case, then
*atom_style* will be set to the required style. Otherwise, the value is left
unchanged (which in the absence of an *atom_style* command in the input script
is the :doc:`default atom_style value <atom_style>`).

Regarding units, the *kim_init* command behaves in different ways depending
on whether or not *unit conversion mode* is activated as indicated by the
optional *unitarg* argument.
If unit conversion mode is **not** active, then *user_units* must
either match the required units of the IM or the IM must be able
to adjust its units to match. (The latter is only possible with some KIM PMs;
SMs can never adjust their units.) If a match is possible, the LAMMPS
:doc:`units <units>` command is called to set the units to
*user_units*. If the match fails, the simulation is terminated with
an error.

Here is an example of a LAMMPS script to compute the cohesive energy
of a face-centered cubic (fcc) lattice for the Ercolessi and Adams (1994)
potential for Al:

.. code-block:: LAMMPS

   kim_init         EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005 metal
   boundary         p p p
   lattice          fcc 4.032
   region           simbox block 0 1 0 1 0 1 units lattice
   create_box       1 simbox
   create_atoms     1 box
   mass             1 26.981539
   kim_interactions Al
   run              0
   variable         Ec equal (pe/count(all))/${_u_energy}
   print            "Cohesive Energy = ${EcJ} eV"

The above script will end with an error in the *kim_init* line if the
IM is changed to another potential for Al that does not work with *metal*
units. To address this *kim_init* offers the *unit_conversion_mode*
as shown below.
If unit conversion mode *is* active, then *kim_init* calls the LAMMPS
:doc:`units <units>` command to set the units to the IM's required or
preferred units. Conversion factors between the IM's units and the *user_units*
are defined for all :doc:`physical quantities <units>` (mass, distance, etc.).
(Note that converting to or from the "lj" unit style is not supported.)
These factors are stored as :doc:`internal style variables <variable>` with
the following standard names:

.. parsed-literal::

   _u_mass
   _u_distance
   _u_time
   _u_energy
   _u_velocity
   _u_force
   _u_torque
   _u_temperature
   _u_pressure
   _u_viscosity
   _u_charge
   _u_dipole
   _u_efield
   _u_density

If desired, the input script can be designed to work with these conversion
factors so that the script will work without change with any OpenKIM IM.
(This approach is used in the
`OpenKIM Testing Framework <https://openkim.org/doc/evaluation/kim-tests/>`_.)
For example, the script given above for the cohesive energy of fcc Al
can be rewritten to work with any IM regardless of units. The following
script constructs an fcc lattice with a lattice parameter defined in
meters, computes the total energy, and prints the cohesive energy in
Joules regardless of the units of the IM.

.. code-block:: LAMMPS

   kim_init         EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005 si unit_conversion_mode
   boundary         p p p
   lattice          fcc 4.032e-10*${_u_distance}
   region           simbox block 0 1 0 1 0 1 units lattice
   create_box       1 simbox
   create_atoms     1 box
   mass             1 4.480134e-26*${_u_mass}
   kim_interactions Al
   run              0
   variable         Ec_in_J equal (pe/count(all))/${_u_energy}
   print            "Cohesive Energy = ${Ec_in_J} J"

Note the multiplication by ${_u_distance} and ${_u_mass} to convert
from SI units (specified in the *kim_init* command) to whatever units the
IM uses (metal in this case), and the division by ${_u_energy}
to convert from the IM's energy units to SI units (Joule). This script
will work correctly for any IM for Al (KIM PM or SM) selected by the
*kim_init* command.

Care must be taken to apply unit conversion to dimensional variables read in
from a file. For example, if a configuration of atoms is read in from a
dump file using the :doc:`read_dump <read_dump>` command, the following can
be done to convert the box and all atomic positions to the correct units:

.. code-block:: LAMMPS

   variable xyfinal equal xy*${_u_distance}
   variable xzfinal equal xz*${_u_distance}
   variable yzfinal equal yz*${_u_distance}
   change_box all x scale ${_u_distance} &
                          y scale ${_u_distance} &
                          z scale ${_u_distance} &
                          xy final ${xyfinal} &
                          xz final ${xzfinal} &
                          yz final ${yzfinal} &
                          remap

.. note::

   Unit conversion will only work if the conversion factors are placed in
   all appropriate places in the input script. It is up to the user to do this
   correctly.

OpenKIM IM Execution (*kim_interactions*)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The second and final step in using an OpenKIM IM is to execute the
*kim_interactions* command. This command must be preceded by a *kim_init*
command and a command that defines the number of atom types *N* (such as
:doc:`create_box <create_box>`).
The *kim_interactions* command has one argument *typeargs*\ . This argument
contains either a list of *N* chemical species, which defines a mapping between
atom types in LAMMPS to the available species in the OpenKIM IM, or the
keyword *fixed_types* for models that have a preset fixed mapping (i.e.
the mapping between LAMMPS atom types and chemical species is defined by
the model and cannot be changed). In the latter case, the user must consult
the model documentation to see how many atom types there are and how they
map to the chemical species.

For example, consider an OpenKIM IM that supports Si and C species.
If the LAMMPS simulation has four atom types, where the first three are Si,
and the fourth is C, the following *kim_interactions* command would be used:

.. code-block:: LAMMPS

   kim_interactions Si Si Si C

Alternatively, for a model with a fixed mapping the command would be:

.. code-block:: LAMMPS

   kim_interactions fixed_types

The *kim_interactions* command performs all the necessary steps to set up
the OpenKIM IM selected in the *kim_init* command. The specific actions depend
on whether the IM is a KIM PM or a KIM SM.  For a KIM PM,
a :doc:`pair_style kim <pair_kim>` command is executed followed by
the appropriate *pair_coeff* command. For example, for the
Ercolessi and Adams (1994) KIM PM for Al set by the following commands:

.. code-block:: LAMMPS

   kim_init EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005 metal
   ...
   ...  box specification lines skipped
   ...
   kim_interactions Al

the *kim_interactions* command executes the following LAMMPS input commands:

.. code-block:: LAMMPS

   pair_style kim EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005
   pair_coeff * * Al

For a KIM SM, the generated input commands may be more complex
and require that LAMMPS is built with the required packages included
for the type of potential being used. The set of commands to be executed
is defined in the SM specification file, which is part of the SM package.
For example, for the Strachan et al. (2003) ReaxFF SM
set by the following commands:

.. code-block:: LAMMPS

   kim_init Sim_LAMMPS_ReaxFF_StrachanVanDuinChakraborty_2003_CHNO__SM_107643900657_000 real
   ...
   ...  box specification lines skipped
   ...
   kim_interactions C H N O

the *kim_interactions* command executes the following LAMMPS input commands:

.. code-block:: LAMMPS

   pair_style reax/c lmp_control safezone 2.0 mincap 100
   pair_coeff * * ffield.reax.rdx C H N O
   fix reaxqeq all qeq/reax 1 0.0 10.0 1.0e-6 param.qeq

Note that the files *lmp_control*, *ffield.reax.rdx* and *param.qeq*
are specific to the Strachan et al. (2003) ReaxFF parameterization
and are archived as part of the SM package in OpenKIM.
Note also that parameters like cutoff radii and charge tolerances,
which have an effect on IM predictions, are also included in the
SM definition ensuring reproducibility.

.. note::

   When using *kim_init* and *kim_interactions* to select
   and set up an OpenKIM IM, other LAMMPS commands
   for the same functions (such as pair_style, pair_coeff, bond_style,
   bond_coeff, fixes related to charge equilibration, etc.) should normally
   not appear in the input script.

Using OpenKIM Web Queries in LAMMPS (*kim_query*)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *kim_query* command performs a web query to retrieve the predictions
of an IM set by *kim_init* for material properties archived in
`OpenKIM <https://openkim.org>`_.

.. note::

   The *kim_query* command must be preceded by a *kim_init* command.

The syntax for the *kim_query* command is as follows:

.. code-block:: LAMMPS

   kim_query variable formatarg query_function queryargs

The result of the query is stored in one or more
:doc:`string style variables <variable>` as determined by the
optional *formatarg* argument :ref:`documented above <formatarg_options>`.
For the "list" setting of *formatarg* (or if *formatarg* is not
specified), the result is returned as a space-separated list of
values in *variable*\ .
The *formatarg* keyword "split" separates the result values into
individual variables of the form *prefix_I*, where *prefix* is set to the
*kim_query* *variable* argument and *I* ranges from 1 to the number of
returned values. The number and order of the returned values is determined
by the type of query performed.  (Note that the "explicit" setting of
*formatarg* is not supported by *kim_query*.)

.. note::

   *kim_query* only supports queries that return a single result or
   an array of values. More complex queries that return a JSON structure
   are not currently supported. An attempt to use *kim_query* in such
   cases will generate an error.

The second required argument *query_function* is the name of the
query function to be called (e.g. *get_lattice_constant_cubic*).
All following :doc:`arguments <Commands_parse>` are parameters handed over to
the web query in the format *keyword=value*\ , where *value* is always
an array of one or more comma-separated items in brackets.
The list of supported keywords and the type and format of their values
depend on the query function used. The current list of query functions
is available on the OpenKIM webpage at
`https://openkim.org/doc/usage/kim-query <https://openkim.org/doc/usage/kim-query>`_.

.. note::

   All query functions require the *model* keyword, which identifies
   the IM whose predictions are being queried. This keyword is automatically
   generated by *kim_query* based on the IM set in *kim_init* and must not
   be specified as an argument to *kim_query*.

.. note::

   Each *query_function* is associated with a default method (implemented
   as a `KIM Test <https://openkim.org/doc/evaluation/kim-tests/>`_)
   used to compute this property. In cases where there are multiple
   methods in OpenKIM for computing a property, a *method* keyword can
   be provided to select the method of choice.  See the
   `query documentation <https://openkim.org/doc/repository/kim-query>`_
   to see which methods are available for a given *query function*\ .

*kim_query* Usage Examples and Further Clarifications
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The data obtained by *kim_query* commands can be used as part of the setup
or analysis phases of LAMMPS simulations. Some examples are given below.

**Define an equilibrium fcc crystal**

.. code-block:: LAMMPS

   kim_init         EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005 metal
   boundary         p p p
   kim_query        a0 get_lattice_constant_cubic crystal=["fcc"] species=["Al"] units=["angstrom"]
   lattice          fcc ${a0}
   ...

The *kim_query* command retrieves from `OpenKIM <https://openkim.org>`_
the equilibrium lattice constant predicted by the Ercolessi and Adams (1994)
potential for the fcc structure and places it in
variable *a0*\ . This variable is then used on the next line to set up the
crystal. By using *kim_query*, the user is saved the trouble and possible
error of tracking this value down, or of having to perform an energy
minimization to find the equilibrium lattice constant.

Note that in *unit_conversion_mode* the results obtained from a
*kim_query* would need to be converted to the appropriate units system.
For example, in the above script, the lattice command would need to be
changed to: "lattice fcc ${a0}\*${_u_distance}".

**Define an equilibrium hcp crystal**

.. code-block:: LAMMPS

   kim_init         EAM_Dynamo_Mendelev_2007_Zr__MO_848899341753_000 metal
   boundary         p p p
   kim_query        latconst split get_lattice_constant_hexagonal crystal=["hcp"] species=["Zr"] units=["angstrom"]
   variable         a0 equal latconst_1
   variable         c0 equal latconst_2
   variable         c_to_a equal ${c0}/${a0}
   lattice          custom ${a0} a1 0.5 -0.866025 0 a2 0.5 0.866025 0 a3 0 0 ${c_to_a} &
                    basis 0.333333 0.666666 0.25 basis 0.666666 0.333333 0.75
   ...

In this case the *kim_query* returns two arguments (since the hexagonal
close packed (hcp) structure has two independent lattice constants).
The *formatarg* keyword "split" places the two values into
the variables *latconst_1* and *latconst_2*. (These variables are
created if they do not already exist.) For convenience the variables
*a0* and *c0* are created in order to make the remainder of the
input script more readable.

**Define a crystal at finite temperature accounting for thermal expansion**

.. code-block:: LAMMPS

   kim_init         EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005 metal
   boundary         p p p
   kim_query        a0 get_lattice_constant_cubic crystal=["fcc"] species=["Al"] units=["angstrom"]
   kim_query        alpha get_linear_thermal_expansion_coefficient_cubic  crystal=["fcc"] species=["Al"] units=["1/K"] temperature=[293.15] temperature_units=["K"]
   variable         DeltaT equal 300
   lattice          fcc ${a0}*${alpha}*${DeltaT}
   ...

As in the previous example, the equilibrium lattice constant is obtained
for the Ercolessi and Adams (1994) potential. However, in this case the
crystal is scaled to the appropriate lattice constant at room temperature
(293.15 K) by using the linear thermal expansion constant predicted by the
potential.

.. note::

   When passing numerical values as arguments (as in the case
   of the temperature in the above example) it is also possible to pass a
   tolerance indicating how close to the value is considered a match.
   If no tolerance is passed a default value is used. If multiple results
   are returned (indicating that the tolerance is too large), *kim_query*
   will return an error. See the
   `query documentation <https://openkim.org/doc/repository/kim-query>`_
   to see which numerical arguments and tolerances are available for a
   given *query function*\ .

**Compute defect formation energy**

.. code-block:: LAMMPS

   kim_init         EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005 metal
   ...
   ... Build fcc crystal containing some defect and compute the total energy
   ... which is stored in the variable *Etot*
   ...
   kim_query        Ec get_cohesive_energy_cubic crystal=["fcc"] species=["Al"] units=["eV"]
   variable         Eform equal ${Etot} - count(all)*${Ec}
   ...

The defect formation energy *Eform* is computed by subtracting from *Etot* the
ideal fcc cohesive energy of the atoms in the system obtained from
`OpenKIM <https://openkim.org>`_ for the Ercolessi and Adams (1994) potential.

.. note::

   *kim_query* commands return results archived in
   `OpenKIM <https://openkim.org>`_. These results are obtained
   using programs for computing material properties
   (KIM Tests and KIM Test Drivers) that were contributed to OpenKIM.
   In order to give credit to Test developers, the number of times results
   from these programs are queried is tracked. No other information about
   the nature of the query or its source is recorded.

Accessing KIM Model Parameters from LAMMPS (*kim_param*)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All IMs are functional forms containing a set of
parameters.  The values of these parameters are typically
selected to best reproduce a training set of quantum mechanical
calculations or available experimental data.  For example, a
Lennard-Jones potential intended to model argon might have the values of
its two parameters, epsilon and sigma, fit to the
dimer dissociation energy or thermodynamic properties at a critical point
of the phase diagram.

Normally a user employing an IM should not modify its parameters since,
as noted above, these are selected to reproduce material properties.
However, there are cases where accessing and modifying IM parameters
is desired, such as for assessing uncertainty, fitting an IM,
or working with an ensemble of IMs. As explained :ref:`above <IM_types>`,
IMs archived in OpenKIM are either Portable Models (PMs) or
Simulator Models (SMs). KIM PMs are complete independent implementations
of an IM, whereas KIM SMs are wrappers to an IM implemented within LAMMPS.
Two different mechanisms are provided for accessing IM parameters in these
two cases:

* For a KIM PM, the *kim_param* command can be used to *get* and *set* the values of the PM's parameters as explained below.
* For a KIM SM, the user should consult the documentation page for the specific IM and follow instructions there for how to modify its parameters (if possible).

The *kim_param get* and *kim_param set* commands provide an interface
to access and change the parameters of a KIM PM that "publishes" its
parameters and makes them publicly available (see the
`KIM API documentation <https://kim-api.readthedocs.io/en/devel/features.html>`_
for details).

.. note::

   The *kim_param get/set* commands must be preceded by *kim_init*.
   The *kim_param set* command must additionally be preceded by a
   *kim_interactions* command (or alternatively by a *pair_style kim*
   and *pair_coeff* commands).  The *kim_param set* command may be used wherever a *pair_coeff* command may occur.

The syntax for the *kim_param* command is as follows:

.. code-block:: LAMMPS

   kim_param get param_name index_range variable formatarg
   kim_param set param_name index_range values

Here, *param_name* is the name of a KIM PM parameter (which is published
by the PM and available for access). The specific string used to identify
a parameter is defined by the PM. For example, for the
`Stillinger--Weber (SW) potential in OpenKIM <https://openkim.org/id/SW_StillingerWeber_1985_Si__MO_405512056662_005>`_,
the parameter names are *A, B, p, q, sigma, gamma, cutoff, lambda, costheta0*\ .

.. note::

   The list of all the parameters that a PM exposes for access/mutation are
   automatically written to the lammps log file when *kim_init* is called.

Each published parameter of a KIM PM takes the form of an array of
numerical values. The array can contain one element for a single-valued
parameter, or a set of values. For example, the
`multispecies SW potential for the Zn-Cd-Hg-S-Se-Te system <https://openkim.org/id/SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_002>`_
has the same parameter names as the
`single-species SW potential <https://openkim.org/id/SW_StillingerWeber_1985_Si__MO_405512056662_005>`_,
but each parameter array contains 21 entries that correspond to the parameter
values used for each pairwise combination of the model's six supported species
(this model does not have parameters specific to individual ternary
combinations of its supported species).

The *index_range* argument may either be an integer referring to
a specific element within the array associated with the parameter
specified by *param_name*, or a pair of integers separated by a colon
that refer to a slice of this array.  In both cases, one-based indexing is
used to refer to the entries of the array.

The result of a *get* operation for a specific *index_range* is stored in
one or more :doc:`LAMMPS string style variables <variable>` as determined
by the optional *formatarg* argument :ref:`documented above. <formatarg_options>`
If not specified, the default for *formatarg* is "explicit" for the
*kim_param* command.

For the case where the result is an array with multiple values
(i.e. *index_range* contains a range), the optional "split" or "explicit"
*formatarg* keywords can be used to separate the results into multiple
variables; see the examples below.
Multiple parameters can be retrieved with a single call to *kim_param get*
by repeating the argument list following *get*\ .

For a *set* operation, the *values* argument contains the new value(s)
for the element(s) of the parameter specified by *index_range*. For the case
where multiple values are being set, *values* contains a set of values
separated by spaces. Multiple parameters can be set with a single call to
*kim_param set* by repeating the argument list following *set*\ .

*kim_param* Usage Examples and Further Clarifications
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Examples of getting and setting KIM PM parameters with further
clarifications are provided below.

**Getting a scalar parameter**

.. code-block:: LAMMPS

   kim_init         SW_StillingerWeber_1985_Si__MO_405512056662_005 metal
   ...
   kim_param        get A 1 VARA

In this case, the value of the SW *A* parameter is retrieved and placed
in the LAMMPS variable *VARA*\ . The variable *VARA* can be used
in the remainder of the input script in the same manner as any other
LAMMPS variable.

**Getting multiple scalar parameters with a single call**

.. code-block:: LAMMPS

   kim_init         SW_StillingerWeber_1985_Si__MO_405512056662_005 metal
   ...
   kim_param        get A 1 VARA B 1 VARB

This retrieves the *A* and *B* parameters of the SW potential and stores
them in the LAMMPS variables *VARA* and *VARB*\ .

**Getting a range of values from a parameter**

There are several options when getting a range of values from a parameter
determined by the *formatarg* argument.

.. code-block:: LAMMPS

   kim_init         SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_002 metal
   ...
   kim_param        get lambda 7:9 LAM_TeTe LAM_TeZn LAM_TeSe

In this case, *formatarg* is not specified and therefore the default
"explicit" mode is used. (The behavior would be the same if the word
*explicit* were added after *LAM_TeSe*.) Elements 7, 8 and 9 of parameter
lambda retrieved by the *get* operation are placed in the LAMMPS variables
*LAM_TeTe*, *LAM_TeZn* and *LAM_TeSe*, respectively.

.. note::

   In the above example, elements 7--9 of the lambda parameter correspond
   to Te-Te, Te-Zm and Te-Se interactions. This can be determined by visiting
   the `model page for the specified potential <https://openkim.org/id/SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_002>`_
   and looking at its parameter file linked to at the bottom of the page
   (file with .param ending) and consulting the README documentation
   provided with the driver for the PM being used. A link to the driver
   is provided at the top of the model page.

.. code-block:: LAMMPS

   kim_init         SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_002 metal
   ...
   kim_param        get lambda 15:17 LAMS list
   variable         LAM_VALUE index ${LAMS}
   label            loop_on_lambda
   ...
   ... do something with current value of lambda
   ...
   next             LAM_VALUE
   jump             SELF loop_on_lambda

In this case, the "list" mode of *formatarg* is used.
The result of the *get* operation is stored in the LAMMPS variable
*LAMS* as a string containing the three retrieved values separated
by spaces, e.g "1.0 2.0 3.0". This can be used in LAMMPS with an
*index* variable to access the values one at a time within a loop
as shown in the example. At each iteration of the loop *LAM_VALUE*
contains the current value of lambda.

.. code-block:: LAMMPS

   kim_init         SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_002 metal
   ...
   kim_param        get lambda 15:17 LAM split

In this case, the "split" mode of *formatarg* is used.
The three values retrieved by the *get* operation are stored in
the three LAMMPS variables *LAM_15*, *LAM_16* and *LAM_17*.
The provided name "LAM" is used as prefix and the location in
the lambda array is appended to create the variable names.

**Setting a scalar parameter**

.. code-block:: LAMMPS

   kim_init         SW_StillingerWeber_1985_Si__MO_405512056662_005 metal
   ...
   kim_interactions Si
   kim_param        set gamma 1 2.6

Here, the SW potential's gamma parameter is set to 2.6.  Note that the *get*
and *set* commands work together, so that a *get* following a *set*
operation will return the new value that was set. For example:

.. code-block:: LAMMPS

   ...
   kim_interactions Si
   kim_param        get gamma 1 ORIG_GAMMA
   kim_param        set gamma 1 2.6
   kim_param        get gamma 1 NEW_GAMMA
   ...
   print            "original gamma = ${ORIG_GAMMA}, new gamma = ${NEW_GAMMA}"

Here, *ORIG_GAMMA* will contain the original gamma value for the SW
potential, while *NEW_GAMMA* will contain the value 2.6.

**Setting multiple scalar parameters with a single call**

.. parsed-literal::

   kim_init         SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_002 metal
   ...
   kim_interactions Cd Te
   variable        VARG equal 2.6
   variable        VARS equal 2.0951
   kim_param       set gamma 1 ${VARG} sigma 3 ${VARS}

In this case, the first element of the *gamma* parameter and
third element of the *sigma* parameter are set to 2.6 and 2.0951,
respectively. This example also shows how LAMMPS variables can
be used when setting parameters.

**Setting a range of values of a parameter**

.. code-block:: LAMMPS

   kim_init         SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_002 metal
   ...
   kim_interactions Cd Te Zn Se Hg S
   kim_param        set sigma 2:6 2.35214 2.23869 2.04516 2.43269 1.80415

In this case, elements 2 through 6 of the parameter *sigma*
are set to the values 2.35214, 2.23869, 2.04516, 2.43269 and 1.80415 in
order.

Citation of OpenKIM IMs
-----------------------

When publishing results obtained using OpenKIM IMs researchers are requested
to cite the OpenKIM project :ref:`(Tadmor) <kim-mainpaper>`, KIM API
:ref:`(Elliott) <kim-api>`, and the specific IM codes used in the simulations,
in addition to the relevant scientific references for the IM.
The citation format for an IM is displayed on its page on
`OpenKIM <https://openkim.org>`_ along with the corresponding BibTex file,
and is automatically added to the LAMMPS *log.cite* file.

Citing the IM software (KIM infrastructure and specific PM or SM codes)
used in the simulation gives credit to the researchers who developed them
and enables open source efforts like OpenKIM to function.

Restrictions
""""""""""""

The set of *kim_commands* is part of the KIM package.  It is only enabled if
LAMMPS is built with that package. A requirement for the KIM package,
is the KIM API library that must be downloaded from the
`OpenKIM website <https://openkim.org/kim-api/>`_ and installed before
LAMMPS is compiled. When installing LAMMPS from binary, the kim-api package
is a dependency that is automatically downloaded and installed. See the KIM
section of the :doc:`Packages details <Packages_details>` for details.

Furthermore, when using *kim_commands* to run KIM SMs, any packages required
by the native potential being used or other commands or fixes that it invokes
must be installed.

Related commands
""""""""""""""""

:doc:`pair_style kim <pair_kim>`

----------

.. _kim-mainpaper:

**(Tadmor)** Tadmor, Elliott, Sethna, Miller and Becker, JOM, 63, 17 (2011).
doi: `https://doi.org/10.1007/s11837-011-0102-6 <https://doi.org/10.1007/s11837-011-0102-6>`_

.. _kim-api:

**(Elliott)** Elliott, Tadmor and Bernstein, `https://openkim.org/kim-api <https://openkim.org/kim-api>`_ (2011)
doi: `https://doi.org/10.25950/FF8F563A <https://doi.org/10.25950/FF8F563A>`_
