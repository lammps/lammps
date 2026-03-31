.. index:: pair_style mbx

pair_style mbx command
======================

Syntax
""""""

.. code-block:: LAMMPS

    pair_style mbx cutoff

* cutoff = real-space cutoff for MBX in Angstroms


Examples
""""""""

.. code-block:: LAMMPS

    pair_style      mbx 9.0
    pair_coeff      * * 1 h2o 1 2 2 json mbx.json
    compute         mbx all pair mbx

    # For a system involving ch4 (atom types C=1, H=2) and
    # water (atom types O=3, H=4)
    pair_style      mbx 9.0
    pair_coeff      * * 2 ch4 1 2 2 2 2 h2o 3 4 4 json mbx.json
    compute         mbx all pair mbx

    # For a system involving water (atom types O=12, H=13) in a hybrid simulation
    pair_style      hybrid/overlay mbx 9.0 lj/cut 9.0 coul/exclude 9.0
    pair_coeff      * * mbx 2 dp1 1*11 h2o 12 13 13 json mbx.json
    pair_coeff      1*11 1*11 coul/exclude
    compute         mbx all pair mbx

    # For a system involving water (atom types O=12, H=13) in a hybrid simulation
    # with special_bonds and coul/exclude to exclude 1-2, 1-3, and 1-4 electrostatics
    # for the charmm framework
    special_bonds   charmm
    pair_style      hybrid/overlay mbx 9.0 lj/cut 9.0 coul/exclude 9.0
    pair_coeff      * * mbx 2 dp1 1*11 h2o 12 13 13 json mbx.json
    pair_coeff      1*11 1*11 coul/exclude
    compute         mbx all pair mbx

See ``examples/PACKAGES/mbx`` for full examples of how to use MBX in LAMMPS.


Description
"""""""""""

.. versionadded:: 11Feb2026

The MBX (Many-Body eXpansion) software is a C++ library that provides
access to many-body energy (MB-nrg) potential energy functions, such as
the MB-pol water model.  Developed over the past decade, these potential
energy functions integrate physics-based and machine-learned many-body
terms trained on electronic structure data calculated at the "gold
standard" coupled-cluster level of theory. :ref:`(Gupta) <Gupta3>`

This pair_style instructs LAMMPS to call the
`MBX library <https://mbxsimulations.com>`_ in order to simulate
MB-nrg models such as MB-pol. The MBX library source code is available at
`https://github.com/paesanilab/MBX <https://github.com/paesanilab/MBX>`_.
MBX is heavily OpenMP parallelized (OMP), and the OMP_NUM_THREADS
environment variable should be properly set to the number of threads desired.
A detailed discussion of the code structure can be found in the
manuscript :ref:`(Riera) <Riera>`, while a detailed description of the
performance scaling can be found in the manuscript :ref:`(Gupta) <Gupta3>`.

The *cutoff* argument specifies the real-space cutoff for MBX in
Angstroms. This real-space cutoff is used for the dispersion interactions of the
MB-nrg monomers, as well as for the electrostatics of the **entire** system.
For periodic systems, a safe value for the real-space cutoff is **9.0 Angstroms**,
and all classical interactions beyond this cutoff will then be handled via particle-mesh
Ewald (PME) within MBX. For non-periodic systems, the cutoff can be set to a
large value, such as 100.0 Angstroms, to ensure that all interactions are
captured in the real-space.

.. warning::

    MBX must currently be used with :doc:`processors <processors>` mapping style `xyz`. If you
    do not, MBX will throw the error:

    .. code-block:: LAMMPS

        [MBX] Inconsistent proc mapping: 'processors * * * map xyz' required for PME solver

For hybrid simulations involving MB-nrg and non-MB-nrg molecules in the
same simulation, one can use :doc:`pair_style hybrid/overlay
<pair_hybrid>` to combine the MB-nrg molecules with other pair styles,
such as :doc:`lj/cut <pair_lj>`. This has been used to simulate
MB-pol water within host frameworks such as metal-organic
frameworks (MOFs) and carbon nanotubes (CNTs).
If using MBX in a hybrid simulation involving :doc:`special_bonds <special_bonds>`,
(such as when using the CHARMM, Amber, OPLS, or ClayFF force fields etc.),
please see the warning below for more details about
using `special_bonds` with MBX `dp1`. See ``examples/PACKAGES/mbx`` for
a complete hybrid example.


If you have questions not answered by this documentation, please
reference the MBX website
`mbxsimulations.com <https://mbxsimulations.com>`_ or reach out to the MBX team at
`https://groups.google.com/g/mbx-users <https://groups.google.com/g/mbx-users>`_


Pair coeff syntax
"""""""""""""""""

MBX is many-body method, and only a single pair_coeff command is needed
to specify the mapping of LAMMPS atom IDs to MBX monomers. The syntax is as follows:

.. code-block:: LAMMPS

    pair_coeff * * num_mon_types mon_name atom_mapping <mon_name2> <atom_mapping2> ... json mbx.json print/settings

* num_mon_types = number of monomer types in the system
* mon_name = name of the monomer type (e.g. h2o, ch4, etc)
* atom mapping = list of LAMMPS atom IDs that correspond to the atoms in the monomer
* *json* arg = specifies the name of the MBX json configuration file, such as mbx.json
* print/settings = optionally print MBX settings to logfile


The *num_mon_types* argument specifies the number of different MB-nrg
monomer types in the system.

For each monomer type, the *mon_name* argument specifies the name of
the monomer, such as `h2o` for water or `ch4` for methane. The *atom
mapping* argument specifies then the mapping of LAMMPS atom IDs to
the atoms in the monomer, such as `1 2 2` for water (O=1, H=2).
For hybrid simulations, the `dp1` (drude particle) monomer
should be used to represent the non-MB-nrg molecules. `dp1` is a
special monomer in MBX in that its *atom_mapping* can be a range of
LAMMPS atom IDs, such as `1*11` to represent atom types 1 through 11.
For a complete list of available monomers in MBX, please see the
`MBX documentation <https://mbxsimulations.com/tutorials/pefs-implemented-and-how-to-cite>`_.

.. warning::

    When using MBX, **all electrostatics are handled internally by MBX.**
    This is important since MB-nrg models such as MB-pol are polarizable
    models that may also use geometrically dependent charges, such as
    the Partridge and Schwenke charges used in MB-pol water.

    Therefore, one should **never** use a coulombic pair style in LAMMPS
    such as `coul/cut` or `coul/long` when also using MBX. This
    would result in double counting of electrostatic interactions.

    When performing a hybrid simulation using dp1,
    note that many frameworks (Amber, CHARMM, OPLS, ClayFF etc.) require the usage of
    special_bonds to exclude some bonded coulomb interactions (1-2, 1-3, and/or 1-4).
    Since MBX is handling all electrostatics, this should therefore be accounted for
    using the :doc:`coul/exclude <pair_coul>` command.

    .. code-block:: LAMMPS

        # For a system involving water (atom types O=12, H=13) in a hybrid simulation
        # with special_bonds and coul/exclude to exclude 1-2, 1-3, and 1-4 electrostatics
        # for the charmm framework
        processors      * * * map xyz
        special_bonds   charmm
        pair_style      hybrid/overlay mbx 9.0 lj/cut 9.0 coul/exclude 9.0
        pair_coeff      * * mbx 2 dp1 1*11 h2o 12 13 13 json mbx.json
        pair_coeff      1*11 1*11 coul/exclude
        compute         mbx all pair mbx

The *json* argument specifies the name of the MBX JSON configuration
file to use, such as `mbx.json`.  If this file is not provided, the pair style
will attempt to use a default configuration. See the `mbx.json documentation
<https://mbxsimulations.com/tutorials/json-file>`_ for more details on
how to create this file.

The *print/settings* argument optionally will print the MBX settings to the LAMMPS
logfile at the start of the simulation. This is optionally used for debugging and
ensuring that the settings are being correctly applied.


Restrictions
""""""""""""

This pair_style is part of the MBX package.  A pair style is only
enabled if LAMMPS was built with its corresponding package.
See the :doc:`Build package <Build_package>` page for more info.

MBX requires the FFTW3 library to be installed. This is needed
as part of the internal PME solver used for long-range electrostatics.

All electrostatic interactions are calculated internally in MBX.
Therefore one should never calculate coulombic interactions in
LAMMPS such as using `coul/cut` or `coul/long` when also using MBX.
See the warning above for more details.

MBX currently only supports :doc:`processors <processors>` mapping style `xyz`.

MBX is primarily tested to work with `units real` and `atom_style full`. If you encounter
issues with other unit or atom styles, please contact the MBX developers.

Related commands
""""""""""""""""

:doc:`pair hybrid/overlay <pair_hybrid>`,
:doc:`pair coul/exclude <pair_coul>`

-----------

.. _Riera:

**(Riera)** M. Riera, C. Knight, E. Bull-Vulpe, X. Zhu, H. Agnew, D. Smith, A. Simmonett, F. Paesani, J. Chem. Phys. 159, 054802 (2023)

.. _Gupta3:

**(Gupta)** S. Gupta, E. Bull-Vulpe, H. Agnew, S. Iyer, X. Zhu, R. Zhou, C. Knight, F. Paesani, J. Chem. Theory Comput. 21, 1938 (2025)

