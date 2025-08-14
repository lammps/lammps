.. index:: fix rmc/partial

fix rmc/partial command
=========================

Syntax
"""""""

.. code-block:: LAMMPS

   fix ID group-ID rmc/partial Nmd Nmc Nad Nas Nd Ns Temp TypeThreshold Ncs restart

* ID, group-ID = fix ID and group ID (as defined in :doc:`fix <fix>`)
* Nmd = number of MD steps in each Reactive MCMD cycle
* Nmc = number of Reactive MC steps in each Reactive MCMD cycle
* Nad = number of atoms in a dopant molecule
* Nas = number of atoms in a semiconductor molecule
* Nd = number of dopant molecules
* Ns = number of semiconductor molecules
* Temp = temperature (Kelvin)
* TypeThreshold = threshold for atom type, atom with type <= TypeThreshold are considered semiconductor atoms, otherwise dopant atoms
* Ncs = number of charge states, choose 2 for integer charge states, or more to include partial charge
* restart = 0 if fresh run, 1 if continuing from previous run. Requires necessary restart files (discussed below)

Examples
"""""""""

.. code-block:: LAMMPS
    
    fix 1 all nvt temp 300 300 100
    fix 2 all rmc/partial 1000 200 20 202 25 100 300 16 2 0  
    # Fresh RMCMD run at 300 K, 2 charge states (0 and +-1)
    # 25 dopant molecules, 100 semiconductor molecules
    # 20 atoms in each dopant molecule, 202 atoms in each semiconductor molecule
    # 1000 MD steps, 200 Reactive MC steps per cycle
    # Atoms of type 16 or lesser belong to semiconductor molecules, others belong to dopant molecules

    fix 2 all rmc/partial 1000 200 20 202 25 100 300 16 6 1
    # Continuing RMCMD run at 300K, 6 partial charge states (0, +-0.2, +-0.4, +-0.6, +-0.8, +-1.0)
    # Other parameters same as above

Files required
"""""""""""""""""""""
In addition to the above command, this fix also requires the following files:

.. parsed-literal::
    
    **semiconductor_neutral.dat**
    `file containing the charges of atoms in a single neutral semiconductor molecule, with one value per line. For eg.`
    0.3924909999999989
    -0.177266
    -0.073638
    -0.080285
    -0.172156
    0.048395
    0.0015773157894736843
    ...
    The order of the charges should match the order of the atoms in the molecule, when arranged from lowest to highest atom ID. This means that the first charge
    corresponds to the first atom in the semiconductor molecule, denoted by the lowest atom ID of all the atoms in that molecule. Note that this order should be
    the same for each semiconductor molecule, i.e if the lowest atom ID in a molecule is subtracted from all the atoms in that molecule, the resulting "local" 
    atom ID should refer to the same atom in each molecule. Additionally, the atom IDs within each molecule must also be contiguous. 

    **semiconductor_charged.dat**
    `file containing the charges of atoms in a single ionized semiconductor molecule (+1 charge), with one value per line. For eg.`
    0.42590700000000004
    -0.196852
    -0.017455
    -0.068876
    -0.118641
    0.06525
    0.0022099473684210526
    The order of the charges should match the order of the atoms in the molecule (same logic as semiconductor_neutral.dat file).

    **dopant_neutral.dat**
    `file containing the charges of atoms in a single neutral dopant molecule, with one value per line. For eg.`
    -0.0043649999999999705
    -0.023369
    -0.067212
    -0.023214
    -0.001734
    -0.005162
    ...
    The order of the charges should match the order of the atoms (same logic as semiconductor_neutral.dat file).

    **dopant_charged.dat**
    `file containing the charges of atoms in a single charged dopant molecule (-1 charge), with one value per line. For eg,`
    -0.109287
    -0.023825
    -0.187472
    -0.029575
    -0.118920
    -0.040151
    0.039783
    ...
    The order of the charges should match the order of the atoms in the molecule (same logic as semiconductor_neutral.dat file).

    **dihedral_list.dat**
    `A list of semiconductor dihedrals whose coefficients will be modified during an Reactive Monte Carlo step. The first line shows number of dihedrals to consider,
    the second line should have a dihedral type for each charge state, and subsequent line should have the 4 atom IDs for each dihedral. For eg.`
    700
    6 10
    13535 13536 13564 13560
    2323 2324 2352 2348
    2373 2374 2402 2398
    2348 2349 2377 2373
    ...
    In this example, the first line indicates that there are 700 dihedrals to consider. The second line indicates that the first dihedral type (6) will be used
    for charge state 0, and the second dihedral type (10) will be used for charge state +1. The subsequent lines contain the atom IDs of the atoms in each
    dihedral. The coefficients for dihedral type 6 and type 10 are specified in the input data file.
    
    This straightforwardly expands to the number of charge states you have. For example, if you have 6 charge states, you can specify 6 dihedral types in the
    second line, and then list the atom IDs for each dihedral in subsequent lines.
    700
    6 7 8 9 10 11
    13535 13536 13564 13560
    2323 2324 2352 2348
    2373 2374 2402 2398
    2348 2349 2377 2373
    ...

    If you do not wish to consider dihedral effects, you can create a dihedral_list.dat file with one dihedral and the same type for both states i.e
    1
    6 6
    13535 13536 13564 13560

    **angle_list.dat**
    `A list of semiconductor angles whose coefficients will be modified during an Reactive Monte Carlo step. The first line shows number of angles to consider,
    the second line should have an angle type for each charge state, and subsequent line should have the 3 atom IDs for each angle. For eg.`
    1
    10 10
    1 2 30

    Identical to the dihedral_list.dat file, except with three IDS representing an angle.

    **bond_list.dat**
    `A list of semiconductor bonds whose coefficients will be modified during an Reactive Monte Carlo step. The first line shows number of bonds to consider,
    the second line should have a bond type for each charge state, and subsequent line should have the 2 atom IDs for each bond. For eg.`
    1
    4 4
    965 967

    Identical to the dihedral_list.dat file, except with two IDS representing a bond.

    **barrier_list.dat**
    `A list of the reaction free energies for each charge state. The energies should be written in a single line. For eg.`
    17.304428 13.853026 15.157749 21.218598 32.035572 47.608672

    Here there are six values for six charge states. For a move from 0 to +-1 charge state, the energy barrier used in the MC acceptance step is 
    (47.608672 - 17.304428) = 30.304244 kcal/mol. A much simpler case, with only 0/+1 charge states

    0 30
    would have only two values, one for neutral and one for the ionized state.

Description
"""""""""""

This fix performs Reactive Monte Carlo Molecular Dynamics (RMCMD) simulations for the simultaneous treatment of the doping reaction and morphology evolution
in doped organic semiconductors, based on :ref:`Verma2024<Verma2024>` and :ref:`Raghuraman2025<Raghuraman2025>`. The method alternates between MD steps (for the morphology) and MC steps
where a neutral semiconductor and dopant molecule is chosen randomly, and the atomic charges are replaced to ionize the molecule. The fix is capable of
capturing Integer Charge Transfer (ICT), where a molecule goes from neutral to fully ionized, and Charge Transfer Complexes (CTC), where a molecule can
partially ionize. In addition to switching out charges, the fix can also switch out dihedral types, angle types and bond types of the chosen molecules, enabling
the study of intramolecular conformational effects as a function of doping. The fix should used be alongside an MD fix like nvt or npt.

.. image:: JPG/rmcmd-workflow.png
   :align: center
   :scale: 10%

The figure above shows a schematic of the RMCMD cycle, using NVT MD as the fix for morphology.

A step-by-step guide on running RMCMD will be made available shortly at `thejacksonlab.github.io website <https://thejacksonlab.github.io>`_.

Restrictions
"""""""""""""
This fix requires the LAMMPS package :ref:`MC <PKG-MC>` to be built.

______________

.. _Verma2024:

**(Verma2024)** Verma, A.; Jackson, N. E. Assessing molecular doping efficiency in organic semiconductors with reactive Monte Carlo, Journal of Chemical Physics 
2024, 160 

.. _Raghuraman2025:

**(Raghuraman2025)**  Raghuraman V, Verma A, Jackson N. All-Atom Reactive Monte Carlo Molecular Dynamics for Molecular Doping in Organic Semiconductors. ChemRxiv. 2025; doi:10.26434/chemrxiv-2025-713hw-v2  This content is a preprint and has not been peer-reviewed.

