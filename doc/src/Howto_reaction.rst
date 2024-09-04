
Chemical reactions
==================

The five basic types of chemical reactions are combination, decomposition, single-replacement, double-replacement, and combustion :ref:`(LibreTexts) <howto-reaction-LibreTexts>`. A combination reaction, also known as a synthesis reaction, is a reaction in which two or more substances combine to form a single new substance. One combination reaction is two elements combining to form a compound. The general form of a combination reaction is

.. math:: A + B \rightarrow AB

For example, solid sodium metal reacts with chlorine gas to product solid sodium chloride

.. math:: 2 Na + Cl_2 \rightarrow 2 NaCl



Elementary reaction
^^^^^^^^^^^^^^^^^^^

An *elementary reaction* is a single step reaction with a single transition state and no intermediates. Elementary reactions add up to complex reactions; non-elementary reactions can be described by multiple elementary reaction steps. A set of elementary reactions comprises a reaction mechanism.

Unimolecular reaction
"""""""""""""""""""""

An elementary *unimolecular reaction* occurs when a molecule rearranges itself to produce one or more products

Bimolecular reaction
""""""""""""""""""""

An elementary *bimolecular reaction* involves the collision of two particles. Bimolecular reactions are common in organic reactions such as nucleophilic substitution. There are two types of bimolecular elementary reactions, either the two reactants are the same:

.. math::
  2 A  \rightarrow  products

or the two reactants are different:

.. math::
  A + B  \rightarrow  products

.. admonition:: Example: gas-phase (g) reaction of 1,3-butadiene with ethene to form cyclohexene :ref:`(Moore, Zhou, Garand) <howto-reaction-moore-zhou-garand>`
  :class: Hint

  A suggested workflow for using `fix bond/react` is:

  **(1) identify a reaction to be simulated**

  .. math::
    C_4H_6(g) + C_2H_4(g) \rightarrow C_6H_{10}(g)

  .. image:: img/diels_alder_pes.png
    :align: center
    :width: 61.8%

  In the transition state, bonds that are breaking (three π bonds) are shown in orange and bonds that are forming (one π bond and two σ bonds) are shown in cyan.

  **(2) build a molecule template of the reaction site before the reaction has occurred**

  The pre-reacted molecule template is specified by a molecule command. This molecule template file contains a sample reaction site and its surrounding topology. All atom types in the pre-reacted template must be the same as those of a potential reaction site in the simulation. The initiator atom pairs of the pre-reacted template are specified by atom ID in the map file.

  **(3) build a molecule template of the reaction site after the reaction has occurred**

  The post-reacted molecule template contains a sample of the reaction site and its surrounding topology after the reaction has occurred. It must contain the same number of atoms as the pre-reacted template, unless there are created or deleted atoms (see examples/PACKAGES/reaction for details). A one-to-one correspondence between the atom IDs in the pre- and post-reacted templates is specified in the map file described below (4).













Termolecular reaction
"""""""""""""""""""""

An elementary *termolecular reaction* requires the collision of three particles at the same place and time. This type of reaction is extremely rare because all three reactants must simultaneously collide with each other, with sufficient energy and correct orientation, to produce a reaction. When a reaction involves three reactant molecules, it is much more likely for it to proceed via multiple steps known as a *reaction mechanism* involving elementary unimolecular and/or bimolecular reaction steps.

Reaction mechanism
^^^^^^^^^^^^^^^^^^

A valid multi-step *reaction mechanism* consist of a series of unimolecular and/or bimolecular elementary reaction steps. The sum of the reaction steps should agree with the overall balanced reaction equation. A reaction intermediate is transient species within a multi-step reaction mechanism that is produced in the preceding step and consumed in a subsequent step to ultimately generate the final reaction product. Intermediate reactions are common in the biological world; a prime example can be seen in the metabolism of metabolites and nutrients.



Combustion reaction
"""""""""""""""""""

A *combustion reaction*, which also qualifies as a combination reaction, is a reaction in which a substance reacts with oxygen gas, releasing energy in the form of light and heat. Combustion reactions must involve O\ :sub:`2`\  as one reactant. The combustion of hydrogen gas produces water vapor:

.. math:: 2 H_2 + O_2 \rightarrow 2 H_2O

However this combustion reaction is actually more complicated than 2 diatomic hydrogen molecules simply colliding with a diatomic oxygen molecule to form 2 water molecules. It actually proceeds as a series of steps in the reaction mechanism:

.. math::
  H + O_2 \rightarrow O + OH \\
  H_2 + O \rightarrow H + OH \\
  H_2 + OH \rightarrow H_2O + H \\

..  youtube:: YuqA_uojSJ4
  :align: center

Oxyhydrogen is a mixture of hydrogen (H2) and oxygen (O2) gases, also known as "Knallgas" (from German lit. "bang-gas"). Theoretically, a ratio of 2:1 hydrogen:oxygen is enough to achieve maximum efficiency; in practice a ratio 4:1 or 5:1 is needed to avoid an oxidizing flame. Oxyhydrogen will combust when brought to its autoignition temperature of 843 K at normal atmospheric pressure. The minimum energy required to ignite such a mixture, at lower temperatures, with a spark is about 20 microjoules. At standard temperature and pressure, oxyhydrogen can burn when it is between about 4% and 95% hydrogen by volume.

When ignited, the gas mixture converts to water vapor and releases energy, which sustains the reaction: 241.8 kJ of energy (LHV) for every mole of H2 burned. The amount of heat energy released is independent of the mode of combustion, but the temperature of the flame varies. The maximum temperature of about 3073 K is achieved with an exact stoichiometric 2:1 mixture.

.. admonition:: Example: Hydrogen combustion (REAXFF version)
  :class: Hint

  .. code-block:: LAMMPS
    :caption: `examples/reaxff/hydrogen_combustion/in.hydrogen_combustion`

    units real
    dimension 3
    boundary p p p
    atom_style charge
    newton on

    region box block 0 99 0 99 0 99
    create_box 2 box

    mass 1 1.008
    mass 2 15.999

    create_atoms 1 random 2000 12345 NULL overlap 0.31 maxtry 100
    create_atoms 2 random 1000 12345 NULL overlap 0.66 maxtry 100
    velocity all create 2500 12345

    pair_style      reaxff lmp_control
    pair_coeff      * * ffield.reax.cho H O

    neighbor        2 bin
    neigh_modify    every 10 delay 0 check no

    variable        dt equal 0.1
    timestep        ${dt}
    fix             1 all nvt temp 2500 2500 $(100.0*dt)
    fix             2 all qeq/reax 1 0.0 10.0 1e-6 param.qeq
    fix             4 all reaxff/species 10 10 100 species.out

    thermo          100
    dump            1 all movie 10 hydrogen_combustion.mkv type type size 800 800
    dump_modify     1 acolor * white/red/green/blue/aqua/magenta

    run             $(500*1000/dt)



----------

.. admonition:: Example: Hydrogen combustion (REACTION version)
  :class: Hint

  A suggested workflow for using `fix bond/react` is:

  **(1) identify a reaction to be simulated**
  


Using standard reference data on hydrogen, oxygen, and water :ref:`(NIST WebBook) <howto-reaction-nist-webbook>`

.. parsed-literal::

    Hydrogen, ID: C1333740
      NIST    23112819112D 1   1.00000     0.00000
    Copyright by the U.S. Sec. Commerce on behalf of U.S.A. All rights reserved.
      2  1  0     0  0              1 V2000
        0.0000    0.0000    0.0000 H   0  0  0  0  0  0           0  0  0
        1.0000    0.0000    0.0000 H   0  0  0  0  0  0           0  0  0
      1  2  1  0     0  0
    M  END

.. parsed-literal::

    Oxygen, ID: C7782447
      NIST    23112607342D 1   1.00000     0.00000
    Copyright by the U.S. Sec. Commerce on behalf of U.S.A. All rights reserved.
      2  1  0     0  0              1 V2000
        0.0000    0.0000    0.0000 O   0  0  0  0  0  0           0  0  0
        1.0000    0.0000    0.0000 O   0  0  0  0  0  0           0  0  0
      1  2  2  0     0  0
    M  END

.. parsed-literal::

    Water, ID: C7732185
      NIST    23112607382D 1   1.00000     0.00000
    Copyright by the U.S. Sec. Commerce on behalf of U.S.A. All rights reserved.
      3  2  0     0  0              1 V2000
        0.0000    0.0000    0.0000 H   0  0  0  0  0  0           0  0  0
        0.8638    0.4573    0.0000 O   0  0  0  0  0  0           0  0  0
        1.7785    0.0000    0.0000 H   0  0  0  0  0  0           0  0  0
      1  2  1  0     0  0
      2  3  1  0     0  0
    M  END

we can write the molecule template files:

.. code-block::

    # REACTION PACKAGE COMBUSTION EXAMPLE - Hydrogen molecule template (combustion_H2.molecule_template)

        2 atoms
        1 bonds
        2.01588 mass
            
    Coords

        1    0.0000    0.0000    0.0000
        2    1.0000    0.0000    0.0000

    Types

        1 H
        2 H

    Bonds

        1 1 1 2

.. code-block::

    # REACTION PACKAGE COMBUSTION EXAMPLE - Oxygen molecule template (combustion_O2.molecule_template)

        2 atoms
        1 bonds
        31.9988 mass

    Coords

        1    0.0000    0.0000    0.0000
        2    1.0000    0.0000    0.0000

    Types

        1 O
        2 O

    Bonds

        1 2 1 2


.. code-block::

    # REACTION PACKAGE COMBUSTION EXAMPLE - Water molecule template (combustion_H2O.molecule_template)

        3 atoms
        2 bonds
        1 angles
        18.0153 mass

    Coords

        1    0.0000    0.0000    0.0000
        2    0.8638    0.4573    0.0000
        3    1.7785    0.0000    0.0000

    Types

        1 H
        2 O
        3 H

    Bonds

        1 1 1 2
        2 1 2 3

    Angles

        1 1 1 2 3


**(2) build a molecule template of the reaction site before the reaction has occurred**

The pre-reacted molecule template is specified by a molecule command. This molecule template file contains a sample reaction site and its surrounding topology. All atom types in the pre-reacted template must be the same as those of a potential reaction site in the simulation. The initiator atom pairs of the pre-reacted template are specified by atom ID in the map file.


.. code-block::

    REACTION PACKAGE COMBUSTION EXAMPLE - molecule template pre-reaction 2 H2 and O2 (combustion_pre.molecule_template)

    6 atoms
    3 bonds

    Coords

        1    0.0000    0.0000    0.0000
        2    1.0000    0.0000    0.0000
        3    ???      ???     ???
        4    ???      ???     ???
        5    ???      ???     ???
        6    ???      ???     ???

    Types

        1 H
        2 H
        3 H
        4 H
        5 O
        6 O

    Bonds

        1 H-H      1      2
        2 H-H      3      4
        3 O-O      5      6


**(3) build a molecule template of the reaction site after the reaction has occurred**

The post-reacted molecule template contains a sample of the reaction site and its surrounding topology after the reaction has occurred. It must contain the same number of atoms as the pre-reacted template, unless there are created or deleted atoms (see examples/PACKAGES/reaction for details). A one-to-one correspondence between the atom IDs in the pre- and post-reacted templates is specified in the map file described below (4).


.. code-block::

    REACTION PACKAGE COMBUSTION EXAMPLE - molecule template post-reaction 2 H2O (combustion_post.molecule_template)

    6 atoms
    4 bonds
    2 angles

    Coords

        1    0.0000    0.0000    0.0000
        2    0.8638    0.4573    0.0000
        3    1.7785    0.0000    0.0000
        4    ???    ???    ???
        5    ???    ???    ???
        6    ???    ???    ???

    Types

        1 H
        2 O
        3 H
        4 H
        5 O
        6 H

    Molecules

        1      1
        2      1
        3      1
        4      2
        5      2
        6      2

    Bonds

        1 H-O      1      2
        2 O-H      2      3
        3 H-O      4      5
        4 O-H      5      6

    Angles

        1 H-O-H    1   2   3
        2 H-O-H    4   5   6



**(4) create a map that relates the template-atom-IDs of each atom between pre- and post-reaction molecule templates**

The header of map file contains one mandatory keyword *equivalences*\, which is the number of atoms in the pre- and post-reaction  molecule templates.

The body of the map file contains two mandatory sections. The first mandatory section begins with the keyword *InitiatorIDs*\  listing the two atom IDs of the initiator atom pair in the pre-reacted molecule template. The second mandatory section begins with the keyword *Equivalences*\  listing a one-to-one correspondence between atom IDs of the pre- and post-reacted templates. The first column is an atom ID of the pre-reacted molecule template, and the second column is the corresponding atom ID of the post-reacted molecule template.

.. parsed-literal::

    REACTION PACKAGE COMBUSTION EXAMPLE - map file (combustion.map)

    6 equivalences

    InitiatorIDs

        ???
        ???

    Equivalences

        1  1
        2  3
        3  4
        4  5
        5  2
        6  6

   
**(5) fill a simulation box with molecules and run a simulation with fix bond/react**

.. code-block:: LAMMPS

    # REACTION PACKAGE COMBUSTION EXAMPLE - input script (combustion.in)

    units real
    dimension 3
    boundary p p p
    atom_style full

    region combustion_region block 0.0 10.0 0.0 10.0 0.0 10.0
    create_box 2 combustion_region bond/types 2 angle/types 1 extra/special/per/atom 2

    labelmap atom 1 H 2 O
    mass H 1.008
    mass O 15.999
    molecule H2 combustion_H2.molecule_template
    molecule O2 combustion_O2.molecule_template
    molecule H2O combustion_H2O.molecule_template

    create_atoms 1 random 20 12345 NULL overlap 2.0 maxtry 50
    create_atoms 2 random 10 12345 NULL overlap 2.0 maxtry 50
    velocity all create 310.0 12345

    pair_style lj/cut 2.5
    pair_coeff * * 1.0 1.0 2.5
    
    fix combustion_fix1 all langevin 310.0 310.0 1000 12345
    fix combustion_fix2 all nve

    dump combustion_movie all movie 1 combustion.mpg type type size 512 512
    #dump modify combustion_movie acolor H white
    #dump modify combustion_movie acolor O red

    #molecule combustion_pre combustion_pre.molecule_template
    #molecule combustion_post combustion_post.molecule_template
    #fix combustion_map all bond/react react myrxn1 all 1 0 3.25 combustion_pre combustion_post combustion.map

    run 1000

    undump combustion_movie




Enzyme-substrate reaction
"""""""""""""""""""""""""

A reaction mechanism found in all living systems is the *enzyme-substrate reaction*. In this type of reaction, an enzyme binds to a substrate to produce an enzyme-substrate intermediate, which then forms the final product.

An example is *Glucose-6-phosphate isomerase (GPI)*, an enzyme (EC 5.3.1.9) that converts *glucose-6-phosphate (G6P)* to *fructose-6-phosphate (F6P)* as part of the glycolysis pathway. Since the reaction is reversible, its direction is determined by G6P and F6P concentrations. The mechanism that GPI uses to interconvert glucose 6-phosphate and fructose 6-phosphate consists of three major steps: opening the glucose ring, isomerizing glucose into fructose through an enediol intermediate, and closing the fructose ring. Human GPI pdb 1JLH



























..
  Some chemical reactions have mechanisms that consist of a single bimolecular elementary reaction.


..
  Many reactions have at least one activation energy that must be reached in order for the reaction to go forward.

..
  Chain reactions usually consist of many repeating elementary steps, each of which has a chain carrier. Once started, chain reactions continue until the reactants are exhausted. Fire and explosions are some of the phenomena associated with chain reactions. The chain carriers are some intermediates that appear in the repeating elementary steps. These are usually free radicals.

..
  Once initiated, repeating elementary steps continue until the reactants are exhausted.

..
  Chain Branching Steps

..
  Branching reactions are elementary steps that generate more free radicals than they consume. Branching reactions result in an explosion. For example, in the reaction between hydrogen and oxygen, the following reaction may take place: H⋅+O2→HO⋅+⋅O⋅H⋅+O2→HO⋅+⋅O⋅ where ⋅O⋅⋅O⋅ is a di-radical, because the OO atom has an electronic configuration 2s2 2px2 2py1 2pz1. In this elementary step, three radicals are generated, whereas only one is consumed. The di-radical may react with a H2 molecule to form two radicals. ⋅O⋅+H2→HO⋅+H⋅⋅O⋅+H2→HO⋅+H⋅

..
  Thus, together chain branching reactions increase the number of chain carriers. Branching reactions contribute to the rapid explosion of hydrogen-oxygen mixtures, especially if the mixtures have proper proportions.



----------

REACTION package
----------------


This package implements the REACTER protocol (reacter.org) as :doc:`fix bond/react <fix_bond_react>`. This fix allows for complex topology changes during a running MD simulation, when using classical force fields. Topology changes are defined in pre- and post-reaction molecule templates and can include creation and deletion of bonds, angles, dihedrals, impropers, atom types, bond types, angle types, dihedral types, improper types, and/or atomic charges.

The REACTER protocol is a method for modeling chemical reactions in classical molecular dynamics simulations. It was developed to build physically-realistic initial configurations for amorphous or crosslinked materials. Any number of competing or reversible reaction pathways can be specified, and reacting sites can be stabilized. Other advanced options currently available include reaction constraints (e.g. angle and Arrhenius constraints), deletion of reaction byproducts or other small molecules, creation of new atoms or molecules bonded to existing atoms, and using LAMMPS variables for input parameters.

The REACTER methodology is detailed in: :ref:`(Gissinger, 2017) <howto-reaction-Gissinger-2017>` and :ref:`(Gissinger, 2020) <howto-reaction-Gissinger-2020>`. This package was created by Jacob Gissinger at the NASA Langley Research Center.






----------

.. _howto-reaction-Gissinger-2017:

**(Gissinger, 2017)** Gissinger et al., Modeling chemical reactions in classical molecular dynamics simulations. Polymer 128, 211-217 (2017)
    https://doi.org/10.1016/j.polymer.2017.09.038

.. _howto-reaction-Gissinger-2020:

**(Gissinger, 2020)** Gissinger et al., REACTER: A Heuristic Method for Reactive Molecular Dynamics. Macromolecules 53, 22, 9953-9961 (2020).
    https://doi.org/10.1021/acs.macromol.0c02012

.. _howto-reaction-LibreTexts:

**(LibreTexts)** Valley City State University, Chem 121.
    https://chem.libretexts.org/Courses/Valley_City_State_University/Chem_121/Chapter_5%3A_Introduction_to_Redox_Chemistry/5.3%3A_Types_of_Chemical_Reactions

.. _howto-reaction-nist-webbook:

**(NIST WebBook)**
    Hydrogen: https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2/h1H
    
    Oxygen: https://webbook.nist.gov/cgi/inchi/InChI%3D1S/O2/c1-2
    
    Water: https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O/h1H2

.. _howto-reaction-moore-zhou-garand:

**(Moore, Zhou, Garand)** https://chem.libretexts.org/Bookshelves/General_Chemistry/Interactive_Chemistry_(Moore_Zhou_and_Garand)/03%3A_Unit_Three/3.05%3A_Day_22-_Elementary_Reactions
