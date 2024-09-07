
Chemical reactions
==================

The five basic types of chemical reactions are combination, decomposition, single-replacement, double-replacement, and combustion :ref:`(LibreTexts) <howto-reaction-LibreTexts>`.

  .. image:: img/chemical_reactions.svg
    :align: center
    :width: 61.8%

A combination reaction, also known as a synthesis reaction, is a reaction in which two or more substances combine to form a single new substance. One combination reaction is two elements combining to form a compound. The general form of a combination reaction is

.. math:: A + B \rightarrow AB

For example, solid sodium metal reacts with chlorine gas to product solid sodium chloride

.. math:: 2 Na + Cl_2 \rightarrow 2 NaCl


----------

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

.. admonition:: REACTION example: Diels–Alder reaction


  :class: Hint

  :doc:`fix bond/react <fix_bond_react>` allows for complex topology changes during a running MD simulation, when using classical force fields. Topology changes are defined in pre- and post-reaction molecule templates and can include creation and deletion of bonds, angles, dihedrals, impropers, atom types, bond types, angle types, dihedral types, improper types, and/or atomic charges. :ref:`(Gissinger, 2017) <howto-reaction-Gissinger-2017>` and :ref:`(Gissinger, 2020) <howto-reaction-Gissinger-2020>`. A suggested workflow is:

  **(1) identify a reaction to be simulated**

  gas-phase (g) reaction of 1,3-butadiene with ethene to form cyclohexene :ref:`(Moore, Zhou, Garand) <howto-reaction-moore-zhou-garand>`

  .. math::
    C_4H_6(g) + C_2H_4(g) \rightarrow C_6H_{10}(g)

  .. image:: img/diels_alder_pes.png
    :align: center
    :width: 61.8%

  In the transition state, bonds that are breaking (three π bonds) are shown in orange and bonds that are forming (one π bond and two σ bonds) are shown in cyan.

  **(2) build a molecule template of the reaction site before the reaction has occurred**

  The pre-reacted molecule template is specified by a molecule command. This molecule template file contains a sample reaction site and its surrounding topology. All atom types in the pre-reacted template must be the same as those of a potential reaction site in the simulation. The initiator atom pairs of the pre-reacted template are specified by atom ID in the map file.

  .. code-block::
    :caption: `examples/PACKAGES/reaction/cyclohexene/1-3-butadiene.txt` (PubChem CID 7845)

    10 atoms
    9 bonds

    Types

      1 C
      2 C
      3 C
      4 C
      5 H
      6 H
      7 H
      8 H
      9 H
      10 H

    Coords

      1   -0.6022    0.3972    0.0000
      2    0.6024   -0.3975    0.0000
      3   -1.8315   -0.1305    0.0000
      4    1.8314    0.1308    0.0000
      5   -0.4975    1.4789    0.0001
      6    0.4979   -1.4792    0.0001
      7   -2.7035    0.5151    0.0000
      8   -1.9975   -1.2027    0.0000
      9    2.7036   -0.5143    0.0000
      10   1.9969    1.2030    0.0000

    Bonds

      1  1  2  1  0  0  0  0
      2  1  3  2  0  0  0  0
      3  1  5  1  0  0  0  0
      4  2  4  2  0  0  0  0
      5  2  6  1  0  0  0  0
      6  3  7  1  0  0  0  0
      7  3  8  1  0  0  0  0
      8  4  9  1  0  0  0  0
      9  4 10  1  0  0  0  0

    Charges

      1 -0.15
      10 0.15
      2 -0.15
      3 -0.3
      4 -0.3
      5 0.15
      6 0.15
      7 0.15
      8 0.15
      9 0.15


  .. code-block::
    :caption: `examples/PACKAGES/reaction/cyclohexene/ethene.txt` (PubChem CID 6325)

    6 atoms
    5 bonds

    Types

      1 C
      2 C
      3 H
      4 H
      5 H
      6 H

    Coords

      1 -0.6672  0.0000  0.0000
      2  0.6672  0.0000  0.0000
      3 -1.2213 -0.9290  0.0708
      4 -1.2212  0.9290 -0.0708
      5  1.2213  0.9290 -0.0708
      6  1.2213 -0.9290  0.0708

    Charges

      1 -0.3
      2 -0.3
      3 0.15
      4 0.15
      5 0.15
      6 0.15

    Bonds

      1  2  1  2
      2  1  1  3
      3  1  1  4
      4  1  2  5
      5  1  2  6

  **(3) build a molecule template of the reaction site after the reaction has occurred**

  The post-reacted molecule template contains a sample of the reaction site and its surrounding topology after the reaction has occurred. It must contain the same number of atoms as the pre-reacted template, unless there are created or deleted atoms (see examples/PACKAGES/reaction for details). A one-to-one correspondence between the atom IDs in the pre- and post-reacted templates is specified in the map file described below (4).

  .. code-block::
    :caption: `examples/PACKAGES/reaction/cyclohexene/cyclohexene.txt` (PubChem CID 8079)

    16 atoms
    16 bonds
    18.0153 mass # Molecular Weight 82.14 g/mol

    Types

      1 C
      2 C
      3 C
      4 C
      5 C
      6 C
      7 H
      8 H
      9 H
      10 H
      11 H
      12 H
      13 H
      14 H
      15 H
      16 H

    Coords

      1    0.6964   -1.2528   -0.3007
      2   -0.7059   -1.2475    0.3010
      3    1.4853   -0.0085    0.1150
      4   -1.4851    0.0027   -0.1155
      5    0.6721    1.2505    0.0597
      6   -0.6628    1.2555   -0.0595
      7    0.6205   -1.2901   -1.3951
      8    1.2373   -2.1536    0.0107
      9   -0.6303   -1.2847    1.3954
      10  -1.2535   -2.1442   -0.0101
      11   2.3597    0.0936   -0.5376
      12   1.8605   -0.1328    1.1378
      13  -2.3593    0.1112    0.5364
      14  -1.8604   -0.1189   -1.1386
      15   1.2038    2.1959    0.1161
      16  -1.1873    2.2049   -0.1154

    Charges

      15  0.15
      16  0.15
       3  0.14
       4  0.14
       5 -0.29
       6 -0.29

    Bonds

      1  1  2  1  0  0  0  0
      2  1  3  1  0  0  0  0
      3  1  7  1  0  0  0  0
      4  1  8  1  0  0  0  0
      5  2  4  1  0  0  0  0
      6  2  9  1  0  0  0  0
      7  2 10  1  0  0  0  0
      8  3  5  1  0  0  0  0
      9  3 11  1  0  0  0  0
      10 3 12  1  0  0  0  0
      11 4  6  1  0  0  0  0
      12 4 13  1  0  0  0  0
      13 4 14  1  0  0  0  0
      14 5  6  2  0  0  0  0
      15 5 15  1  0  0  0  0
      16 6 16  1  0  0  0  0


  **(4) create a map that relates the template-atom-IDs of each atom between pre- and post-reaction molecule templates**

  The header of map file contains one mandatory keyword *equivalences*\, which is the number of atoms in the pre- and post-reaction  molecule templates.

  The body of the map file contains two mandatory sections. The first mandatory section begins with the keyword *InitiatorIDs*\  listing the two atom IDs of the initiator atom pair in the pre-reacted molecule template. The second mandatory section begins with the keyword *Equivalences*\  listing a one-to-one correspondence between atom IDs of the pre- and post-reacted templates. The first column is an atom ID of the pre-reacted molecule template, and the second column is the corresponding atom ID of the post-reacted molecule template.

  .. code-block::
    :caption: `examples/PACKAGES/reaction/cyclohexene/cyclohexene_map.txt`

  **(5) fill a simulation box with molecules and run a simulation with fix bond/react**

  .. code-block::
    :caption: `examples/PACKAGES/reaction/cyclohexene/in.cyclohexene`







Termolecular reaction
"""""""""""""""""""""

An elementary *termolecular reaction* requires the collision of three particles at the same place and time. This type of reaction is extremely rare because all three reactants must simultaneously collide with each other, with sufficient energy and correct orientation, to produce a reaction. When a reaction involves three reactant molecules, it is much more likely for it to proceed via multiple steps known as a *reaction mechanism* involving elementary unimolecular and/or bimolecular reaction steps.

----------

Reaction mechanism
^^^^^^^^^^^^^^^^^^

A valid multi-step *reaction mechanism* consist of a series of unimolecular and/or bimolecular elementary reaction steps. The sum of the reaction steps should agree with the overall balanced reaction equation. A reaction intermediate is transient species within a multi-step reaction mechanism that is produced in the preceding step and consumed in a subsequent step to ultimately generate the final reaction product. Intermediate reactions are common in the biological world; a prime example can be seen in the metabolism of metabolites and nutrients.



Combustion reaction
"""""""""""""""""""

A *combustion reaction*, which also qualifies as a combination reaction, is a reaction in which a substance reacts with oxygen gas, releasing energy in the form of light and heat. Combustion reactions must involve O\ :sub:`2`\  as one reactant. The combustion of hydrogen gas produces water vapor:

.. math:: 2 H_2 + O_2 \rightarrow 2 H_2O

However this combustion reaction is actually more complicated than 2 diatomic hydrogen molecules simply colliding with a diatomic oxygen molecule to form 2 water molecules. It actually proceeds as a series of steps in the reaction mechanism:

.. math::
  H_2 \rightarrow H + H \\
  O_2 \rightarrow O + O \\
  H + O_2 \rightarrow O + OH \\
  H_2 + O \rightarrow H + OH \\
  H_2 + OH \rightarrow H_2O + H \\

..  youtube:: YuqA_uojSJ4
  :align: center
  :width: 90%

Oxyhydrogen is a mixture of hydrogen (H2) and oxygen (O2) gases, also known as "Knallgas" (from German lit. "bang-gas"). Theoretically, a ratio of 2:1 hydrogen:oxygen is enough to achieve maximum efficiency; in practice a ratio 4:1 or 5:1 is needed to avoid an oxidizing flame. Oxyhydrogen will combust when brought to its autoignition temperature of 843 K at normal atmospheric pressure. The minimum energy required to ignite such a mixture, at lower temperatures, with a spark is about 20 microjoules. At standard temperature and pressure, oxyhydrogen can burn when it is between about 4% and 95% hydrogen by volume.

When ignited, the gas mixture converts to water vapor and releases energy, which sustains the reaction: 241.8 kJ of energy (LHV) for every mole of H2 burned. The amount of heat energy released is independent of the mode of combustion, but the temperature of the flame varies. The maximum temperature of about 3073 K is achieved with an exact stoichiometric 2:1 mixture.

.. admonition:: REAXFF example: Hydrogen combustion
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

.. admonition:: REACTION example: Hydrogen combustion
  :class: Hint

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

.. admonition:: Example: Glucose-6-phosphate isomerase (GPI)
  :class: Hint

  *Glucose-6-phosphate isomerase (GPI)* is an enzyme (EC 5.3.1.9) that converts *glucose-6-phosphate (G6P)* to *fructose-6-phosphate (F6P)* as part of the glycolysis pathway. Since the reaction is reversible, its direction is determined by G6P and F6P concentrations. The mechanism that GPI uses to interconvert glucose 6-phosphate and fructose 6-phosphate consists of three major steps: opening the glucose ring, isomerizing glucose into fructose through an enediol intermediate, and closing the fructose ring. Functional GPI is a 64-kDa dimer composed of two identical monomers.[6][7] The two monomers interact notably through the two protrusions in a hugging embrace. The active site of each monomer is formed by a cleft between the two domains and the dimer interface. Human GPI pdb 1JLH. https://www.rcsb.org/structure/1jlh



Polymerization reaction
"""""""""""""""""""""""

.. admonition:: REACTION example: polymerization of nylon 6,6 :ref:`(Gissinger, 2020) <howto-reaction-Gissinger-2020>`
  :class: Hint

  .. code-block:: LAMMPS
    :caption: `examples/PACKAGES/reaction/nylon\,6-6_melt/in.large_nylon_melt`

    # 35,000 atom nylon melt example
    units real
    boundary p p p
    atom_style full
    kspace_style pppm 1.0e-4
    pair_style lj/class2/coul/long 8.5
    angle_style class2
    bond_style class2
    dihedral_style class2
    improper_style class2
    special_bonds lj/coul 0 0 1
    pair_modify tail yes mix sixthpower

    read_data large_nylon_melt.data.gz &
      extra/bond/per/atom 5  &
      extra/angle/per/atom 15 &
      extra/dihedral/per/atom 15 &
      extra/improper/per/atom 25 &
      extra/special/per/atom 25

    velocity all create 800.0 4928459 dist gaussian

    molecule mol1 rxn1_stp1_unreacted.molecule_template
    molecule mol2 rxn1_stp1_reacted.molecule_template
    molecule mol3 rxn1_stp2_unreacted.molecule_template
    molecule mol4 rxn1_stp2_reacted.molecule_template

    fix myrxns all bond/react stabilization yes statted_grp .03 &
      react rxn1 all 1 0.0 2.9 mol1 mol2 rxn1_stp1_map &
      react rxn2 all 1 0.0 5.0 mol3 mol4 rxn1_stp2_map

    # stable at 800K
    fix 1 statted_grp_REACT nvt temp 800 800 100

    thermo 50
    thermo_style custom step temp press density f_myrxns[*] # cumulative reaction counts
    run 200


  ..  youtube:: -dlyapmW7uI
    :align: center
    :width: 99%





























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
