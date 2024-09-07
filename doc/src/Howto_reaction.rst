
Chemical reactions
==================

The five basic types of chemical reactions are combination, decomposition, single-replacement, double-replacement, and combustion :ref:`(LibreTexts) <howto-reaction-LibreTexts>`.

  .. image:: img/chemical_reactions.svg
    :align: center
    :width: 62%

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

.. admonition:: REACTION example: Diels–Alder
  :class: Hint

  :doc:`fix bond/react <fix_bond_react>` allows for complex topology changes during a running MD simulation, when using classical force fields. Topology changes are defined in pre- and post-reaction molecule templates and can include creation and deletion of bonds, angles, dihedrals, impropers, atom types, bond types, angle types, dihedral types, improper types, and/or atomic charges. :ref:`(Gissinger, 2017) <howto-reaction-Gissinger-2017>` and :ref:`(Gissinger, 2020) <howto-reaction-Gissinger-2020>`.

  A suggested workflow is:

  **(1) identify a reaction to be simulated**

  Gas-phase reaction of 1,3-butadiene with ethene to form cyclohexene :ref:`(Moore, Zhou, Garand) <howto-reaction-moore-zhou-garand>`:

  .. math::
    C_4H_6(g) + C_2H_4(g) \rightarrow C_6H_{10}(g)

  In the transition state, bonds that are breaking (three π bonds) are shown in orange and bonds that are forming (one π bond and two σ bonds) are shown in cyan.

  .. image:: img/diels_alder_pes.png
    :align: center
    :width: 62%


  **(2) build a molecule template of the reaction site before the reaction has occurred**

  The pre-reacted molecule template is specified by a molecule command. This molecule template file contains a sample reaction site and its surrounding topology. All atom types in the pre-reacted template must be the same as those of a potential reaction site in the simulation. The initiator atom pairs of the pre-reacted template are specified by atom ID in the map file.

  .. literalinclude:: ../../examples/PACKAGES/reaction/cyclohexene/1,3-butadiene.txt
    :caption: **examples/PACKAGES/reaction/cyclohexene/1,3-butadiene.txt**
    :class: code

  .. literalinclude:: ../../examples/PACKAGES/reaction/cyclohexene/ethene.txt
    :caption: **examples/PACKAGES/reaction/cyclohexene/ethene.txt**
    :class: code

  **(3) build a molecule template of the reaction site after the reaction has occurred**

  The post-reacted molecule template contains a sample of the reaction site and its surrounding topology after the reaction has occurred. It must contain the same number of atoms as the pre-reacted template, unless there are created or deleted atoms (see examples/PACKAGES/reaction for details). A one-to-one correspondence between the atom IDs in the pre- and post-reacted templates is specified in the map file described below (4).

  .. literalinclude:: ../../examples/PACKAGES/reaction/cyclohexene/cyclohexene.txt
    :caption: **examples/PACKAGES/reaction/cyclohexene/cyclohexene.txt**
    :class: code

  **(4) create a map that relates the template-atom-IDs of each atom between pre- and post-reaction molecule templates**

  The header of map file contains one mandatory keyword *equivalences*\, which is the number of atoms in the pre- and post-reaction molecule templates.

  The body of the map file contains two mandatory sections. The first mandatory section begins with the keyword *InitiatorIDs*\  listing the two atom IDs of the initiator atom pair in the pre-reacted molecule template. The second mandatory section begins with the keyword *Equivalences*\  listing a one-to-one correspondence between atom IDs of the pre- and post-reacted templates. The first column is an atom ID of the pre-reacted molecule template, and the second column is the corresponding atom ID of the post-reacted molecule template.

  Small molecules (i.e., ones that have all their atoms contained within the reaction templates) never have edge atoms.

  .. literalinclude:: ../../examples/PACKAGES/reaction/cyclohexene/cyclohexene_map.txt
    :caption: **examples/PACKAGES/reaction/cyclohexene/cyclohexene_map.txt**
    :class: code

  **(5) fill a simulation box with molecules and run a simulation with fix bond/react**

  .. literalinclude:: ../../examples/PACKAGES/reaction/cyclohexene/in.cyclohexene
    :caption: **examples/PACKAGES/reaction/cyclohexene/in.cyclohexene**
    :class: code







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
  :width: 62%

Oxyhydrogen is a mixture of hydrogen (H2) and oxygen (O2) gases, also known as "Knallgas" (from German lit. "bang-gas"). Theoretically, a ratio of 2:1 hydrogen:oxygen is enough to achieve maximum efficiency; in practice a ratio 4:1 or 5:1 is needed to avoid an oxidizing flame. Oxyhydrogen will combust when brought to its autoignition temperature of 843 K at normal atmospheric pressure. The minimum energy required to ignite such a mixture, at lower temperatures, with a spark is about 20 microjoules. At standard temperature and pressure, oxyhydrogen can burn when it is between about 4% and 95% hydrogen by volume.

When ignited, the gas mixture converts to water vapor and releases energy, which sustains the reaction: 241.8 kJ of energy (LHV) for every mole of H2 burned. The amount of heat energy released is independent of the mode of combustion, but the temperature of the flame varies. The maximum temperature of about 3073 K is achieved with an exact stoichiometric 2:1 mixture.

.. admonition:: REAXFF example: Hydrogen combustion
  :class: Hint

  .. literalinclude:: ../../examples/reaxff/hydrogen_combustion/H2.txt
    :caption: **examples/reaxff/hydrogen_combustion/H2.txt**
    :class: code

  .. literalinclude:: ../../examples/reaxff/hydrogen_combustion/O2.txt
    :caption: **examples/reaxff/hydrogen_combustion/O2.txt**
    :class: code

  .. literalinclude:: ../../examples/reaxff/hydrogen_combustion/in.hydrogen_combustion
    :caption: **examples/reaxff/hydrogen_combustion/in.hydrogen_combustion**
    :class: code
    :language: LAMMPS


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

  ..  youtube:: -dlyapmW7uI
    :align: center
    :width: 62%

  .. literalinclude:: ../../examples/PACKAGES/reaction/nylon,6-6_melt/in.large_nylon_melt
    :caption: **examples/PACKAGES/reaction/nylon,6-6_melt/in.large_nylon_melt**
    :class: code
    :language: LAMMPS
    
  .. literalinclude:: ../../examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp1_unreacted.molecule_template
    :caption: **examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp1_unreacted.molecule_template**
    :class: code

  .. literalinclude:: ../../examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp1_reacted.molecule_template
    :caption: **examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp1_reacted.molecule_template**
    :class: code

  .. literalinclude:: ../../examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp2_unreacted.molecule_template
    :caption: **examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp2_unreacted.molecule_template**
    :class: code

  .. literalinclude:: ../../examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp2_reacted.molecule_template
    :caption: **examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp2_reacted.molecule_template**
    :class: code

  .. literalinclude:: ../../examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp1_map
    :caption: **examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp1_map**
    :class: code

  .. literalinclude:: ../../examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp2_map
    :caption: **examples/PACKAGES/reaction/nylon,6-6_melt/rxn1_stp2_map**
    :class: code
































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
