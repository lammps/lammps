
.. Yuan documentation master file, created by
   sphinx-quickstart on Sat Jan 30 14:06:22 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   tc387: Multiple text additions/changes, Feb 2 2021
.. index:: fix fix_charge_regulation

fix_charge_regulation command
=============================
Syntax
""""""

.. parsed-literal::
   
    fix ID group-ID charge_regulation cation_type anion_type keyword value(s)

* ID, group-ID are documented in fix command
* charge_regulation = style name of this fix command
* cation_type = atom type of free cations
* anion_type = atom type of free anions
  
* zero or more keyword/value pairs may be appended

  .. parsed-literal::
     
     keyword = *pH*, *pKa*, *pKb*, *pIp*, *pIm*, *pKs*, *acid_type*, *base_type*, *lunit_nm*, *temp*, *tempfixid*, *nevery*, *nmc*, *xrd*, *seed*, *tag*, *group*, *onlysalt*, *pmcmoves* 
     *pH* value = pH of the solution
     *pKa* value = acid dissociation constant 
     *pKb* value = base dissociation constant
     *pIp* value = chemical potential of free cations
     *pIm* value = chemical potential of free anions
     *pKs* value = solution self-dissociation constant
     *acid_type* = atom type of acid groups
     *base_type*  = atom type of base groups
     *lunit_nm* value = unit length used by LAMMPS (# in the units of nanometers)
     *temp* value = temperature 
     *tempfixid* value = fix ID of temperature thermostat
     *nevery* value = invoke this fix every nevery steps
     *nmc* value = number of charge regulation MC moves to attempt every nevery steps
     *xrd* value = cutoff distance for acid/base reaction
     *seed* value = random # seed (positive integer)
     *tag* value = yes or no (yes: The code assign unique tags to inserted ions; no: The tag of all inserted ions is "0")
     *group* value = group-ID, inserted ions are assigned to group group-ID. Can be used multiple times to assign inserted ions to multiple groups.
     *onlysalt* values = flag charge_cation charge_anion. 
        flag = yes or no (yes: the fix performs only ion insertion/deletion, no: perform acid/base dissociation and ion insertion/deletion)
        charge_cation, charge_anion = value of cation/anion charge, must be an integer (only specify if flag = yes)
     *pmcmoves* values = pmcA pmcB pmcI -  MC move fractions for acid ionization (pmcA), base ionization (pmcB) and free ion exchange (pmcI) 

Examples
""""""""
.. code-block:: LAMMPS

    fix chareg all charge_regulation 1 2 acid_type 3 base_type 4 pKa 5 pKb 7 lb 1.0 nevery 200 nexchange 200 seed 123 tempfixid fT 

    fix chareg all charge_regulation 1 2 pIp 3 pIm 3 tempfixid fT tag yes onlysalt yes 2 -1

Description
"""""""""""
This fix performs Monte Carlo (MC) sampling of charge regulation and exchange of ions with a reservoir as discussed in :ref:`(Curk1) <Curk1>` and :ref:`(Curk2) <Curk2>`.  
The implemented method is largely analogous to the grand-reaction ensemble method in :ref:`(Landsgesell) <Landsgesell>`.
The implementation is parallelized, straightforward to use, compatible with existing LAMMPS functionalities, and applicable to any system utilizing discreet, ionizable groups or surface sites.
The fix requires a LAMMPS atom style with a “charge” attribute, for example *charge* or *full*. 

Specifically, the following three types of general reactions are implemented, including :math:`\mathrm{A} \rightleftharpoons \mathrm{A}^-+\mathrm{X}^+`, :math:`\mathrm{B} \rightleftharpoons \mathrm{B}^++\mathrm{X}^-`,
and :math:`\emptyset \rightleftharpoons Z^-\mathrm{X}^{Z^+}+Z^+\mathrm{X}^{-Z^-}`
where the particles include acid groups (neutral acid molecule :math:`\mathrm{A}` and negatively charged ionization state :math:`\mathrm{A}^-`), base groups (neutral base molecule :math:`\mathrm{B}` and positively charged ionization state :math:`\mathrm{B}^+`), free cations (:math:`\mathrm{X}^{Z^+}` with valency :math:`{Z^+}`), and free anions (:math:`\mathrm{X}^{-Z^-}` with valency :math:`-{Z^-}`).
In the former two types of reactions, Monte Carlo moves alter the charge value of specific atoms and simultaneously insert a counterion into the system to preserve charge neutrality, which models the dissociation/association process.
The last type of reaction allows the grand canonical MC exchange of ion pairs with a reservoir, which can be used to control the ionic strength of the solution.
The scheme is limited to integer charges only, 
and any atoms with non-integer charges will not be considered. 
In our implementation "acid" refers to particles that can attain charge :math:`q=[0,-1]` and "base" to particles with :math:`q=[0,1]`,
whereas the MC exchange of salt ions allows any integer charge values of :math:`{Z^+}` and :math:`{Z^-}`.
Multiple reactions can be added to the same simulation system.

Below we give several practical examples for modeling charge regulation eﬀects in solvated systems.
An acid ionization reaction (:math:`\mathrm{A} \rightleftharpoons \mathrm{A}^-+\mathrm{H}^+`) can be defined via a single line in the input file

.. code-block:: LAMMPS

    fix acid_reaction all charge_regulation 2 3 acid_type 1 pH 7.0 pKa 5.0 pIp 7.0 pIm 7.0
where the fix attempts to charge :math:`\mathrm{A}` (discharge :math:`\mathrm{A}^-`) to :math:`\mathrm{A}^-` (:math:`\mathrm{A}`) and insert (delete) a proton of atom type 2. 
Besides, the fix implements self-ionization reaction of water :math:`\emptyset \rightleftharpoons \mathrm{H}^++\mathrm{OH}^-`.
However, this approach is highly inefficient at :math:`\mathrm{pH} \approx 7` when the concentration of both protons and hydroxyl ions is low, resulting in a relatively low acceptance rate of MC moves.
A far more efficient and correct solution to this issue is to allow salt ions to 
participate in ionization reactions, which can be easily achieved via 

.. code-block:: LAMMPS

    fix acid_reaction all charge_regulation 2 3 acid_type 1 pH 7.0 pKa 5.0 pIp 2.0 pIm 2.0
where particles of atom type 2 contain both protons and free salt cations. See :ref:`(Curk1) <Curk1>` 
and :ref:`(Landsgesell) <Landsgesell>` for more details.

Similarly, a base ionization reaction (:math:`\mathrm{B} \rightleftharpoons \mathrm{B}^++\mathrm{OH}^-`) 
can be defined via 

.. code-block:: LAMMPS

    fix base_reaction all charge_regulation 2 3 base_type 4 pH 7.0 pKb 6.0 pIp 7.0 pIm 7.0
where the fix will attempt to charge :math:`\mathrm{B}` (discharge :math:`\mathrm{B}^+`) to :math:`\mathrm{B}^+` (:math:`\mathrm{B}`) and insert (delete) a hydroxyl ion  :math:`\mathrm{OH}^-` of atom type 3.
If in the command line neither the acid_type nor the base_type is specified, for example 

.. code-block:: LAMMPS

    fix salt_reaction all charge_regulation 2 3 pH 7.0 pIp 2.0 pIm 2.0
the fix simply inserts or deletes an ion pair of a free cation (atom type 2) and a free anion (atom type 3)
as is done in a conventional grand-canonical MC simulation.

The fix is compatible with LAMMPS sub-packages such as *molecule* or *rigid*. That said, the acid and base particles can be part of larger molecules or rigid bodies. Free ions that are inserted to or deleted from the system must be deﬁned as single particles (no bonded interactions allowed) and cannot be part of larger molecules or rigid bodies. If *molecule* package is used, all inserted ions have a molecule ID equal to zero.

Note that LAMMPS implicitly assumes a constant number of particles (degrees of freedom). Since using this fix alters the total number of particles during the simulation, any thermostat used by LAMMPS, such as NVT or Langevin, must use a dynamic calculation of system temperature. This can be achieved by specifying a dynamic temperature compute (e.g. dtemp) and using it with the desired thermostat, e.g. a Langevin thermostat:

.. code-block:: LAMMPS

    compute dtemp all temp
    compute_modify dtemp dynamic yes 
    fix fT all langevin 1.0 1.0 1.0 123 
    fix_modify fT temp dtemp

The chemical potential units (e.g. pH) are in the standard log10 representation assuming reference concentration :math:`\rho_0 = {mol}/{l}`. 
Therefore, to perform the internal unit conversion, the length (in nanometers) of the LAMMPS unit length 
must be specified via *lunit_nm* (default is set to the Bjerrum length in water at room temprature *lunit_nm* = 0.72nm). For example, in the dilute ideal solution limit, the concentration of free ions 
will be :math:`c_\mathrm{I} = 10^{-\mathrm{pIp}}{mol}/{l}`.

The temperature used in MC acceptance probability is set by  *temp*. This temperature should be the same as the temperature set by the molecular dynamics thermostat. For most purposes, it is probably the best to use *tempfixid* keyword which sets the temperature equal to the chosen MD thermostat temperature, in the example above we assumed the thermostat fix-ID is *fT*. The inserted particles attain a random velocity corresponding to the specified temperature. Suing *tempfixid* overrides any fixed temperature set by *temp*.   

The *xrd* keyword is can be used to restrict the inserted/deleted counterions to a specific radial distance from the chosen acid or base particle. This can be used to simulate more realist reaction dynamics. If *xrd* = 0 or *xrd* > *L* / 2, where *L* is the smallest box dimension, the radial restriction is automatically turned off and particles can be inserted or deleted anywhere in the box. 

If the *tag yes* is used, every inserted atom gets a unique tag ID, otherwise, the tag of every inserted atom is set to 0. *tag yes* might cause an integer overflow in very long simulations since the tags are unique to every particle and thus increase with every successful particle insertion. 

The fix only attempts to perform particle charging MC moves if *acid_type* or *base_type* is defined. Otherwise fix only performs free ion insertion/deletion. For example, if *acid_type* is not defined, *pmcA* is automatically set to 0. The vector *pmcmoves* is automatically normalized, for example, if set to *pmcmoves* 0 0.33 0.33, the vector would be normalized to [0,0.5,0.5]. 

The *only_salt* option can be used to perform multivalent grand-canonical ion-exchange moves. If *only_salt yes* is used, no charge exchange is performed, only ion insertion/deletion (*pmcmoves* is set to [0,0,1]), but ions can be multivalent. In the example above, an MC move would consist of three ion insertion/deletion to preserve the charge neutrality of the system.

The *group* keyword can be used to add inserted particles to a specific group-ID. All inserted particles are automatically added to group *all*.


Output
""""""
This fix computes a global vector of length 8, which can be accessed by various output commands. The vector values are the following global cumulative quantities:

* 1 = cumulative MC attempts
* 2 = cumulative MC successes
* 3 = current # of neutral acid atoms 
* 4 = current # of -1 charged acid atoms 
* 5 = current # of neutral base atoms 
* 6 = current # of +1 charged acid atoms 
* 7 = current # of free cations 
* 8 = current # of free anions


Restrictions
""""""""""""
This fix is part of the USER-MISC package. It is only enabled if LAMMPS was built with that package.
See the :doc:`Build package <Build_package>` doc page for more info.

The *atom_style* used must contain the *charge* property, for example, the style could be *charge* or *full* style. Only usable for 3D simulations. Atoms specified as free ions cannot be part of rigid bodies or molecules and cannot have bonding interactions.

Note: Regions restrictions are not yet implemented.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nh>`,
:doc:`fix langevin <fix_nh>`

Default
"""""""
pH = 7.0; pKa = 100.0; pKb = 100.0; pIp = 100; pIm = 100; pKs=14.0; acid_type = -1; base_type = -1; lunit_nm = 0.72; temp = 1.0; nevery = 100; nmc = 100; xrd = 0; seed = 2345; tag = no; onlysalt = no, pmcmoves = 0.33 0.33 0.33, group-ID = all

----------

.. _Curk1:

**(Curk1)** T. Curk, J. Yuan, and E. Luijten, "Coarse-grained simulation of charge regulation using LAMMPS", preprint (2021).

.. _Curk2:

**(Curk2)** T. Curk and E. Luijten, "Charge-regulation effects in nanoparticle self-assembly", PRL (2021)

.. _Landsgesell:

**(Landsgesell)** J. Landsgesell, P. Hebbeker, O. Rud, R. Lunkad, P. Kosovan, and C. Holm, “Grand-reaction method for simulations of ionization equilibria coupled to ion partitioning,” Macromolecules 53, 3007–3020 (2020).
