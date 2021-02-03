
.. Yuan documentation master file, created by
   sphinx-quickstart on Sat Jan 30 14:06:22 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   tc387: Multiple text additions/changes, Feb 2 2021
.. index:: fix gcmc

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
This fix performs Monte Carlo (MC) sampling of charge regulation and exchange of ions with a reservoir [1,2].  Monte Carlo moves alter the charge value of specific atoms and simultaneously insert a counterion into the system to preserve charge neutrality. Besides, the ionic strength of the simulated solution is controlled by the grand canonical MC insertion/deletion of ion pairs with a reservoir. The implemented method is largely analogous to the grand-reaction ensemble method [3].

The implementation is parallelized, straightforward to use, compatible with existing LAMMPS functionalities, and applicable to any system utilizing discreet, ionizable groups or surface sites.
The fix requires a LAMMPS atom style with a “charge” attribute, for example *charge* or *full*. 
The scheme is limited to integer charges, and therefore any atoms with non-integer charges will not be considered. 
All particles participating in the dissociation/association process are limited to monovalent charges, :math:`q = [+1, 0, −1]`, while salt ions of multivalency and any integer charge values is allowed. In our implementation "acid" refers to particles which can attain charge :math:`q=[0,-1]` and "base" to particles with :math:`q=[0,1]`. 

MC moves alter the charge value of acid or base particles, which models the dissociation/association process.
A standard acid ionization reaction (:math:`\mathrm{A} \rightleftharpoons \mathrm{A}^-+\mathrm{H}^+`) can be defined via a single line in the input file

.. code-block:: LAMMPS

    fix acid_reaction all charge_regulation 2 3 acid_type 1 pH 7.0 pKa 5.0 pIp 7.0 pIm 7.0
where the fix will attempt to charge :math:`\mathrm{A}` (discharge :math:`\mathrm{A}^-`) to :math:`\mathrm{A}^-` (:math:`\mathrm{A}`) and insert (delete) a monovalent proton of atom type 2. 
In addition, the fix also inserts and deletes ion pair of a proton :math:`\mathrm{H}^+` (atom type 2) and a hydroxyl ion :math:`\mathrm{OH}^-` (atom type 3).
This is adopted in the reaction ensemble approach [4] which is known to be highly inefficient 
at :math:`\mathrm{pH} \approx 7` when the concentration of both protons and hydroxyl ions is low, resulting in a relatively low acceptance rate of MC moves.
An efficient and correct solution to overcome this issue is to allow salt ions to 
participate in ionization reactions, which can easily achieved via for example, 

.. code-block:: LAMMPS

    fix acid_reaction all charge_regulation 2 3 acid_type 1 pH 7.0 pKa 5.0 pIp 2.0 pIm 2.0
where particles of atom type 2 contain both protons and free salt cations and particles of atom type 3 contain both hydroxyl ions and free salt anions (note the change in :math:`\mathrm{pIp}` and :math:`\mathrm{pIm}`) . See Refs. [1,3] for more details.

Similarly, a standard base ionization reaction (:math:`\mathrm{B} \rightleftharpoons \mathrm{B}^++\mathrm{OH}^-`) 
can be defined via 

.. code-block:: LAMMPS

    fix base_reaction all charge_regulation 2 3 base_type 4 pH 7.0 pKb 6.0 pIp 7.0 pIm 7.0
where the fix will attempt to charge :math:`\mathrm{B}` (discharge :math:`\mathrm{B}^+`) to :math:`\mathrm{B}^+` (:math:`\mathrm{B}`) and insert (delete) a monovalent anion (e.g. hydroxyl ion) of atom type 3, 
as well as to insert and delete an ion pair of a proton :math:`\mathrm{H}^+` (atom type 2) and a hydroxyl ion :math:`\mathrm{OH}^-` (atom type 3).

If neither the acid_type nor the base_type is specified, for example 

.. code-block:: LAMMPS

    fix salt_reaction all charge_regulation 2 3 pH 7.0 pIp 7.0 pIm 7.0
the fix simply inserts or deletes an ion pair of a proton :math:`\mathrm{H}^+` (atom type 2) and a hydroxyl ion :math:`\mathrm{OH}^-` (atom type 3) at the chemical potential :math:`\mathrm{pIp}=\mathrm{pIm}=7.0`, which corresponds to the self-ionization of water (:math:`\emptyset \rightleftharpoons \mathrm{H}^++\mathrm{OH}^-`) in this case.

The fix is compatible with LAMMPS sub-packages such as *molecule* or *rigid*. That said, the acid and base particles can be part of larger molecules or rigid bodies. Free ions that are inserted to or deleted from the system must be deﬁned as single particles (no bonded interactions allowed) and cannot be part of larger molecules or rigid bodies. If *molecule* package is used, all inserted ions have a molecule ID equal to zero.

Note that LAMMPS implicitly assumes a constant number of particles (degrees of freedom). Since using this fix alters the total number of particles changes during the simulation, any thermostat used by LAMMPS, such as NVT or Langevin, must use a dynamic calculation of system temperature. This can be achieved by specifying a dynamic temperature compute (e.g. dtemp) and using it with the desired thermostat, e.g. for a Langevin thermostat:

compute dtemp all temp

compute_modify dtemp dynamic yes 

fix fT all langevin 1.0 1.0 1.0 123 

fix_modify fT temp dtemp


The chemical potential units (e.g. pH) are in the standard log10 representation assuming reference concentration :math:`\rho_0 = {mol}/{l}`. 
Therefore, to perform internal unit conversion, the length (in nanometers) of the LAMMPS unit length 
must be specified via *lunit_nm* (default is set to the Bjerrum length in water at room temprature *lunit_nm* = 0.72nm). For example, in the dilute ideal solution limit, the concentration of free ions 
will be :math:`c_I = 10^{-pIp}{mol}/{l}`.

The temperature used in MC acceptance probability is set by  *temp*. This temperature should be the same as the temperature set by the molecular dynamics thermostat. For most purposes, it is probably the best to use *tempfixid* keyword which sets the temperature equal to the chosen MD thermostat temperature, in the example above we assumed the thermostat fix-ID is *fT*. The inserted particles attain a random velocity corresponding to the specified temperature. Suing *tempfixid* overrides any fixed temperature set by *temp*.   

The *xrd* keyword is can be used to restrict the inserted/deleted counterions to a specific radial distance from the chosen acid or base particle. This can be used to simulate more realist reaction dynamics. If *xrd* = 0 or *xrd* > *L* / 2, where *L* is the smallest box dimension, the radial restriction is automatically turned off and particles can be inserted or deleted anywhere in the box. 

If the *tag yes* is used, every inserted atom gets a unique tag ID, otherwise, the tag of every inserted atom is set to 0. *tag yes* might cause an integer overflow in very long simulations since the tags are unique to every particle and thus increase with every successful particle insertion. 

The fix only attempts to perform particle charging MC moves if *acid_type* or *base_type* is defined. Otherwise fix only performs free ion insertion/deletion. For example, if *acid_type* is not defined, *pmcA* is automatically set to 0. The vector *pmcmoves* is automatically normalized, for example, if set to *pmcmoves* 0 0.33 0.33, the vector would be normalized to [0,0.5,0.5]. 

The *only_salt* option can be used to perform multivalent grand-canonical ion-exchange moves. If *only_salt yes* is used, no charge exchange is performed, only ion insertion/deletion (*pmcmoves* is set to [0,0,1]), but ions can be multivalent. In the example above, an MC move would consist of three ion insertion/deletion to preserve charge neutrality of the system.

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
This fix is part of the USER-MISC package. It is only enabled if LAMMPS was built with that package. See the Build package doc page for more info.

The *atom_style* used must contain the *charge* property, for example, the style could be *charge* or *full* style. Only usable for 3D simulations. Atoms specified as free ions cannot be part of rigid bodies or molecules and cannot have bonding interactions.

Note: Regions restrictions are not yet implemented.

Default
"""""""
pH = 7.0; pKa = 100.0; pKb = 100.0; pIp = 100; pIm = 100; pKs=14.0; acid_type = -1; base_type = -1; lunit_nm = 0.72; temp = 1.0; nevery = 100; nmc = 100; xrd = 0; seed = 2345; tag = no; onlysalt = no, pmcmoves = 0.33 0.33 0.33, group-ID = all

References
""""""""""

[1] T. Curk, J. Yuan, and E. Luijten, "Coarse-grained simulation of charge regulation using LAMMPS", preprint (2021).

[2] T. Curk and E. Luijten, "Charge-regulation effects in nanoparticle self-assembly", PRL (2021)

[3] J. Landsgesell, P. Hebbeker, O. Rud, R. Lunkad, P. Kosovan, and C. Holm, “Grand-reaction method for simulations of ionization equilibria coupled to ion partitioning,” Macromolecules 53, 3007–3020 (2020).

[4] W. R. Smith and B. Triska, “The reaction ensemble method for the computer simulation of chemical and phase equilibria. I. Theory and basic examples,” J. Chem. Phys. 100, 3019–3027 (1994).
