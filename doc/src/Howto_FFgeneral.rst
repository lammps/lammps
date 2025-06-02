Some general force field considerations
=======================================

A compact summary of the concepts, definitions, and properties of force
fields with explicit bonded interactions (like the ones discussed in
this HowTo) is given in :ref:`(Gissinger) <Typelabel2>`.

A force field has 2 parts: the formulas that define its potential
functions and the coefficients used for a particular system.  To assign
parameters it is first required to assign atom types.  Those are not
only based on the elements, but also on the chemical environment due to
the atoms bound to them.  This often follows the chemical concept of
*functional groups*.  Example: a carbon atom bound with a single bond to
a single OH-group (alcohol) would be a different atom type than a carbon
atom bound to a methyl CH3 group (aliphatic carbon).  The atom types
usually then determine the non-bonded Lennard-Jones parameters and the
parameters for bonds, angles, dihedrals, and impropers.  On top of that,
partial charges have to be applied.  Those are usually independent of
the atom types and are determined either for groups of atoms called
residues with some fitting procedure based on quantum mechanical
calculations, or based on some increment system that add or subtract
increments from the partial charge of an atom based on the types of
the neighboring atoms.

Force fields differ in the strategies they employ to determine the
parameters and charge distribution in how generic or specific they are
which in turn has an impact on the accuracy (compare for example
CGenFF to CHARMM and GAFF to Amber).  Because of the different
strategies, it is not a good idea to use a mix of parameters from
different force field *families* (like CHARMM, Amber, or GROMOS)
and that extends to the parameters for the solvent, especially
water.  The publication describing the parameterization of a force
field will describe which water model to use.  Changing the water
model usually leads to overall worse results (even if it may improve
on the water itself).

In addition, one has to consider that *families* of force fields like
CHARMM, Amber, OPLS, or GROMOS have evolved over time and thus provide
different *revisions* of the force field parameters.  These often
corresponds to changes in the functional form or the parameterization
strategies.  This may also result in changes required for simulation
settings like the preferred cutoff or how Coulomb interactions are
computed (cutoff, smoothed/shifted cutoff, or long-range with Ewald
summation or equivalent).  Unless explicitly stated in the publication
describing the force field, the Coulomb interaction cannot be chosen at
will but must match the revision of the force field.  That said,
liberties may be taken during the initial equilibration of a system to
speed up the process, but not for production simulations.

----------

.. _Typelabel2:

**(Gissinger)** J. R. Gissinger, I. Nikiforov, Y. Afshar, B. Waters, M. Choi, D. S. Karls, A. Stukowski, W. Im, H. Heinz, A. Kohlmeyer, and E. B. Tadmor, J Phys Chem B, 128, 3282-3297 (2024).

