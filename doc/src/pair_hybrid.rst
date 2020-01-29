.. index:: pair\_style hybrid

pair\_style hybrid command
==========================

pair\_style hybrid/kk command
=============================

pair\_style hybrid/overlay command
==================================

pair\_style hybrid/overlay/kk command
=====================================

Syntax
""""""


.. parsed-literal::

   pair_style hybrid style1 args style2 args ...
   pair_style hybrid/overlay style1 args style2 args ...

* style1,style2 = list of one or more pair styles and their arguments

Examples
""""""""


.. parsed-literal::

   pair_style hybrid lj/cut/coul/cut 10.0 eam lj/cut 5.0
   pair_coeff 1\*2 1\*2 eam niu3
   pair_coeff 3 3 lj/cut/coul/cut 1.0 1.0
   pair_coeff 1\*2 3 lj/cut 0.5 1.2

   pair_style hybrid/overlay lj/cut 2.5 coul/long 2.0
   pair_coeff \* \* lj/cut 1.0 1.0
   pair_coeff \* \* coul/long

Description
"""""""""""

The *hybrid* and *hybrid/overlay* styles enable the use of multiple
pair styles in one simulation.  With the *hybrid* style, exactly one
pair style is assigned to each pair of atom types.  With the
*hybrid/overlay* style, one or more pair styles can be assigned to
each pair of atom types.  The assignment of pair styles to type pairs
is made via the :doc:`pair_coeff <pair_coeff>` command.

Here are two examples of hybrid simulations.  The *hybrid* style could
be used for a simulation of a metal droplet on a LJ surface.  The
metal atoms interact with each other via an *eam* potential, the
surface atoms interact with each other via a *lj/cut* potential, and
the metal/surface interaction is also computed via a *lj/cut*
potential.  The *hybrid/overlay* style could be used as in the 2nd
example above, where multiple potentials are superposed in an additive
fashion to compute the interaction between atoms.  In this example,
using *lj/cut* and *coul/long* together gives the same result as if
the *lj/cut/coul/long* potential were used by itself.  In this case,
it would be more efficient to use the single combined potential, but
in general any combination of pair potentials can be used together in
to produce an interaction that is not encoded in any single pair\_style
file, e.g. adding Coulombic forces between granular particles.

All pair styles that will be used are listed as "sub-styles" following
the *hybrid* or *hybrid/overlay* keyword, in any order.  Each
sub-style's name is followed by its usual arguments, as illustrated in
the example above.  See the doc pages of individual pair styles for a
listing and explanation of the appropriate arguments.

Note that an individual pair style can be used multiple times as a
sub-style.  For efficiency this should only be done if your model
requires it.  E.g. if you have different regions of Si and C atoms and
wish to use a Tersoff potential for pure Si for one set of atoms, and
a Tersoff potential for pure C for the other set (presumably with some
3rd potential for Si-C interactions), then the sub-style *tersoff*
could be listed twice.  But if you just want to use a Lennard-Jones or
other pairwise potential for several different atom type pairs in your
model, then you should just list the sub-style once and use the
pair\_coeff command to assign parameters for the different type pairs.

.. note::

   There is one exception to this option to list an individual
   pair style multiple times: GPU-enabled pair styles in the GPU package.
   This is because the GPU package currently assumes that only one
   instance of a pair style is being used.

In the pair\_coeff commands, the name of a pair style must be added
after the I,J type specification, with the remaining coefficients
being those appropriate to that style.  If the pair style is used
multiple times in the pair\_style command, then an additional numeric
argument must also be specified which is a number from 1 to M where M
is the number of times the sub-style was listed in the pair style
command.  The extra number indicates which instance of the sub-style
these coefficients apply to.

For example, consider a simulation with 3 atom types: types 1 and 2
are Ni atoms, type 3 are LJ atoms with charges.  The following
commands would set up a hybrid simulation:


.. parsed-literal::

   pair_style hybrid eam/alloy lj/cut/coul/cut 10.0 lj/cut 8.0
   pair_coeff \* \* eam/alloy nialhjea Ni Ni NULL
   pair_coeff 3 3 lj/cut/coul/cut 1.0 1.0
   pair_coeff 1\*2 3 lj/cut 0.8 1.3

As an example of using the same pair style multiple times, consider a
simulation with 2 atom types.  Type 1 is Si, type 2 is C.  The
following commands would model the Si atoms with Tersoff, the C atoms
with Tersoff, and the cross-interactions with Lennard-Jones:


.. parsed-literal::

   pair_style hybrid lj/cut 2.5 tersoff tersoff
   pair_coeff \* \* tersoff 1 Si.tersoff Si NULL
   pair_coeff \* \* tersoff 2 C.tersoff NULL C
   pair_coeff 1 2 lj/cut 1.0 1.5

If pair coefficients are specified in the data file read via the
:doc:`read_data <read_data>` command, then the same rule applies.
E.g. "eam/alloy" or "lj/cut" must be added after the atom type, for
each line in the "Pair Coeffs" section, e.g.


.. parsed-literal::

   Pair Coeffs

   1 lj/cut/coul/cut 1.0 1.0
   ...

Note that the pair\_coeff command for some potentials such as
:doc:`pair_style eam/alloy <pair_eam>` includes a mapping specification
of elements to all atom types, which in the hybrid case, can include
atom types not assigned to the *eam/alloy* potential.  The NULL
keyword is used by many such potentials (eam/alloy, Tersoff, AIREBO,
etc), to denote an atom type that will be assigned to a different
sub-style.

For the *hybrid* style, each atom type pair I,J is assigned to exactly
one sub-style.  Just as with a simulation using a single pair style,
if you specify the same atom type pair in a second pair\_coeff command,
the previous assignment will be overwritten.

For the *hybrid/overlay* style, each atom type pair I,J can be
assigned to one or more sub-styles.  If you specify the same atom type
pair in a second pair\_coeff command with a new sub-style, then the
second sub-style is added to the list of potentials that will be
calculated for two interacting atoms of those types.  If you specify
the same atom type pair in a second pair\_coeff command with a
sub-style that has already been defined for that pair of atoms, then
the new pair coefficients simply override the previous ones, as in the
normal usage of the pair\_coeff command.  E.g. these two sets of
commands are the same:


.. parsed-literal::

   pair_style lj/cut 2.5
   pair_coeff \* \* 1.0 1.0
   pair_coeff 2 2 1.5 0.8

   pair_style hybrid/overlay lj/cut 2.5
   pair_coeff \* \* lj/cut 1.0 1.0
   pair_coeff 2 2 lj/cut 1.5 0.8

Coefficients must be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as described above, or in the
data file or restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, or by mixing as described
below.

For both the *hybrid* and *hybrid/overlay* styles, every atom type
pair I,J (where I <= J) must be assigned to at least one sub-style via
the :doc:`pair_coeff <pair_coeff>` command as in the examples above, or
in the data file read by the :doc:`read_data <read_data>`, or by mixing
as described below.

If you want there to be no interactions between a particular pair of
atom types, you have 3 choices.  You can assign the type pair to some
sub-style and use the :doc:`neigh_modify exclude type <neigh_modify>`
command.  You can assign it to some sub-style and set the coefficients
so that there is effectively no interaction (e.g. epsilon = 0.0 in a
LJ potential).  Or, for *hybrid* and *hybrid/overlay* simulations, you
can use this form of the pair\_coeff command in your input script:


.. parsed-literal::

   pair_coeff        2 3 none

or this form in the "Pair Coeffs" section of the data file:


.. parsed-literal::

   3 none

If an assignment to *none* is made in a simulation with the
*hybrid/overlay* pair style, it wipes out all previous assignments of
that atom type pair to sub-styles.

Note that you may need to use an :doc:`atom_style <atom_style>` hybrid
command in your input script, if atoms in the simulation will need
attributes from several atom styles, due to using multiple pair
potentials.


----------


Different force fields (e.g. CHARMM vs AMBER) may have different rules
for applying weightings that change the strength of pairwise
interactions between pairs of atoms that are also 1-2, 1-3, and 1-4
neighbors in the molecular bond topology, as normally set by the
:doc:`special_bonds <special_bonds>` command.  Different weights can be
assigned to different pair hybrid sub-styles via the :doc:`pair_modify special <pair_modify>` command. This allows multiple force fields
to be used in a model of a hybrid system, however, there is no consistent
approach to determine parameters automatically for the interactions
between the two force fields, this is only recommended when particles
described by the different force fields do not mix.

Here is an example for mixing CHARMM and AMBER: The global *amber*
setting sets the 1-4 interactions to non-zero scaling factors and
then overrides them with 0.0 only for CHARMM:


.. parsed-literal::

   special_bonds amber
   pair_hybrid lj/charmm/coul/long 8.0 10.0 lj/cut/coul/long 10.0
   pair_modify pair lj/charmm/coul/long special lj/coul 0.0 0.0 0.0

The this input achieves the same effect:


.. parsed-literal::

   special_bonds 0.0 0.0 0.1
   pair_hybrid lj/charmm/coul/long 8.0 10.0 lj/cut/coul/long 10.0
   pair_modify pair lj/cut/coul/long special lj 0.0 0.0 0.5
   pair_modify pair lj/cut/coul/long special coul 0.0 0.0 0.83333333
   pair_modify pair lj/charmm/coul/long special lj/coul 0.0 0.0 0.0

Here is an example for mixing Tersoff with OPLS/AA based on
a data file that defines bonds for all atoms where for the
Tersoff part of the system the force constants for the bonded
interactions have been set to 0. Note the global settings are
effectively *lj/coul 0.0 0.0 0.5* as required for OPLS/AA:


.. parsed-literal::

   special_bonds lj/coul 1e-20 1e-20 0.5
   pair_hybrid tersoff lj/cut/coul/long 12.0
   pair_modify pair tersoff special lj/coul 1.0 1.0 1.0

For use with the various :doc:`compute \*/tally <compute_tally>`
computes, the :doc:`pair_modify compute/tally <pair_modify>`
command can be used to selectively turn off processing of
the compute tally styles, for example, if those pair styles
(e.g. many-body styles) do not support this feature.

See the :doc:`pair_modify <pair_modify>` doc page for details on
the specific syntax, requirements and restrictions.


----------


The potential energy contribution to the overall system due to an
individual sub-style can be accessed and output via the :doc:`compute pair <compute_pair>` command.


----------


.. note::

   Several of the potentials defined via the pair\_style command in
   LAMMPS are really many-body potentials, such as Tersoff, AIREBO, MEAM,
   ReaxFF, etc.  The way to think about using these potentials in a
   hybrid setting is as follows.

A subset of atom types is assigned to the many-body potential with a
single :doc:`pair_coeff <pair_coeff>` command, using "\* \*" to include
all types and the NULL keywords described above to exclude specific
types not assigned to that potential.  If types 1,3,4 were assigned in
that way (but not type 2), this means that all many-body interactions
between all atoms of types 1,3,4 will be computed by that potential.
Pair\_style hybrid allows interactions between type pairs 2-2, 1-2,
2-3, 2-4 to be specified for computation by other pair styles.  You
could even add a second interaction for 1-1 to be computed by another
pair style, assuming pair\_style hybrid/overlay is used.

But you should not, as a general rule, attempt to exclude the
many-body interactions for some subset of the type pairs within the
set of 1,3,4 interactions, e.g. exclude 1-1 or 1-3 interactions.  That
is not conceptually well-defined for many-body interactions, since the
potential will typically calculate energies and foces for small groups
of atoms, e.g. 3 or 4 atoms, using the neighbor lists of the atoms to
find the additional atoms in the group.  It is typically non-physical
to think of excluding an interaction between a particular pair of
atoms when the potential computes 3-body or 4-body interactions.

However, you can still use the pair\_coeff none setting or the
:doc:`neigh_modify exclude <neigh_modify>` command to exclude certain
type pairs from the neighbor list that will be passed to a many-body
sub-style.  This will alter the calculations made by a many-body
potential, since it builds its list of 3-body, 4-body, etc
interactions from the pair list.  You will need to think carefully as
to whether it produces a physically meaningful result for your model.

For example, imagine you have two atom types in your model, type 1 for
atoms in one surface, and type 2 for atoms in the other, and you wish
to use a Tersoff potential to compute interactions within each
surface, but not between surfaces.  Then either of these two command
sequences would implement that model:


.. parsed-literal::

   pair_style hybrid tersoff
   pair_coeff \* \* tersoff SiC.tersoff C C
   pair_coeff 1 2 none

   pair_style tersoff
   pair_coeff \* \* SiC.tersoff C C
   neigh_modify exclude type 1 2

Either way, only neighbor lists with 1-1 or 2-2 interactions would be
passed to the Tersoff potential, which means it would compute no
3-body interactions containing both type 1 and 2 atoms.

Here is another example, using hybrid/overlay, to use 2 many-body
potentials together, in an overlapping manner.  Imagine you have CNT
(C atoms) on a Si surface.  You want to use Tersoff for Si/Si and Si/C
interactions, and AIREBO for C/C interactions.  Si atoms are type 1; C
atoms are type 2.  Something like this will work:


.. parsed-literal::

   pair_style hybrid/overlay tersoff airebo 3.0
   pair_coeff \* \* tersoff SiC.tersoff.custom Si C
   pair_coeff \* \* airebo CH.airebo NULL C

Note that to prevent the Tersoff potential from computing C/C
interactions, you would need to modify the SiC.tersoff file to turn
off C/C interaction, i.e. by setting the appropriate coefficients to
0.0.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.

Since the *hybrid* and *hybrid/overlay* styles delegate computation to
the individual sub-styles, the suffix versions of the *hybrid* and
*hybrid/overlay* styles are used to propagate the corresponding suffix
to all sub-styles, if those versions exist. Otherwise the
non-accelerated version will be used.

The individual accelerated sub-styles are part of the GPU, USER-OMP
and OPT packages, respectively.  They are only enabled if LAMMPS was
built with those packages.  See the :doc:`Build package <Build_package>`
doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

Any pair potential settings made via the
:doc:`pair_modify <pair_modify>` command are passed along to all
sub-styles of the hybrid potential.

For atom type pairs I,J and I != J, if the sub-style assigned to I,I
and J,J is the same, and if the sub-style allows for mixing, then the
coefficients for I,J can be mixed.  This means you do not have to
specify a pair\_coeff command for I,J since the I,J type pair will be
assigned automatically to the sub-style defined for both I,I and J,J
and its coefficients generated by the mixing rule used by that
sub-style.  For the *hybrid/overlay* style, there is an additional
requirement that both the I,I and J,J pairs are assigned to a single
sub-style.  See the "pair\_modify" command for details of mixing rules.
See the See the doc page for the sub-style to see if allows for
mixing.

The hybrid pair styles supports the :doc:`pair_modify <pair_modify>`
shift, table, and tail options for an I,J pair interaction, if the
associated sub-style supports it.

For the hybrid pair styles, the list of sub-styles and their
respective settings are written to :doc:`binary restart files <restart>`, so a :doc:`pair_style <pair_style>` command does
not need to specified in an input script that reads a restart file.
However, the coefficient information is not stored in the restart
file.  Thus, pair\_coeff commands need to be re-specified in the
restart input script.

These pair styles support the use of the *inner*\ , *middle*\ , and
*outer* keywords of the :doc:`run_style respa <run_style>` command, if
their sub-styles do.

Restrictions
""""""""""""


When using a long-range Coulombic solver (via the
:doc:`kspace_style <kspace_style>` command) with a hybrid pair\_style,
one or more sub-styles will be of the "long" variety,
e.g. *lj/cut/coul/long* or *buck/coul/long*\ .  You must insure that the
short-range Coulombic cutoff used by each of these long pair styles is
the same or else LAMMPS will generate an error.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
