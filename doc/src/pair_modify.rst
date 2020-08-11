.. index:: pair_modify

pair_modify command
===================

Syntax
""""""

.. code-block:: LAMMPS

   pair_modify keyword values ...

* one or more keyword/value pairs may be listed
* keyword = *pair* or *shift* or *mix* or *table* or *table/disp* or *tabinner*
  or *tabinner/disp* or *tail* or *compute* or *nofdotr* or *special* or
  *compute/tally*

  .. parsed-literal::

       *pair* value = sub-style N
         sub-style = sub-style of :doc:`pair hybrid <pair_hybrid>`
         N = which instance of sub-style (1 to M), only specify if sub-style is used multiple times
       *mix* value = *geometric* or *arithmetic* or *sixthpower*
       *shift* value = *yes* or *no*
       *table* value = N
         2\^N = # of values in table
       *table/disp* value = N
         2\^N = # of values in table
       *tabinner* value = cutoff
         cutoff = inner cutoff at which to begin table (distance units)
       *tabinner/disp* value = cutoff
         cutoff = inner cutoff at which to begin table (distance units)
       *tail* value = *yes* or *no*
       *compute* value = *yes* or *no*
       *nofdotr* value = none
       *special* values = which wt1 wt2 wt3
          which = *lj/coul* or *lj* or *coul*
          w1,w2,w3 = 1-2, 1-3, 1-4 weights from 0.0 to 1.0 inclusive
       *compute/tally* value = *yes* or *no*

Examples
""""""""

.. code-block:: LAMMPS

   pair_modify shift yes mix geometric
   pair_modify tail yes
   pair_modify table 12
   pair_modify pair lj/cut compute no
   pair_modify pair tersoff compute/tally no
   pair_modify pair lj/cut/coul/long 1 special lj/coul 0.0 0.0 0.0
   pair_modify pair lj/cut/coul/long special lj 0.0 0.0 0.5 special coul 0.0 0.0 0.8333333

Description
"""""""""""

Modify the parameters of the currently defined pair style.  If the
pair style is :doc:`hybrid or hybrid/overlay <pair_hybrid>`, then the
specified parameters are by default modified for all the hybrid sub-styles.

.. note::

   The behavior for hybrid pair styles can be changed by using the *pair*
   keyword, which allows selection of a specific sub-style to apply all
   remaining keywords to.
   The *special* and *compute/tally* keywords can **only** be
   used in conjunction with the *pair* keyword.  See further details about
   these 3 keywords below.

The *mix* keyword affects pair coefficients for interactions between
atoms of type I and J, when I != J and the coefficients are not
explicitly set in the input script.  Note that coefficients for I = J
must be set explicitly, either in the input script via the
:doc:`pair_coeff <pair_coeff>` command or in the "Pair Coeffs" section of the
:doc:`data file <read_data>`.  For some pair styles it is not
necessary to specify coefficients when I != J, since a "mixing" rule
will create them from the I,I and J,J settings.  The pair_modify
*mix* value determines what formulas are used to compute the mixed
coefficients.  In each case, the cutoff distance is mixed the same way
as sigma.

Note that not all pair styles support mixing and some mix options
are not available for certain pair styles. Also, there are additional
restrictions when using :doc:`pair style hybrid or hybrid/overlay <pair_hybrid>`.
See the doc page for individual pair styles for those restrictions.  Note also that the
:doc:`pair_coeff <pair_coeff>` command also can be used to directly set
coefficients for a specific I != J pairing, in which case no mixing is
performed.

- mix *geometric*

  .. math::

     \epsilon_{ij} = & \sqrt{\epsilon_i  \epsilon_j} \\
     \sigma_{ij}   = & \sqrt{\sigma_i  \sigma_j}

- mix *arithmetic*

  .. math::

    \epsilon_{ij} = & \sqrt{\epsilon_i  \epsilon_j} \\
    \sigma_{ij}   = & \frac{1}{2} (\sigma_i + \sigma_j)

- mix *sixthpower*

  .. math::

    \epsilon_{ij} = & \frac{2 \sqrt{\epsilon_i \epsilon_j} \sigma_i^3 \sigma_j^3}{\sigma_i^6 + \sigma_j^6} \\
    \sigma_{ij} =   & \left(\frac{1}{2} (\sigma_i^6 + \sigma_j^6) \right)^{\frac{1}{6}}

The *shift* keyword determines whether a Lennard-Jones potential is
shifted at its cutoff to 0.0.  If so, this adds an energy term to each
pairwise interaction which will be included in the thermodynamic
output, but does not affect pair forces or atom trajectories.  See the
doc page for individual pair styles to see which ones support this
option.

The *table* and *table/disp* keywords apply to pair styles with a
long-range Coulombic term or long-range dispersion term respectively;
see the doc page for individual styles to see which potentials support
these options.  If N is non-zero, a table of length 2\^N is
pre-computed for forces and energies, which can shrink their
computational cost by up to a factor of 2.  The table is indexed via a
bit-mapping technique :ref:`(Wolff) <Wolff1>` and a linear
interpolation is performed between adjacent table values.  In our
experiments with different table styles (lookup, linear, spline), this
method typically gave the best performance in terms of speed and
accuracy.

The choice of table length is a tradeoff in accuracy versus speed.  A
larger N yields more accurate force computations, but requires more
memory which can slow down the computation due to cache misses.  A
reasonable value of N is between 8 and 16.  The default value of 12
(table of length 4096) gives approximately the same accuracy as the
no-table (N = 0) option.  For N = 0, forces and energies are computed
directly, using a polynomial fit for the needed erfc() function
evaluation, which is what earlier versions of LAMMPS did.  Values
greater than 16 typically slow down the simulation and will not
improve accuracy; values from 1 to 8 give unreliable results.

The *tabinner* and *tabinner/disp* keywords set an inner cutoff above
which the pairwise computation is done by table lookup (if tables are
invoked), for the corresponding Coulombic and dispersion tables
discussed with the *table* and *table/disp* keywords.  The smaller the
cutoff is set, the less accurate the table becomes (for a given number
of table values), which can require use of larger tables.  The default
cutoff value is sqrt(2.0) distance units which means nearly all
pairwise interactions are computed via table lookup for simulations
with "real" units, but some close pairs may be computed directly
(non-table) for simulations with "lj" units.

When the *tail* keyword is set to *yes*\ , certain pair styles will
add a long-range VanderWaals tail "correction" to the energy and
pressure.  These corrections are bookkeeping terms which do not affect
dynamics, unless a constant-pressure simulation is being performed.
See the doc page for individual styles to see which support this
option.  These corrections are included in the calculation and
printing of thermodynamic quantities (see the :doc:`thermo_style
<thermo_style>` command).  Their effect will also be included in
constant NPT or NPH simulations where the pressure influences the
simulation box dimensions (e.g. the :doc:`fix npt <fix_nh>` and
:doc:`fix nph <fix_nh>` commands).  The formulas used for the
long-range corrections come from equation 5 of :ref:`(Sun) <Sun>`.

.. note::

   The tail correction terms are computed at the beginning of each
   run, using the current atom counts of each atom type.  If atoms are
   deleted (or lost) or created during a simulation, e.g. via the
   :doc:`fix gcmc <fix_gcmc>` command, the correction factors are not
   re-computed.  If you expect the counts to change dramatically, you
   can break a run into a series of shorter runs so that the
   correction factors are re-computed more frequently.

Several additional assumptions are inherent in using tail corrections,
including the following:

* The simulated system is a 3d bulk homogeneous liquid. This option
  should not be used for systems that are non-liquid, 2d, have a slab
  geometry (only 2d periodic), or inhomogeneous.
* G(r), the radial distribution function (rdf), is unity beyond the
  cutoff, so a fairly large cutoff should be used (i.e. 2.5 sigma for
  an LJ fluid), and it is probably a good idea to verify this
  assumption by checking the rdf.  The rdf is not exactly unity beyond
  the cutoff for each pair of interaction types, so the tail
  correction is necessarily an approximation.

  The tail corrections are computed at the beginning of each
  simulation run.  If the number of atoms changes during the run,
  e.g. due to atoms leaving the simulation domain, or use of the
  :doc:`fix gcmc <fix_gcmc>` command, then the corrections are not
  updated to reflect the changed atom count.  If this is a large
  effect in your simulation, you should break the long run into
  several short runs, so that the correction factors are re-computed
  multiple times.

* Thermophysical properties obtained from calculations with this
  option enabled will not be thermodynamically consistent with the
  truncated force-field that was used.  In other words, atoms do not
  feel any LJ pair interactions beyond the cutoff, but the energy and
  pressure reported by the simulation include an estimated
  contribution from those interactions.

The *compute* keyword allows pairwise computations to be turned off,
even though a :doc:`pair_style <pair_style>` is defined.  This is not
useful for running a real simulation, but can be useful for debugging
purposes or for performing a :doc:`rerun <rerun>` simulation, when you
only wish to compute partial forces that do not include the pairwise
contribution.

Two examples are as follows.  First, this option allows you to perform
a simulation with :doc:`pair_style hybrid <pair_hybrid>` with only a
subset of the hybrid sub-styles enabled.  Second, this option allows
you to perform a simulation with only long-range interactions but no
short-range pairwise interactions.  Doing this by simply not defining
a pair style will not work, because the :doc:`kspace_style
<kspace_style>` command requires a Kspace-compatible pair style be
defined.

The *nofdotr* keyword allows to disable an optimization that computes
the global stress tensor from the total forces and atom positions
rather than from summing forces between individual pairs of atoms.

----------

The *pair* keyword can only be used with the :doc:`hybrid and
hybrid/overlay <pair_hybrid>` pair styles.  If used, it must appear
first in the list of keywords.

Its meaning is that all the following parameters will only be modified
for the specified sub-style.  If the sub-style is defined multiple
times, then an additional numeric argument *N* must also be specified,
which is a number from 1 to M where M is the number of times the
sub-style was listed in the :doc:`pair_style hybrid <pair_hybrid>`
command.  The extra number indicates which instance of the sub-style
the remaining keywords will be applied to.

The *special* and *compute/tally* keywords can **only** be used in
conjunction with the *pair* keyword and they must directly follow it.
I.e. any other keyword, must appear after *pair*, *special*, and
*compute/tally*.

The *special* keyword overrides the global :doc:`special_bonds <special_bonds>`
1-2, 1-3, 1-4 exclusion settings (weights) for the sub-style selected
by the *pair* keyword.

Similar to the :doc:`special_bonds <special_bonds>` command, it takes
4 arguments.  The *which* argument can be *lj* to change only the
non-Coulomb weights (e.g. Lennard-Jones or Buckingham), *coul* to change
only the Coulombic settings, or *lj/coul* to change both to the same
values.  The *wt1,wt2,wt3* values are numeric weights from 0.0 to 1.0
inclusive, for the 1-2, 1-3, and 1-4 bond topology neighbors, respectively.
The *special* keyword can be used multiple times, e.g. to set the *lj*
and *coul* settings to different values.

.. note::

   The *special* keyword is not compatible with pair styles from the
   GPU or the USER-INTEL package and attempting to use it will cause
   an error.

.. note::

   Weights of exactly 0.0 or 1.0 in the :doc:`special_bonds <special_bonds>`
   command have implications on the neighbor list construction, which
   means that they cannot be overridden by using the *special* keyword.
   One workaround for this restriction is to use the :doc:`special_bonds <special_bonds>`
   command with weights like 1.0e-10 or 0.999999999 instead of 0.0 or 1.0,
   respectively, which enables to reset each them to any value between 0.0
   and 1.0 inclusively.  Otherwise you can set **all** global weights to
   an arbitrary number between 0.0 or 1.0, like 0.5, and then you have
   to override **all** *special* settings for **all** sub-styles which use
   the 1-2, 1-3, and 1-4 exclusion weights in their force/energy computation.

The *compute/tally* keyword disables or enables registering :doc:`compute
\*/tally <compute_tally>` computes for the sub-style specified by
the *pair* keyword.  Use *no* to disable, or *yes* to enable.

.. note::

   The "pair_modify pair compute/tally" command must be issued
   **before** the corresponding compute style is defined.

----------

Restrictions
""""""""""""

You cannot use *shift* yes with *tail* yes, since those are
conflicting options.  You cannot use *tail* yes with 2d simulations.
You cannot use *special* with pair styles from the GPU or
USER-INTEL package.

Related commands
""""""""""""""""

:doc:`pair_style <pair_style>`, :doc:`pair_style hybrid <pair_hybrid>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`thermo_style <thermo_style>`,
:doc:`compute \*/tally <compute_tally>`

Default
"""""""

The option defaults are mix = geometric, shift = no, table = 12,
tabinner = sqrt(2.0), tail = no, and compute = yes.

Note that some pair styles perform mixing, but only a certain style of
mixing.  See the doc pages for individual pair styles for details.

----------

.. _Wolff1:

**(Wolff)** Wolff and Rudd, Comp Phys Comm, 120, 200-32 (1999).

.. _Sun:

**(Sun)** Sun, J Phys Chem B, 102, 7338-7364 (1998).
