.. index:: fix adapt

fix adapt command
=================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID adapt N attribute args ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* adapt = style name of this fix command
* N = adapt simulation settings every this many timesteps
* one or more attribute/arg pairs may be appended
* attribute = *pair* or *bond* or *angle* or *kspace* or *atom*

  .. parsed-literal::

       *pair* args = pstyle pparam I J v_name
         pstyle = pair style name (e.g., lj/cut)
         pparam = parameter to adapt over time
         I,J = type pair(s) to set parameter for
         v_name = variable with name that calculates value of pparam
       *bond* args = bstyle bparam I v_name
         bstyle = bond style name (e.g., harmonic)
         bparam = parameter to adapt over time
         I = type bond to set parameter for
         v_name = variable with name that calculates value of bparam
       *angle* args = astyle aparam I v_name
         astyle = angle style name (e.g., harmonic)
         aparam = parameter to adapt over time
         I = type angle to set parameter for
         v_name = variable with name that calculates value of aparam
       *kspace* arg = v_name
         v_name = variable with name that calculates scale factor on :math:`k`-space terms
       *atom* args = atomparam v_name
         atomparam = *charge* or *diameter* or *diameter/disc* = parameter to adapt over time
         v_name = variable with name that calculates value of atomparam

* zero or more keyword/value pairs may be appended
* keyword = *scale* or *reset* or *mass*

  .. parsed-literal::

     *scale* value = *no* or *yes*
       *no* = the variable value is the new setting
       *yes* = the variable value multiplies the original setting
     *reset* value = *no* or *yes*
       *no* = values will remain altered at the end of a run
       *yes* = reset altered values to their original values at the end of a run
     *mass* value = *no* or *yes*
       *no* = mass is not altered by changes in diameter
       *yes* = mass is altered by changes in diameter

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all adapt 1 pair soft a 1 1 v_prefactor
   fix 1 all adapt 1 pair soft a 2* 3 v_prefactor
   fix 1 all adapt 1 pair lj/cut epsilon * * v_scale1 pair coul/cut scale 3 3 v_scale2 scale yes reset yes
   fix 1 all adapt 10 atom diameter v_size

   variable ramp_up equal "ramp(0.01,0.5)"
   fix stretch all adapt 1 bond harmonic r0 1 v_ramp_up

Description
"""""""""""

Change or adapt one or more specific simulation attributes or settings over
time as a simulation runs.  Pair potential and :math:`k`-space and atom
attributes which can be varied by this fix are discussed below.  Many other
fixes can also be used to time-vary simulation parameters (e.g., the
:doc:`fix deform <fix_deform>` command will change the simulation box
size/shape and the :doc:`fix move <fix_move>` command will change atom
positions and velocities in a prescribed manner).  Also note that many commands
allow variables as arguments for specific parameters, if described in that
manner on their doc pages.  An equal-style variable can calculate a
time-dependent quantity, so this is another way to vary a simulation parameter
over time.

If :math:`N` is specified as 0, the specified attributes are only changed
once, before the simulation begins.  This is all that is needed if the
associated variables are not time-dependent.  If :math:`N > 0`, then changes
are made every :math:`N` steps during the simulation, presumably with a
variable that is time-dependent.

Depending on the value of the *reset* keyword, attributes changed by
this fix will or will not be reset back to their original values at
the end of a simulation.  Even if *reset* is specified as *yes*, a
restart file written during a simulation will contain the modified
settings.

If the *scale* keyword is set to *no*, which is the default, then
the value of the altered parameter will be whatever the variable
generates.  If the *scale* keyword is set to *yes*, then the value
of the altered parameter will be the initial value of that parameter
multiplied by whatever the variable generates (i.e., the variable is
now a "scale factor" applied in (presumably) a time-varying fashion to
the parameter).

Note that whether scale is *no* or *yes*, internally, the parameters
themselves are actually altered by this fix.  Make sure you use the
*reset yes* option if you want the parameters to be restored to their
initial values after the run.

----------

The *pair* keyword enables various parameters of potentials defined by
the :doc:`pair_style <pair_style>` command to be changed, if the pair
style supports it.  Note that the :doc:`pair_style <pair_style>` and
:doc:`pair_coeff <pair_coeff>` commands must be used in the usual manner
to specify these parameters initially; the fix adapt command simply
overrides the parameters.

The *pstyle* argument is the name of the pair style.  If
:doc:`pair_style hybrid or hybrid/overlay <pair_hybrid>` is used,
*pstyle* should be a sub-style name.  If there are multiple
sub-styles using the same pair style, then *pstyle* should be specified
as "style:N", where *N* is which instance of the pair style you wish to
adapt (e.g., the first or second).  For example, *pstyle* could be
specified as "soft" or "lubricate" or "lj/cut:1" or "lj/cut:2".  The
*pparam* argument is the name of the parameter to change.  This is the
current list of pair styles and parameters that can be varied by this
fix.  See the doc pages for individual pair styles and their energy
formulas for the meaning of these parameters:

+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`born <pair_born>`                                                      | a,b,c                                            | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`born/coul/long, born/coul/msm <pair_born>`                             | coulombic_cutoff                                 | type global |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`born/gauss <pair_born_gauss>`                                          | biga0,biga1,r0                                   | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`buck, buck/coul/cut  <pair_buck>`                                      | a,c                                              | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`buck/coul/long, buck/coul/msm <pair_buck>`                             | a,c,coulombic_cutoff                             | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`buck/mdf <pair_mdf>`                                                   | a,c                                              | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`coul/cut, coul/cut/global <pair_coul>`                                 | scale                                            | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`coul/cut/soft <pair_fep_soft>`                                         | lambda                                           | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`coul/debye <pair_coul>`                                                | scale                                            | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`coul/dsf <pair_coul>`                                                  | coulombic_cutoff                                 | type global |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`coul/long, coul/msm <pair_coul>`                                       | coulombic_cutoff, scale                          | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`coul/long/soft <pair_fep_soft>`                                        | scale, lambda, coulombic_cutoff                  | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`coul/slater/long <pair_coul_slater>`                                   | scale                                            | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`coul/streitz <pair_coul>`                                              | scale                                            | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`eam, eam/alloy, eam/fs <pair_eam>`                                     | scale                                            | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`gauss <pair_gauss>`                                                    | a                                                | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`harmonic/cut <pair_harmonic_cut>`                                      | k, cutoff                                        | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`kim <pair_kim>`                                                        | scale                                            | type global |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lennard/mdf <pair_mdf>`                                                | A,B                                              | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/class2 <pair_class2>`                                               | epsilon,sigma                                    | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/class2/coul/cut, lj/class2/coul/long <pair_class2>`                 | epsilon,sigma,coulombic_cutoff                   | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/cut <pair_lj>`                                                      | epsilon,sigma                                    | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/cut/coul/cut, lj/cut/coul/long, lj/cut/coul/msm <pair_lj_cut_coul>` | epsilon,sigma,coulombic_cutoff                   | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/cut/coul/cut/soft, lj/cut/coul/long/soft <pair_fep_soft>`           | epsilon,sigma,lambda,coulombic_cutoff            | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/cut/coul/dsf <pair_lj_cut_coul>`                                    | cutoff                                           | type global |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/cut/tip4p/cut <pair_lj_cut_tip4p>`                                  | epsilon,sigma,coulombic_cutoff                   | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/cut/soft <pair_fep_soft>`                                           | epsilon,sigma,lambda                             | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/expand <pair_lj_expand>`                                            | epsilon,sigma,delta                              | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/mdf <pair_mdf>`                                                     | epsilon,sigma                                    | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lj/sf/dipole/sf <pair_dipole>`                                         | epsilon,sigma,scale                              | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`lubricate <pair_lubricate>`                                            | mu                                               | global      |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`meam <pair_meam>`                                                      | scale                                            | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`mie/cut <pair_mie>`                                                    | epsilon,sigma,gamma_repulsive,gamma_attractive   | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`morse, morse/smooth/linear <pair_morse>`                               | D0,R0,alpha                                      | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`morse/soft <pair_morse>`                                               | D0,R0,alpha,lambda                               | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`nm/cut <pair_nm>`                                                      | E0,R0,m,n                                        | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`nm/cut/coul/cut, nm/cut/coul/long <pair_nm>`                           | E0,R0,m,n,coulombic_cutoff                       | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`pace, pace/extrapolation <pair_pace>`                                  | scale                                            | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`snap <pair_snap>`                                                      | scale                                            | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`spin/dmi <pair_spin_dmi>`                                              | coulombic_cutoff                                 | type global |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`spin/exchange <pair_spin_exchange>`                                    | coulombic_cutoff                                 | type global |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`spin/magelec <pair_spin_magelec>`                                      | coulombic_cutoff                                 | type global |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`spin/neel <pair_spin_neel>`                                            | coulombic_cutoff                                 | type global |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`soft <pair_soft>`                                                      | a                                                | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`table <pair_table>`                                                    | table_cutoff                                     | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`ufm <pair_ufm>`                                                        | epsilon,sigma                                    | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+
| :doc:`wf/cut <pair_wf_cut>`                                                  | epsilon,sigma,nu,mu                              | type pairs  |
+------------------------------------------------------------------------------+--------------------------------------------------+-------------+

.. note::

   It is easy to add new pairwise potentials and their parameters
   to this list.  All it typically takes is adding an extract() method to
   the pair\_\*.cpp file associated with the potential.

Some parameters are global settings for the pair style (e.g., the
viscosity setting "mu" for :doc:`pair_style lubricate <pair_lubricate>`).
Other parameters apply to atom type pairs within the pair style (e.g., the
prefactor :math:`a` for :doc:`pair_style soft <pair_soft>`).

Note that for many of the potentials, the parameter that can be varied
is effectively a prefactor on the entire energy expression for the
potential (e.g., the lj/cut epsilon).  The parameters listed as "scale"
are exactly that, since the energy expression for the
:doc:`coul/cut <pair_coul>` potential (for example) has no labeled
prefactor in its formula.  To apply an effective prefactor to some
potentials, multiple parameters need to be altered.  For example, the
:doc:`Buckingham potential <pair_buck>` needs both the :math:`A` and
:math:`C` terms altered together.  To scale the Buckingham potential, you
should thus list the pair style twice, once for :math:`A` and once for
:math:`C`.

If a type pair parameter is specified, the :math:`I` and :math:`J` settings
should be specified to indicate which type pairs to apply it to.  If a global
parameter is specified, the :math:`I` and :math:`J` settings still need to be
specified, but are ignored.

Similar to the :doc:`pair_coeff command <pair_coeff>`, :math:`I` and :math:`J`
can be specified in one of two ways.  Explicit numeric values can be used for
each, as in the first example above.  :math:`I \le J` is required.  LAMMPS sets
the coefficients for the symmetric :math:`J,I` interaction to the same values.

A wild-card asterisk can be used in place of or in conjunction with
the :math:`I,J` arguments to set the coefficients for multiple pairs of atom
types.  This takes the form "\*" or "\*n" or "m\*" or "m\*n".  If :math:`N`
is the number of atom types, then an asterisk with no numeric values
means all types from 1 to :math:`N`.  A leading asterisk means all types from
1 to n (inclusive).  A trailing asterisk means all types from m to :math:`N`
(inclusive).  A middle asterisk means all types from m to n
(inclusive).  Note that only type pairs with :math:`I \le J` are considered; if
asterisks imply type pairs where :math:`J < I`, they are ignored.

IMPORTANT NOTE: If :doc:`pair_style hybrid or hybrid/overlay
<pair_hybrid>` is being used, then the *pstyle* will be a sub-style
name.  You must specify :math:`I,J` arguments that correspond to type pair
values defined (via the :doc:`pair_coeff <pair_coeff>` command) for
that sub-style.

The *v_name* argument for keyword *pair* is the name of an
:doc:`equal-style variable <variable>` which will be evaluated each time
this fix is invoked to set the parameter to a new value.  It should be
specified as v_name, where name is the variable name.  Equal-style
variables can specify formulas with various mathematical functions, and
include :doc:`thermo_style <thermo_style>` command keywords for the
simulation box parameters and timestep and elapsed time.  Thus it is
easy to specify parameters that change as a function of time or span
consecutive runs in a continuous fashion.  For the latter, see the
*start* and *stop* keywords of the :doc:`run <run>` command and the
*elaplong* keyword of :doc:`thermo_style custom <thermo_style>` for
details.

For example, these commands would change the prefactor coefficient of
the :doc:`pair_style soft <pair_soft>` potential from 10.0 to 30.0 in a
linear fashion over the course of a simulation:

.. code-block:: LAMMPS

   variable prefactor equal ramp(10,30)
   fix 1 all adapt 1 pair soft a * * v_prefactor

----------

The *bond* keyword uses the specified variable to change the value of
a bond coefficient over time, very similar to how the *pair* keyword
operates. The only difference is that now a bond coefficient for a
given bond type is adapted.

A wild-card asterisk can be used in place of or in conjunction with the
bond type argument to set the coefficients for multiple bond types.
This takes the form "\*" or "\*n" or "m\*" or "m\*n".  If :math:`N` is
the number of bond types, then an asterisk with no numeric values means
all types from 1 to :math:`N`.  A leading asterisk means all types from
1 to n (inclusive).  A trailing asterisk means all types from m to
:math:`N` (inclusive).  A middle asterisk means all types from m to n
(inclusive).

Currently *bond* does not support bond_style hybrid nor bond_style
hybrid/overlay as bond styles. The bond styles that currently work
with fix_adapt are

+------------------------------------+-------+-----------------+
| :doc:`class2 <bond_class2>`        | r0         | type bonds |
+------------------------------------+-------+-----------------+
| :doc:`fene <bond_fene>`            | k,r0       | type bonds |
+------------------------------------+-------+-----------------+
| :doc:`fene/nm <bond_fene>`         | k,r0       | type bonds |
+------------------------------------+-------+-----------------+
| :doc:`gromos <bond_gromos>`        | k,r0       | type bonds |
+------------------------------------+-------+-----------------+
| :doc:`harmonic <bond_harmonic>`    | k,r0       | type bonds |
+------------------------------------+-------+-----------------+
| :doc:`morse <bond_morse>`          | r0         | type bonds |
+------------------------------------+-------+-----------------+
| :doc:`nonlinear <bond_nonlinear>`  | epsilon,r0 | type bonds |
+------------------------------------+-------+-----------------+

----------

.. versionadded:: 4May2022

The *angle* keyword uses the specified variable to change the value of
an angle coefficient over time, very similar to how the *pair* keyword
operates. The only difference is that now an angle coefficient for a
given angle type is adapted.

A wild-card asterisk can be used in place of or in conjunction with the
angle type argument to set the coefficients for multiple angle types.
This takes the form "\*" or "\*n" or "m\*" or "m\*n".  If :math:`N` is
the number of angle types, then an asterisk with no numeric values means
all types from 1 to :math:`N`.  A leading asterisk means all types from
1 to n (inclusive).  A trailing asterisk means all types from m to
:math:`N` (inclusive).  A middle asterisk means all types from m to n
(inclusive).

Currently *angle* does not support angle_style hybrid nor angle_style
hybrid/overlay as angle styles. The angle styles that currently work
with fix_adapt are

+------------------------------------+-------+-----------------+
| :doc:`harmonic <angle_harmonic>`    | k,theta0 | type angles |
+------------------------------------+-------+-----------------+
| :doc:`cosine <angle_cosine>`        | k        | type angles |
+------------------------------------+-------+-----------------+

Note that internally, theta0 is stored in radians, so the variable
this fix uses to reset theta0 needs to generate values in radians.

----------

The *kspace* keyword used the specified variable as a scale factor on
the energy, forces, virial calculated by whatever :math:`k`-space solver is
defined by the :doc:`kspace_style <kspace_style>` command.  If the
variable has a value of 1.0, then the solver is unaltered.

The *kspace* keyword works this way whether the *scale* keyword
is set to *no* or *yes*\ .

----------

The *atom* keyword enables various atom properties to be changed.  The
*aparam* argument is the name of the parameter to change.  This is the
current list of atom parameters that can be varied by this fix:

* charge = charge on particle
* diameter or diameter/disc = diameter of particle

The *v_name* argument of the *atom* keyword is the name of an
:doc:`equal-style variable <variable>` which will be evaluated each
time this fix is invoked to set, or scale the parameter to a new
value.  It should be specified as v_name, where name is the variable
name.  See the discussion above describing the formulas associated
with equal-style variables.  The new value is assigned to the
corresponding attribute for all atoms in the fix group.

If the atom parameter is *diameter* and per-atom density and per-atom
mass are defined for particles (e.g., :doc:`atom_style granular
<atom_style>`), then the mass of each particle is, by default, also
changed when the diameter changes. The mass is set from the particle
volume for 3d systems (density is assumed to stay constant). For 2d,
the default is for LAMMPS to model particles with a radius attribute
as spheres. However, if the atom parameter is *diameter/disc*, then the
mass is set from the particle area (the density is assumed to be in
mass/distance\ :math:`^2` units). The mass of the particle may also be kept
constant if the *mass* keyword is set to *no*. This can be useful to account
for diameter changes that do not involve mass changes (e.g., thermal
expansion).

For example, these commands would shrink the diameter of all granular
particles in the "center" group from 1.0 to 0.1 in a linear fashion
over the course of a 1000-step simulation:

.. code-block:: LAMMPS

   variable size equal ramp(1.0,0.1)
   fix 1 center adapt 10 atom diameter v_size

----------

This fix can be used in long simulations which are restarted one or
more times to continuously adapt simulation parameters, but it must be
done carefully.  There are two issues to consider.  The first is how
to adapt the parameters in a continuous manner from one simulation to
the next.  The second is how, if desired, to reset the parameters to
their original values at the end of the last restarted run.

Note that all the parameters changed by this fix are written into a
restart file in their current changed state.  A new restarted
simulation does not know the original time=0 values, unless the
input script explicitly resets the parameters (after the restart file
is read) to the original values.

Also note that the time-dependent variable(s) used in the restart
script should typically be written as a function of time elapsed since
the original simulation began.

With this in mind, if the *scale* keyword is set to *no* (the default)
in a restarted simulation, original parameters are not needed.  The
adapted parameters should seamlessly continue their variation relative
to the preceding simulation.

If the *scale* keyword is set to *yes*, then the input script should
typically reset the parameters being adapted to their original values,
so that the scaling formula specified by the variable will operate
correctly.  An exception is if the *atom* keyword is being used with
*scale yes*.  In this case, information is added to the restart file
so that per-atom properties in the new run will automatically be
scaled relative to their original values.  This will only work if the
fix adapt command specified in the restart script has the same ID as
the one used in the original script.

In a restarted run, if the *reset* keyword is set to *yes*, and the
run ends in this script (as opposed to just writing more restart
files), parameters will be restored to the values they were at the
beginning of the run command in the restart script, which as
explained above, may or may not be the original values of the
parameters.  Again, an exception is if the *atom* keyword is being
used with *reset yes* (in all the runs). In that case, the original
per-atom parameters are stored in the restart file, and will be
restored when the restarted run finally completes.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If the *atom* keyword is used and the *scale* or *reset* keyword is
set to *yes*, then this fix writes information to a restart file so
that in a restarted run scaling can continue in a seamless manner
and/or the per-atom values can be restored, as explained above.

None of the :doc:`fix_modify <fix_modify>` options are relevant to
this fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.  No parameter
of this fix can be used with the *start/stop* keywords of the
:doc:`run <run>` command.  This fix is not invoked during :doc:`energy
minimization <minimize>`.

For :doc:`rRESPA time integration <run_style>`, this fix changes
parameters on the outermost rRESPA level.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute ti <compute_ti>`

Default
"""""""

The option defaults are scale = no, reset = no, mass = yes.
