.. index:: fix adapt/fep

fix adapt/fep command
=====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID adapt/fep N attribute args ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* adapt/fep = style name of this fix command
* N = adapt simulation settings every this many timesteps
* one or more attribute/arg pairs may be appended
* attribute = *pair* or *kspace* or *atom*

  .. parsed-literal::

       *pair* args = pstyle pparam I J v_name
         pstyle = pair style name, e.g. lj/cut
         pparam = parameter to adapt over time
         I,J = type pair(s) to set parameter for
         v_name = variable with name that calculates value of pparam
       *kspace* arg = v_name
         v_name = variable with name that calculates scale factor on K-space terms
       *atom* args = aparam v_name
         aparam = parameter to adapt over time
         I = type(s) to set parameter for
         v_name = variable with name that calculates value of aparam

* zero or more keyword/value pairs may be appended
* keyword = *scale* or *reset* or *after*

  .. parsed-literal::

       *scale* value = *no* or *yes*
         *no* = the variable value is the new setting
         *yes* = the variable value multiplies the original setting
       *reset* value = *no* or *yes*
         *no* = values will remain altered at the end of a run
         *yes* = reset altered values to their original values at the end
         of a run
       *after* value = *no* or *yes*
         *no* = parameters are adapted at timestep N
         *yes* = parameters are adapted one timestep after N



Examples
""""""""


.. parsed-literal::

   fix 1 all adapt/fep 1 pair soft a 1 1 v_prefactor
   fix 1 all adapt/fep 1 pair soft a 2\* 3 v_prefactor
   fix 1 all adapt/fep 1 pair lj/cut epsilon \* \* v_scale1 coul/cut scale 3 3 v_scale2 scale yes reset yes
   fix 1 all adapt/fep 10 atom diameter 1 v_size

Description
"""""""""""

Change or adapt one or more specific simulation attributes or settings
over time as a simulation runs.

This is an enhanced version of the :doc:`fix adapt <fix_adapt>` command
with two differences,

* It is possible to modify the charges of chosen atom types only,
  instead of scaling all the charges in the system.
* There is a new option *after* for better compatibility with "fix
  ave/time".


This version is suited for free energy calculations using
:doc:`compute ti <compute_ti>` or :doc:`compute fep <compute_fep>`.

If *N* is specified as 0, the specified attributes are only changed
once, before the simulation begins.  This is all that is needed if the
associated variables are not time-dependent.  If *N* > 0, then changes
are made every *N* steps during the simulation, presumably with a
variable that is time-dependent.

Depending on the value of the *reset* keyword, attributes changed by
this fix will or will not be reset back to their original values at
the end of a simulation.  Even if *reset* is specified as *yes*\ , a
restart file written during a simulation will contain the modified
settings.

If the *scale* keyword is set to *no*\ , then the value the parameter is
set to will be whatever the variable generates.  If the *scale*
keyword is set to *yes*\ , then the value of the altered parameter will
be the initial value of that parameter multiplied by whatever the
variable generates.  I.e. the variable is now a "scale factor" applied
in (presumably) a time-varying fashion to the parameter.  Internally,
the parameters themselves are actually altered; make sure you use the
*reset yes* option if you want the parameters to be restored to their
initial values after the run.

If the *after* keyword is set to *yes*\ , then the parameters are
changed one timestep after the multiple of N. In this manner, if a fix
such as "fix ave/time" is used to calculate averages at every N
timesteps, all the contributions to the average will be obtained with
the same values of the parameters.


----------


The *pair* keyword enables various parameters of potentials defined by
the :doc:`pair_style <pair_style>` command to be changed, if the pair
style supports it.  Note that the :doc:`pair_style <pair_style>` and
:doc:`pair_coeff <pair_coeff>` commands must be used in the usual manner
to specify these parameters initially; the fix adapt command simply
overrides the parameters.

The *pstyle* argument is the name of the pair style.  If :doc:`pair_style hybrid or hybrid/overlay <pair_hybrid>` is used, *pstyle* should be
a sub-style name.  For example, *pstyle* could be specified as "soft"
or "lubricate".  The *pparam* argument is the name of the parameter to
change.  This is the current list of pair styles and parameters that
can be varied by this fix.  See the doc pages for individual pair
styles and their energy formulas for the meaning of these parameters:

+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`born <pair_born>`                                             | a,b,c                   | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`buck <pair_buck>`                                             | a,c                     | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`buck/mdf <pair_mdf>`                                          | a,c                     | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`coul/cut <pair_coul>`                                         | scale                   | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`coul/cut/soft <pair_fep_soft>`                                | lambda                  | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`coul/long, coul/msm <pair_coul>`                              | scale                   | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`coul/long/soft <pair_fep_soft>`                               | scale, lambda           | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`eam <pair_eam>`                                               | scale                   | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`gauss <pair_gauss>`                                           | a                       | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lennard/mdf <pair_mdf>`                                       | a,b                     | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/class2 <pair_class2>`                                      | epsilon,sigma           | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/class2/coul/cut, lj/class2/coul/long <pair_class2>`        | epsilon,sigma           | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/cut <pair_lj>`                                             | epsilon,sigma           | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/cut/soft <pair_fep_soft>`                                  | epsilon,sigma,lambda    | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/cut/coul/cut, lj/cut/coul/long, lj/cut/coul/msm <pair_lj>` | epsilon,sigma           | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/cut/coul/cut/soft, lj/cut/coul/long/soft <pair_fep_soft>`  | epsilon,sigma,lambda    | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/cut/tip4p/cut, lj/cut/tip4p/long <pair_lj>`                | epsilon,sigma           | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/cut/tip4p/long/soft <pair_fep_soft>`                       | epsilon,sigma,lambda    | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/expand <pair_lj_expand>`                                   | epsilon,sigma,delta     | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/mdf <pair_mdf>`                                            | epsilon,sigma           | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`lj/sf/dipole/sf <pair_dipole>`                                | epsilon,sigma,scale     | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`mie/cut <pair_mie>`                                           | epsilon,sigma,gamR,gamA | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`morse, morse/smooth/linear <pair_morse>`                      | d0,r0,alpha             | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`morse/soft <pair_morse>`                                      | d0,r0,alpha,lambda      | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`nm/cut <pair_nm>`                                             | e0,r0,nn,mm             | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`nm/cut/coul/cut, nm/cut/coul/long <pair_nm>`                  | e0,r0,nn,mm             | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`ufm <pair_ufm>`                                               | epsilon,sigma,scale     | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+
| :doc:`soft <pair_soft>`                                             | a                       | type pairs |
+---------------------------------------------------------------------+-------------------------+------------+

.. note::

   It is easy to add new potentials and their parameters to this
   list.  All it typically takes is adding an extract() method to the
   pair\_\*.cpp file associated with the potential.

Note that for many of the potentials, the parameter that can be varied
is effectively a prefactor on the entire energy expression for the
potential, e.g. the lj/cut epsilon.  The parameters listed as "scale"
are exactly that, since the energy expression for the
:doc:`coul/cut <pair_coul>` potential (for example) has no labeled
prefactor in its formula.  To apply an effective prefactor to some
potentials, multiple parameters need to be altered.  For example, the
:doc:`Buckingham potential <pair_buck>` needs both the A and C terms
altered together.  To scale the Buckingham potential, you should thus
list the pair style twice, once for A and once for C.

If a type pair parameter is specified, the *I* and *J* settings should
be specified to indicate which type pairs to apply it to.  If a global
parameter is specified, the *I* and *J* settings still need to be
specified, but are ignored.

Similar to the :doc:`pair_coeff command <pair_coeff>`, I and J can be
specified in one of two ways.  Explicit numeric values can be used for
each, as in the 1st example above.  I <= J is required.  LAMMPS sets
the coefficients for the symmetric J,I interaction to the same values.

A wild-card asterisk can be used in place of or in conjunction with
the I,J arguments to set the coefficients for multiple pairs of atom
types.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the
number of atom types, then an asterisk with no numeric values means
all types from 1 to N.  A leading asterisk means all types from 1 to n
(inclusive).  A trailing asterisk means all types from n to N
(inclusive).  A middle asterisk means all types from m to n
(inclusive).  Note that only type pairs with I <= J are considered; if
asterisks imply type pairs where J < I, they are ignored.

IMPROTANT NOTE: If :doc:`pair_style hybrid or hybrid/overlay <pair_hybrid>` is being used, then the *pstyle* will
be a sub-style name.  You must specify I,J arguments that correspond
to type pair values defined (via the :doc:`pair_coeff <pair_coeff>`
command) for that sub-style.

The *v\_name* argument for keyword *pair* is the name of an
:doc:`equal-style variable <variable>` which will be evaluated each time
this fix is invoked to set the parameter to a new value.  It should be
specified as v\_name, where name is the variable name.  Equal-style
variables can specify formulas with various mathematical functions,
and include :doc:`thermo_style <thermo_style>` command keywords for the
simulation box parameters and timestep and elapsed time.  Thus it is
easy to specify parameters that change as a function of time or span
consecutive runs in a continuous fashion.  For the latter, see the
*start* and *stop* keywords of the :doc:`run <run>` command and the
*elaplong* keyword of :doc:`thermo_style custom <thermo_style>` for
details.

For example, these commands would change the prefactor coefficient of
the :doc:`pair_style soft <pair_soft>` potential from 10.0 to 30.0 in a
linear fashion over the course of a simulation:


.. parsed-literal::

   variable prefactor equal ramp(10,30)
   fix 1 all adapt 1 pair soft a \* \* v_prefactor


----------


The *kspace* keyword used the specified variable as a scale factor on
the energy, forces, virial calculated by whatever K-Space solver is
defined by the :doc:`kspace_style <kspace_style>` command.  If the
variable has a value of 1.0, then the solver is unaltered.

The *kspace* keyword works this way whether the *scale* keyword
is set to *no* or *yes*\ .


----------


The *atom* keyword enables various atom properties to be changed.  The
*aparam* argument is the name of the parameter to change.  This is the
current list of atom parameters that can be varied by this fix:

* charge = charge on particle
* diameter = diameter of particle

The *I* argument indicates which atom types are affected. A wild-card
asterisk can be used in place of or in conjunction with the I argument
to set the coefficients for multiple atom types.

The *v\_name* argument of the *atom* keyword is the name of an
:doc:`equal-style variable <variable>` which will be evaluated each time
this fix is invoked to set the parameter to a new value.  It should be
specified as v\_name, where name is the variable name.  See the
discussion above describing the formulas associated with equal-style
variables.  The new value is assigned to the corresponding attribute
for all atoms in the fix group.

If the atom parameter is *diameter* and per-atom density and per-atom
mass are defined for particles (e.g. :doc:`atom_style granular <atom_style>`), then the mass of each particle is also
changed when the diameter changes (density is assumed to stay
constant).

For example, these commands would shrink the diameter of all granular
particles in the "center" group from 1.0 to 0.1 in a linear fashion
over the course of a 1000-step simulation:


.. parsed-literal::

   variable size equal ramp(1.0,0.1)
   fix 1 center adapt 10 atom diameter \* v_size

For :doc:`rRESPA time integration <run_style>`, this fix changes
parameters on the outermost rRESPA level.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute fep <compute_fep>`, :doc:`fix adapt <fix_adapt>`, :doc:`compute ti <compute_ti>`, :doc:`pair_style \*/soft <pair_fep_soft>`

Default
"""""""

The option defaults are scale = no, reset = no, after = no.
