.. index:: pair\_style brownian

pair\_style brownian command
============================

pair\_style brownian/omp command
================================

pair\_style brownian/poly command
=================================

pair\_style brownian/poly/omp command
=====================================

Syntax
""""""


.. parsed-literal::

   pair_style style mu flaglog flagfld cutinner cutoff t_target seed flagHI flagVF

* style = *brownian* or *brownian/poly*
* mu = dynamic viscosity (dynamic viscosity units)
* flaglog = 0/1 log terms in the lubrication approximation on/off
* flagfld = 0/1 to include/exclude Fast Lubrication Dynamics effects
* cutinner = inner cutoff distance (distance units)
* cutoff = outer cutoff for interactions (distance units)
* t\_target = target temp of the system (temperature units)
* seed = seed for the random number generator (positive integer)
* flagHI (optional) = 0/1 to include/exclude 1/r hydrodynamic interactions
* flagVF (optional) = 0/1 to include/exclude volume fraction corrections in the long-range isotropic terms

Examples
""""""""


.. parsed-literal::

   pair_style brownian 1.5 1 1 2.01 2.5 2.0 5878567 (assuming radius = 1)
   pair_coeff 1 1 2.05 2.8
   pair_coeff \* \*

Description
"""""""""""

Styles *brownian* and *brownian/poly* compute Brownian forces and
torques on finite-size spherical particles.  The former requires
monodisperse spherical particles; the latter allows for polydisperse
spherical particles.

These pair styles are designed to be used with either the :doc:`pair\_style lubricate <pair_lubricate>` or :doc:`pair\_style lubricateU <pair_lubricateU>` commands to provide thermostatting
when dissipative lubrication forces are acting.  Thus the parameters
*mu*\ , *flaglog*\ , *flagfld*\ , *cutinner*\ , and *cutoff* should be
specified consistent with the settings in the lubrication pair styles.
For details, refer to either of the lubrication pair styles.

The *t\_target* setting is used to specify the target temperature of
the system.  The random number *seed* is used to generate random
numbers for the thermostatting procedure.

The *flagHI* and *flagVF* settings are optional.  Neither should be
used, or both must be defined.


----------


The following coefficients must be defined for each pair of atoms
types via the :doc:`pair\_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands, or by mixing as described below:

* cutinner (distance units)
* cutoff (distance units)

The two coefficients are optional.  If neither is specified, the two
cutoffs specified in the pair\_style command are used.  Otherwise both
must be specified.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed in :doc:`this section <Speed>` of
the manual.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See :doc:`this section <Speed>` of the manual for more
instructions on how to use the accelerated styles effectively.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, the two cutoff distances for this
pair style can be mixed.  The default mix value is *geometric*\ .  See
the "pair\_modify" command for details.

This pair style does not support the :doc:`pair\_modify <pair_modify>`
shift option for the energy of the pair interaction.

The :doc:`pair\_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair\_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run\_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


These styles are part of the COLLOID package.  They are only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Only spherical monodisperse particles are allowed for pair\_style
brownian.

Only spherical particles are allowed for pair\_style brownian/poly.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`, :doc:`pair\_style lubricate <pair_lubricate>`, :doc:`pair\_style lubricateU <pair_lubricateU>`

Default
"""""""

The default settings for the optional args are flagHI = 1 and flagVF =
1.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
