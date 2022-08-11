.. index:: pair_style brownian
.. index:: pair_style brownian/omp
.. index:: pair_style brownian/poly
.. index:: pair_style brownian/poly/omp

pair_style brownian command
===========================

Accelerator Variants: *brownian/omp*

pair_style brownian/poly command
================================

Accelerator Variants: *brownian/poly/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style mu flaglog flagfld cutinner cutoff t_target seed flagHI flagVF

* style = *brownian* or *brownian/poly*
* mu = dynamic viscosity (dynamic viscosity units)
* flaglog = 0/1 log terms in the lubrication approximation on/off
* flagfld = 0/1 to include/exclude Fast Lubrication Dynamics effects
* cutinner = inner cutoff distance (distance units)
* cutoff = outer cutoff for interactions (distance units)
* t_target = target temp of the system (temperature units)
* seed = seed for the random number generator (positive integer)
* flagHI (optional) = 0/1 to include/exclude 1/r hydrodynamic interactions
* flagVF (optional) = 0/1 to include/exclude volume fraction corrections in the long-range isotropic terms

Examples
""""""""

.. code-block:: LAMMPS

   pair_style brownian 1.5 1 1 2.01 2.5 2.0 5878567 # (assuming radius = 1)
   pair_coeff 1 1 2.05 2.8
   pair_coeff * *

Description
"""""""""""

Styles *brownian* and *brownian/poly* compute Brownian forces and
torques on finite-size spherical particles.  The former requires
monodisperse spherical particles; the latter allows for polydisperse
spherical particles.

These pair styles are designed to be used with either the :doc:`pair_style lubricate <pair_lubricate>` or :doc:`pair_style lubricateU <pair_lubricateU>` commands to provide thermostatting
when dissipative lubrication forces are acting.  Thus the parameters
*mu*, *flaglog*, *flagfld*, *cutinner*, and *cutoff* should be
specified consistent with the settings in the lubrication pair styles.
For details, refer to either of the lubrication pair styles.

The *t_target* setting is used to specify the target temperature of
the system.  The random number *seed* is used to generate random
numbers for the thermostatting procedure.

The *flagHI* and *flagVF* settings are optional.  Neither should be
used, or both must be defined.

----------

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* cutinner (distance units)
* cutoff (distance units)

The two coefficients are optional.  If neither is specified, the two
cutoffs specified in the pair_style command are used.  Otherwise both
must be specified.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the two cutoff distances for this
pair style can be mixed.  The default mix value is *geometric*\ .  See
the "pair_modify" command for details.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

These styles are part of the COLLOID package.  They are only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Only spherical monodisperse particles are allowed for pair_style
brownian.

Only spherical particles are allowed for pair_style brownian/poly.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style lubricate <pair_lubricate>`, :doc:`pair_style lubricateU <pair_lubricateU>`

Default
"""""""

The default settings for the optional args are flagHI = 1 and flagVF =
1.
