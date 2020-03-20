.. index:: fix colvars

fix colvars command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID colvars configfile keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* colvars = style name of this fix command
* configfile = the configuration file for the colvars module
* keyword = *input* or *output* or *seed* or *tstat*

  .. parsed-literal::

       *input* arg = colvars.state file name or prefix or NULL (default: NULL)
       *output* arg = output filename prefix (default: out)
       *seed* arg = seed for random number generator (default: 1966)
       *unwrap* arg = *yes* or *no*
         use unwrapped coordinates in collective variables (default: yes)
       *tstat* arg = fix id of a thermostat or NULL (default: NULL)

Examples
""""""""

.. code-block:: LAMMPS

   fix mtd all colvars peptide.colvars.inp seed 2122 input peptide.colvars.state output peptide
   fix abf all colvars colvars.inp tstat 1

Description
"""""""""""

This fix interfaces LAMMPS to the collective variables "Colvars"
library, which allows to calculate potentials of mean force
(PMFs) for any set of colvars, using different sampling methods:
currently implemented are the Adaptive Biasing Force (ABF) method,
metadynamics, Steered Molecular Dynamics (SMD) and Umbrella Sampling
(US) via a flexible harmonic restraint bias.

This documentation describes only the fix colvars command itself and
LAMMPS specific parts of the code.  The full documentation of the
colvars library is available as `this supplementary PDF document <PDF/colvars-refman-lammps.pdf>`_

The Colvars library is developed at `https://github.com/colvars/colvars <https://github.com/colvars/colvars>`_
A detailed discussion of its implementation is in :ref:`(Fiorin) <Fiorin>`.

There are some example scripts for using this package with LAMMPS in the
examples/USER/colvars directory.

----------

The only mandatory argument to the fix is the filename to the colvars
input file that contains the input that is independent from the MD
program in which the colvars library has been integrated.

The *group-ID* entry is ignored. The collective variable module will
always apply to the entire system and there can only be one instance
of the colvars fix at a time. The colvars fix will only communicate
the minimum information necessary and the colvars library supports
multiple, completely independent collective variables, so there is
no restriction to functionality by limiting the number of colvars fixes.

The *input* keyword allows to specify a state file that would contain
the restart information required in order to continue a calculation from
a prerecorded state. Fix colvars records it state in :doc:`binary restart <restart>`
files, so when using the :doc:`read_restart <read_restart>` command,
this is usually not needed.

The *output* keyword allows to specify the output prefix. All output
files generated will use this prefix followed by the ".colvars." and
a word like "state" or "traj".

The *seed* keyword contains the seed for the random number generator
that will be used in the colvars module.

The *unwrap* keyword controls whether wrapped or unwrapped coordinates
are passed to the colvars library for calculation of the collective
variables and the resulting forces. The default is *yes*\ , i.e. to use
the image flags to reconstruct the absolute atom positions.
Setting this to *no* will use the current local coordinates that are
wrapped back into the simulation cell at each re-neighboring instead.

The *tstat* keyword can be either NULL or the label of a thermostatting
fix that thermostats all atoms in the fix colvars group. This will be
used to provide the colvars module with the current thermostat target
temperature.

**Restart, fix_modify, output, run start/stop, minimize info:**

This fix writes the current status of the colvars module into
:doc:`binary restart files <restart>`. This is in addition to the text
mode status file that is written by the colvars module itself and the
kind of information in both files is identical.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy change from the biasing force added by the fix
to the system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.

The *fix_modify configfile <config file>* option allows to add settings
from an additional config file to the colvars module. This option can
only be used, after the system has been initialized with a :doc:`run <run>`
command.

The *fix_modify config <quoted string>* option allows to add settings
from inline strings. Those have to fit on a single line when enclosed
in a pair of double quotes ("), or can span multiple lines when bracketed
by a pair of triple double quotes (""", like python embedded documentation).

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the cumulative
energy change due to this fix.  The scalar value calculated by this
fix is "extensive".

Restrictions
""""""""""""

This fix is part of the USER-COLVARS package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

There can only be one colvars fix active at a time. Since the interface
communicates only the minimum amount of information and colvars module
itself can handle an arbitrary number of collective variables, this is
not a limitation of functionality.

Related commands
""""""""""""""""

:doc:`fix smd <fix_smd>`, :doc:`fix spring <fix_spring>`,
:doc:`fix plumed <fix_plumed>`

Default
"""""""

The default options are input = NULL, output = out, seed = 1966, unwrap yes,
and tstat = NULL.

----------

.. _Fiorin:

**(Fiorin)** Fiorin, Klein, Henin, Mol. Phys., DOI:10.1080/00268976.2013.813594
