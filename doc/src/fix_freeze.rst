.. index:: fix freeze
.. index:: fix freeze/kk

fix freeze command
==================

Accelerator Variants: *freeze/kk*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID freeze

* ID, group-ID are documented in :doc:`fix <fix>` command
* freeze = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 bottom freeze

Description
"""""""""""

Zero out the force and torque on a granular particle.  This is useful
for preventing certain particles from moving in a simulation.  The
:doc:`granular pair styles <pair_gran>` also detect if this fix has been
defined and compute interactions between frozen and non-frozen
particles appropriately, as if the frozen particle has infinite mass.
A similar functionality for normal (point) particles can be obtained
using :doc:`fix setforce <fix_setforce>`.

----------

.. include:: accel_styles.rst

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix computes a global 3-vector of forces, which can be accessed
by various :doc:`output commands <Howto_output>`.  This is the total
force on the group of atoms before the forces on individual atoms are
changed by the fix.  The vector values calculated by this fix are
"extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the GRANULAR package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

There can only be a single freeze fix defined.  This is because other
the :doc:`granular pair styles <pair_gran>` treat frozen particles
differently and need to be able to reference a single group to which
this fix is applied.

Related commands
""""""""""""""""

:doc:`atom_style sphere <atom_style>`, :doc:`fix setforce <fix_setforce>`

Default
"""""""

none
