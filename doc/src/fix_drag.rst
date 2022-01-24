.. index:: fix drag

fix drag command
================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID drag x y z fmag delta

* ID, group-ID are documented in :doc:`fix <fix>` command
* drag = style name of this fix command
* x,y,z = coord to drag atoms towards
* fmag = magnitude of force to apply to each atom (force units)
* delta = cutoff distance inside of which force         is not applied (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   fix center small-molecule drag 0.0 10.0 0.0 5.0 2.0

Description
"""""""""""

Apply a force to each atom in a group to drag it towards the point
(x,y,z).  The magnitude of the force is specified by fmag.  If an atom
is closer than a distance delta to the point, then the force is not
applied.

Any of the x,y,z values can be specified as NULL which means do not
include that dimension in the distance calculation or force
application.

This command can be used to steer one or more atoms to a new location
in the simulation.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is adding its forces. Default is the outermost level.

This fix computes a global 3-vector of forces, which can be accessed
by various :doc:`output commands <Howto_output>`.  This is the total
force on the group of atoms by the drag force.  The vector values
calculated by this fix are "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the EXTRA-FIX package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix spring <fix_spring>`, :doc:`fix spring/self <fix_spring_self>`,
:doc:`fix spring/rg <fix_spring_rg>`, :doc:`fix smd <fix_smd>`

Default
"""""""

none
