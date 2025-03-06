.. index:: fix manifoldforce

fix manifoldforce command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID manifoldforce manifold manifold-args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* manifold = name of the manifold
* manifold-args = parameters for the manifold

Examples
""""""""

.. code-block:: LAMMPS

   fix constrain all manifoldforce sphere 5.0

Description
"""""""""""

This fix subtracts each time step from the force the component along
the normal of the specified :doc:`manifold <Howto_manifold>`.  This can be
used in combination with :doc:`minimize <minimize>` to remove overlap
between particles while keeping them (roughly) constrained to the
given manifold, e.g. to set up a run with :doc:`fix nve/manifold/rattle <fix_nve_manifold_rattle>`.  I have found that
only *hftn* and *quickmin* with a very small time step perform
adequately though.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is invoked during :doc:`energy minimization <minimize>`.

----------

Restrictions
""""""""""""

This fix is part of the MANIFOLD package. It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Only use this with *min_style hftn* or *min_style quickmin*. If not,
the constraints will not be satisfied very well at all. A warning is
generated if the *min_style* is incompatible but no error.

----------

Related commands
""""""""""""""""

:doc:`fix nve/manifold/rattle <fix_nve_manifold_rattle>`, :doc:`fix nvt/manifold/rattle <fix_nvt_manifold_rattle>`
