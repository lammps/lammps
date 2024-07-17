.. index:: fix surface/global

fix surface/global command
===============

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID surface/global

* ID, group-ID are documented in :doc:`fix <fix>` command
* surface/global = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all surface/global

Description
"""""""""""

Mention Howto granular surfaces page.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.  No global or per-atom quantities are stored by
this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`fix surface/local <fix_surface_local>`

Default
"""""""

none
