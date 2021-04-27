.. index:: fix move

fix mdi/engine command
================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID mdi/engine

* ID, group-ID are documented in :doc:`fix <fix>` command
* mdi/engine = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all mdi/engine

Description
"""""""""""

This fix is used alongside the :doc:`mdi_engine <mdi_engine>` command
to enable LAMMPS to use the
`MDI Library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`_
to run as an MDI engine.
The fix provides hooks that are enable external drivers to communicate
with LAMMPS at various points within LAMMPS timesteps.

It is not generally necessary to add this fix to a LAMMPS input file,
even when using the :doc:`mdi_engine <mdi_engine>` command; if the
:doc:`mdi_engine <mdi_engine>` command is executed and this fix is not
present, it will automatically be added to the end of the fix list and
applied to all atoms for the duration of the command.  It is only
necessary to add this fix to an input file for cases in which the user
would like to modify order or group-ID of the fix.

For more information about running LAMMPS as an MDI engine, see the
:doc:`mdi_engine <mdi_engine>` command.


Restrictions
""""""""""""
This command is part of the USER-MDI package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`mdi_engine <mdi_engine>`

Default
"""""""

none
