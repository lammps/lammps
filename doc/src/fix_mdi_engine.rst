.. index:: fix mdi/engine

fix mdi/engine command
======================

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

This fix is used along with the :doc:`mdi/engine <mdi_engine>` command
to enable LAMMPS to use the `MDI Library
<https://molssi-mdi.github.io/MDI_Library/html/index.html>`_ to run as
an MDI engine.  The fix provides hooks that enable MDI driver codes to
communicate with LAMMPS at various points within a LAMMPS timestep.

It is not generally necessary to add this fix to a LAMMPS input file,
even when using the :doc:`mdi/engine <mdi_engine>` command.  If the
:doc:`mdi/engine <mdi_engine>` command is executed and this fix is not
present, it will automatically be added and applied as a new fix for
all atoms for the duration of the command.  Thus it is only necessary
to add this fix to an input file when you want to modify the group-ID
or the ordering of this fix relative to other fixes in the input script.

For more information about running LAMMPS as an MDI engine, see the
:doc:`mdi/engine <mdi_engine>` command and the :doc:`Howto mdi
<Howto_mdi>` doc page.

Restrictions
""""""""""""

This command is part of the MDI package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`mdi/engine <mdi_engine>`

Default
"""""""

none
