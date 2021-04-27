.. index:: mdi_engine

mdi_engine command
===============

Syntax
""""""

.. parsed-literal::

   mdi_engine

Description
"""""""""""

This command causes LAMMPS to use the
`MDI Library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`_
to run as an MDI Engine, responding to commands from an external
MDI Driver.
General information about launching codes that communicate using the
MDI Library can be found in the
`corresponding page <https://molssi-mdi.github.io/MDI_Library/html/library_page.html#library_launching_sec>`_
of the MDI Library's documentation.

----------

Other commands can be executed both before and after this command,
but typically this command should be used at the end of a LAMMPS
input script.

Restrictions
""""""""""""

This command is part of the USER-MDI package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix mdi/engine <fix_mdi_engine>`

Default
"""""""

none
