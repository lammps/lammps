.. index:: echo

echo command
============

Syntax
""""""


.. parsed-literal::

   echo style

* style = *none* or *screen* or *log* or *both*

Examples
""""""""


.. parsed-literal::

   echo both
   echo log

Description
"""""""""""

This command determines whether LAMMPS echoes each input script
command to the screen and/or log file as it is read and processed.  If
an input script has errors, it can be useful to look at echoed output
to see the last command processed.

The :doc:`command-line switch <Run_options>` -echo can be used in place
of this command.

Restrictions
""""""""""""
 none

**Related commands:** none

Default
"""""""


.. parsed-literal::

   echo log


