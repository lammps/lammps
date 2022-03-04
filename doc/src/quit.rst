.. index:: quit

quit command
============

Syntax
""""""

.. parsed-literal::

   quit status

status = numerical exit status (optional)

Examples
""""""""

.. code-block:: LAMMPS

   quit
   if "$n > 10000" then "quit 1"

Description
"""""""""""

This command causes LAMMPS to exit, after shutting down all output
cleanly.

It can be used as a debug statement in an input script, to terminate
the script at some intermediate point.

It can also be used as an invoked command inside the "then" or "else"
portion of an :doc:`if <if>` command.

The optional status argument is an integer which signals the return
status to a program calling LAMMPS.  A return status of 0 usually
indicates success.  A status != 0 is failure, where the specified
value can be used to distinguish the kind of error, e.g. where in the
input script the quit was invoked.  If not specified, a status of 0 is
returned.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`if <if>`

Default
"""""""

none
