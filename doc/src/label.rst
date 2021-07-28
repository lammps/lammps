.. index:: label

label command
=============

Syntax
""""""

.. parsed-literal::

   label ID

* ID = string used as label name

Examples
""""""""

.. code-block:: LAMMPS

   label xyz
   label loop

Description
"""""""""""

Label this line of the input script with the chosen ID.  Unless a jump
command was used previously, this does nothing.  But if a
:doc:`jump <jump>` command was used with a label argument to begin
invoking this script file, then all command lines in the script prior
to this line will be ignored.  I.e. execution of the script will begin
at this line.  This is useful for looping over a section of the input
script as discussed in the :doc:`jump <jump>` command.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

none


Default
"""""""

none
