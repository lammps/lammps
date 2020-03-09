.. index:: compute improper

compute improper command
========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID improper

* ID, group-ID are documented in :doc:`compute <compute>` command
* improper = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all improper

Description
"""""""""""

Define a computation that extracts the improper energy calculated by
each of the improper sub-styles used in the :doc:`improper_style hybrid <improper_hybrid>` command.  These values are made
accessible for output or further processing by other commands.  The
group specified for this command is ignored.

This compute is useful when using :doc:`improper_style hybrid <improper_hybrid>` if you want to know the portion of the
total energy contributed by one or more of the hybrid sub-styles.

**Output info:**

This compute calculates a global vector of length N where N is the
number of sub\_styles defined by the :doc:`improper_style hybrid <improper_style>` command.  which can be accessed by indices
1-N.  These values can be used by any command that uses global scalar
or vector values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

The vector values are "extensive" and will be in energy
:doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute pe <compute_pe>`, :doc:`compute pair <compute_pair>`

**Default:** none
