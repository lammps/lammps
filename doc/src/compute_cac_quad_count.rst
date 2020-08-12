.. index:: compute cac/quad/count

compute cac/quad/count command
===============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID cac/quad/count

* ID, group-ID are documented in :doc:`compute <compute>` command
* cac/quad/count = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS
   compute Eweight all cac/quad/count

Description
"""""""""""

Define a computation that calculates the computational weight required
by each atom and finite element in CAC model. This is performed by using
the neighbor list sizes associated with the force computation.

.. note::

   The first time this compute is invoked it will simply count the 
   number of integration points and atoms owned.

----------

**Output info:**

This compute calculates a per atom scalar (computational weight) that is
associated with each CAC element and atom in the model.These values are intended
to be used by the variable command in order to provide them to the :doc:`fix balance command <fix_balance>`

Restrictions
""""""""""""
 requires a CAC atom style to be defined

Related commands
""""""""""""""""

:doc:`fix_balance <fix_balance>`, :doc:`variable <variable>`

**Default:** none
