.. index:: fix heat/flow

fix heat/flow command
==========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID heat/flow style values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* heat/flow = style name of this fix command
* one style with corresponding value(s) needs to be listed

  .. parsed-literal::

     style = *constant* or *type*
       *constant* = cp
         cp = value of specifc heat (energy/(mass * temperature) units)
       *type* = cp1 ... cpN
         cpN = value of specifc heat for type N (energy/(mass * temperature) units)

*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all heat/flow constant 1.0
   fix 1 all heat/flow type 1.0 0.5

Description
"""""""""""

Perform plain time integration to update temperature for atoms in the
group each timestep. The specific heat of atoms can be defined using either
the *constant* or *type* keywords. For style *constant*, the specific heat
is a constant value *cp* for all atoms. For style *type*, *N* different values
of the specific heat are defined, one for each of the *N* types of atoms.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.
None of the :doc:`fix_modify <fix_modify>` options are relevant to this fix.
No global or per-atom quantities are stored by this fix for access by various
:doc:`output commands <Howto_output>`. No parameter of this fix can be used
with the *start/stop* keywords of the :doc:`run <run>` command.  This fix is
not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix requires that atoms store temperature and heat flow
as defined by the :doc:`fix property/atom <fix_property_atom>` command.

Related commands
""""""""""""""""

:doc:`pair granular <pair_granular>`, :doc:`fix property/atom <fix_property_atom>`

Default
"""""""

none
