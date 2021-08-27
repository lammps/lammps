.. index:: fix update/special/bonds

fix update/special/bonds command
======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID update/special/bonds

* ID, group-ID are documented in :doc:`fix <fix>` command
* update/special/bonds = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all update/special/bonds
   
Description
"""""""""""

This fix is used to update the 1-2 special bond list for BPM bond styles. 
This feature is used to censor pair forces between bonded particles.
See the :doc:`BPM how to <Howto_bpm>` for more information.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix requires :doc:`newton <newton>` bond off. 

Related commands
""""""""""""""""

:doc:`bond bpm/rotational <bond_bpm_rotational>`

Default
"""""""

none

