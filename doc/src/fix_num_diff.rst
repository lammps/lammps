.. index:: fix numdiff

fix numdiff command
====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID numdiff Nevery Delta

* ID, group-ID are documented in :doc:`fix <fix>` command
* numdiff = style name of this fix command
* Nevery = calculate force by finite difference every this many timesteps
* Delta = finite difference displacement length (distance units)
  

Examples
""""""""


.. parsed-literal::

   fix 1 all numdiff 1 0.0001
   fix 1 all numdiff 10 1e-6
   fix 1 all numdiff 100 0.01

Description
"""""""""""

Calculate forces through finite difference of energy versus position.
The resulting per-atom forces can be used by :doc:`dump custom <dump>`.

The group specified with the command means only atoms within the group
have their averages computed.  Results are set to 0.0 for atoms not in
the group.


----------


The *Nevery*\ argument specifies on what timesteps the force will 
be used calculated by finite difference.

The *Delta*\ argument specifies the positional displacement each
atom will undergo.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global scalar or vector quantities are
stored by this fix for access by various :doc:`output commands <Howto_output>`.

This fix produces a per-atom array which can be accessed by
various :doc:`output commands <Howto_output>`.  .  This fix produces
a per-atom array of the forces calculated by finite difference. The
per-atom values should only be accessed on timesteps that are multiples
of *Nfreq* since that is when the finite difference forces are calculated.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`dynamical_matrix <dynamical_matrix>`,

**Default:** none
