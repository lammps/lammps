.. index:: fix neb/spin

fix neb/spin command
====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID neb/spin Kspring

* ID, group-ID are documented in :doc:`fix <fix>` command
* neb/spin = style name of this fix command

.. parsed-literal::

   Kspring = spring constant for parallel nudging force
   (force/distance units or force units, see parallel keyword)

Examples
""""""""

fix 1 active neb/spin 1.0

Description
"""""""""""

Add nudging forces to spins in the group for a multi-replica
simulation run via the :doc:`neb/spin <neb_spin>` command to perform a 
geodesic nudged elastic band (GNEB) calculation for finding the 
transition state.
Hi-level explanations of GNEB are given with the 
:doc:`neb/spin <neb_spin>` command and on the 
:doc:`Howto replica <Howto_replica>` doc page.  
The fix neb/spin command must be used with the "neb/spin" command and 
defines how inter-replica nudging forces are computed.  A GNEB 
calculation is divided in two stages. In the first stage n replicas 
are relaxed toward a MEP until convergence.  In the second stage, the 
climbing image scheme is enabled, so that the replica having the highest 
energy relaxes toward the saddle point (i.e. the point of highest energy 
along the MEP), and a second relaxation is performed.

The nudging forces are calculated as explained in
:ref:`(BessarabB) <BessarabB>`).
See this reference for more explanation about their expression.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix\_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
as invoked by the :doc:`minimize <minimize>` command via the
:doc:`neb/spin <neb_spin>` command.

Restrictions
""""""""""""


This command can only be used if LAMMPS was built with the SPIN
package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`neb\_spin <neb_spin>`

Default
"""""""

none


----------


.. _BessarabB:



**(BessarabB)** Bessarab, Uzdin, Jonsson, Comp Phys Comm, 196,
335-347 (2015).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
