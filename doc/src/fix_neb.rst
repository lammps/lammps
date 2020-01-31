.. index:: fix neb

fix neb command
===============

Syntax
""""""


.. parsed-literal::

   fix ID group-ID neb Kspring keyword value

* ID, group-ID are documented in :doc:`fix <fix>` command
* neb = style name of this fix command
* Kspring = spring constant for parallel nudging force (force/distance units or force units, see parallel keyword)
* zero or more keyword/value pairs may be appended
* keyword = *parallel* or *perp* or *end*

.. parsed-literal::

     *parallel* value = *neigh* or *ideal*
       *neigh* = parallel nudging force based on distance to neighbor replicas (Kspring = force/distance units)
       *ideal* = parallel nudging force based on interpolated ideal position (Kspring = force units)
     *perp* value = *Kspring2*
       *Kspring2* = spring constant for perpendicular nudging force (force/distance units)
     *end* values = estyle Kspring3
       *estyle* = *first* or *last* or *last/efirst* or *last/efirst/middle*
         *first* = apply force to first replica
         *last* = apply force to last replica
         *last/efirst* = apply force to last replica and set its target energy to that of first replica
         *last/efirst/middle* = same as *last/efirst* plus prevent middle replicas having lower energy than first replica
       *Kspring3* = spring constant for target energy term (1/distance units)

Examples
""""""""


.. parsed-literal::

   fix 1 active neb 10.0
   fix 2 all neb 1.0 perp 1.0 end last
   fix 2 all neb 1.0 perp 1.0 end first 1.0 end last 1.0
   fix 1 all neb 1.0 parallel ideal end last/efirst 1

Description
"""""""""""

Add nudging forces to atoms in the group for a multi-replica
simulation run via the :doc:`neb <neb>` command to perform a nudged
elastic band (NEB) calculation for finding the transition state.
Hi-level explanations of NEB are given with the :doc:`neb <neb>` command
and on the :doc:`Howto replica <Howto_replica>` doc page.  The fix neb
command must be used with the "neb" command and defines how
inter-replica nudging forces are computed.  A NEB calculation is
divided in two stages. In the first stage n replicas are relaxed
toward a MEP until convergence.  In the second stage, the climbing
image scheme (see :ref:`(Henkelman2) <Henkelman2>`) is enabled, so that the
replica having the highest energy relaxes toward the saddle point
(i.e. the point of highest energy along the MEP), and a second
relaxation is performed.

A key purpose of the nudging forces is to keep the replicas equally
spaced.  During the NEB calculation, the 3N-length vector of
interatomic force Fi = -Grad(V) for each replica I is altered.  For
all intermediate replicas (i.e. for 1 < I < N, except the climbing
replica) the force vector becomes:


.. parsed-literal::

   Fi = -Grad(V) + (Grad(V) dot T') T' + Fnudge_parallel + Fnudge_perp

T' is the unit "tangent" vector for replica I and is a function of Ri,
Ri-1, Ri+1, and the potential energy of the 3 replicas; it points
roughly in the direction of (Ri+i - Ri-1); see the
:ref:`(Henkelman1) <Henkelman1>` paper for details.  Ri are the atomic
coordinates of replica I; Ri-1 and Ri+1 are the coordinates of its
neighbor replicas.  The term (Grad(V) dot T') is used to remove the
component of the gradient parallel to the path which would tend to
distribute the replica unevenly along the path.  Fnudge\_parallel is an
artificial nudging force which is applied only in the tangent
direction and which maintains the equal spacing between replicas (see
below for more information).  Fnudge\_perp is an optional artificial
spring which is applied in a direction perpendicular to the tangent
direction and which prevent the paths from forming acute kinks (see
below for more information).

In the second stage of the NEB calculation, the interatomic force Fi
for the climbing replica (the replica of highest energy after the
first stage) is changed to:


.. parsed-literal::

   Fi = -Grad(V) + 2 (Grad(V) dot T') T'

and the relaxation procedure is continued to a new converged MEP.


----------


The keyword *parallel* specifies how the parallel nudging force is
computed.  With a value of *neigh*\ , the parallel nudging force is
computed as in :ref:`(Henkelman1) <Henkelman1>` by connecting each
intermediate replica with the previous and the next image:


.. parsed-literal::

   Fnudge_parallel = *Kspring* \* (\|Ri+1 - Ri\| - \|Ri - Ri-1\|)

Note that in this case the specified *Kspring* is in force/distance
units.

With a value of *ideal*\ , the spring force is computed as suggested in
ref`(WeinanE) <WeinanE>` 


.. parsed-literal::

   Fnudge_parallel = -\ *Kspring* \* (RD-RDideal) / (2 \* meanDist)

where RD is the "reaction coordinate" see :doc:`neb <neb>` section, and
RDideal is the ideal RD for which all the images are equally spaced.
I.e. RDideal = (I-1)\*meanDist when the climbing replica is off, where
I is the replica number).  The meanDist is the average distance
between replicas.  Note that in this case the specified *Kspring* is
in force units.

Note that the *ideal* form of nudging can often be more effective at
keeping the replicas equally spaced.


----------


The keyword *perp* specifies if and how a perpendicular nudging force
is computed.  It adds a spring force perpendicular to the path in
order to prevent the path from becoming too strongly kinked.  It can
significantly improve the convergence of the NEB calculation when the
resolution is poor.  I.e. when few replicas are used; see
:ref:`(Maras) <Maras1>` for details.

The perpendicular spring force is given by


.. parsed-literal::

   Fnudge_perp = *Kspring2* \* F(Ri-1,Ri,Ri+1) (Ri+1 + Ri-1 - 2 Ri)

where *Kspring2* is the specified value.  F(Ri-1 Ri R+1) is a smooth
scalar function of the angle Ri-1 Ri Ri+1.  It is equal to 0.0 when
the path is straight and is equal to 1 when the angle Ri-1 Ri Ri+1 is
acute.  F(Ri-1 Ri R+1) is defined in :ref:`(Jonsson) <Jonsson>`.

If *Kspring2* is set to 0.0 (the default) then no perpendicular spring
force is added.


----------


By default, no additional forces act on the first and last replicas
during the NEB relaxation, so these replicas simply relax toward their
respective local minima.  By using the key word *end*\ , additional
forces can be applied to the first and/or last replicas, to enable
them to relax toward a MEP while constraining their energy E to the
target energy ETarget.

If ETarget>E, the interatomic force Fi for the specified replica becomes:


.. parsed-literal::

   Fi = -Grad(V) + (Grad(V) dot T' + (E-ETarget)\*Kspring3) T',  *when* Grad(V) dot T' < 0
   Fi = -Grad(V) + (Grad(V) dot T' + (ETarget- E)\*Kspring3) T', *when* Grad(V) dot T' > 0

The "spring" constant on the difference in energies is the specified
*Kspring3* value.

When *estyle* is specified as *first*\ , the force is applied to the
first replica.  When *estyle* is specified as *last*\ , the force is
applied to the last replica.  Note that the *end* keyword can be used
twice to add forces to both the first and last replicas.

For both these *estyle* settings, the target energy *ETarget* is set
to the initial energy of the replica (at the start of the NEB
calculation).

If the *estyle* is specified as *last/efirst* or *last/efirst/middle*\ ,
force is applied to the last replica, but the target energy *ETarget*
is continuously set to the energy of the first replica, as it evolves
during the NEB relaxation.

The difference between these two *estyle* options is as follows.  When
*estyle* is specified as *last/efirst*\ , no change is made to the
inter-replica force applied to the intermediate replicas (neither
first or last).  If the initial path is too far from the MEP, an
intermediate replica may relax "faster" and reach a lower energy than
the last replica.  In this case the intermediate replica will be
relaxing toward its own local minima.  This behavior can be prevented
by specifying *estyle* as *last/efirst/middle* which will alter the
inter-replica force applied to intermediate replicas by removing the
contribution of the gradient to the inter-replica force.  This will
only be done if a particular intermediate replica has a lower energy
than the first replica.  This should effectively prevent the
intermediate replicas from over-relaxing.

After converging a NEB calculation using an *estyle* of
*last/efirst/middle*\ , you should check that all intermediate replicas
have a larger energy than the first replica. If this is not the case,
the path is probably not a MEP.

Finally, note that the last replica may never reach the target energy
if it is stuck in a local minima which has a larger energy than the
target energy.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
as invoked by the :doc:`minimize <minimize>` command via the
:doc:`neb <neb>` command.

Restrictions
""""""""""""


This command can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`neb <neb>`

Default
"""""""

The option defaults are parallel = neigh, perp = 0.0, ends is not
specified (no inter-replica force on the end replicas).


----------


.. _Henkelman1:



**(Henkelman1)** Henkelman and Jonsson, J Chem Phys, 113, 9978-9985 (2000).

.. _Henkelman2:



**(Henkelman2)** Henkelman, Uberuaga, Jonsson, J Chem Phys, 113,
9901-9904 (2000).

.. _WeinanE:



**(WeinanE)** E, Ren, Vanden-Eijnden, Phys Rev B, 66, 052301 (2002).

.. _Jonsson:



**(Jonsson)** Jonsson, Mills and Jacobsen, in Classical and Quantum
Dynamics in Condensed Phase Simulations, edited by Berne, Ciccotti,
and Coker World Scientific, Singapore, 1998, p 385.

.. _Maras1:



**(Maras)** Maras, Trushin, Stukowski, Ala-Nissila, Jonsson,
Comp Phys Comm, 205, 13-21 (2016).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
