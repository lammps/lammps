.. index:: fix align/neighbor

fix align/neighbor command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID align/neighbor mode magnitude cutoff keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* align/neighbor = style name of this fix command
* mode = *dipole*
* magnitude = prefactor of the alignment torque (energy units)
* cutoff = distance within which particles align (distance units)
* zero or more keyword/value pairs may be appended
* keyword = *symmetry*

  .. parsed-literal::

       *symmetry* value = *polar* or *nematic*
         *polar* = align orientations in parallel
         *nematic* = align orientations in parallel or antiparallel, whichever is closer

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all align/neighbor dipole 2.0 1.5
   fix 1 rods align/neighbor dipole 0.5 1.2 symmetry nematic

Description
"""""""""""

.. versionadded:: TBD

Add a torque to each atom in the group that rotates its orientation
toward the orientations of the other atoms in the group within the
*cutoff* distance.  This is the local alignment interaction of models of
collective motion ("flocking") like the Vicsek model :ref:`(Vicsek1995)
<Vicsek1995>` for polar alignment, and of models of self-propelled rods
and active nematics :ref:`(Ginelli2010) <Ginelli2010>` for nematic
alignment, formulated as a continuous-time torque :ref:`(Peruani2008)
<Peruani2008>` so that it can be combined with the time integrators of
LAMMPS.  See :ref:`(Chate2020) <Chate2020>` for a review of such models
and the :doc:`Howto active matter <Howto_active>` page for an overview
of active matter models in LAMMPS.

The torque on atom *i* is

.. math::

   \mathbf{\tau}_i = K \sum_{j} (\mathbf{e}_i \times \mathbf{e}_j)
   \qquad \mathrm{(polar)}

.. math::

   \mathbf{\tau}_i = K \sum_{j} (\mathbf{e}_i \cdot \mathbf{e}_j)\,(\mathbf{e}_i \times \mathbf{e}_j)
   \qquad \mathrm{(nematic)}

where :math:`K` is the *magnitude*, :math:`\mathbf{e}_i` is the
orientation of atom *i*, and the sum runs over all atoms *j* of the
group within the cutoff distance of atom *i*.  For unit orientation
vectors with an angle :math:`\theta_{ij}` between them, the polar torque
has the magnitude :math:`K \sin\theta_{ij}` and rotates
:math:`\mathbf{e}_i` toward :math:`\mathbf{e}_j`; the nematic torque has
the magnitude :math:`\frac{K}{2} \sin 2\theta_{ij}` and rotates
:math:`\mathbf{e}_i` toward :math:`\mathbf{e}_j` or
:math:`-\mathbf{e}_j`, whichever is closer.  The torques of a pair of
atoms are equal and opposite, and equivalent to those resulting from the
pair energies :math:`-K\, \mathbf{e}_i \cdot \mathbf{e}_j` (polar) and
:math:`-\frac{K}{2} (\mathbf{e}_i \cdot \mathbf{e}_j)^2` (nematic, as in
the Lebwohl-Lasher model of liquid crystals :ref:`(Lebwohl1972)
<Lebwohl1972>`).  These energies are not computed.  A negative
*magnitude* results in anti-alignment.  Since the torques are summed
over the neighbors rather than averaged, the tendency to align grows
with the local density.

The torque only changes the orientations when it is applied by a time
integrator that rotates the orientation vectors.  With the overdamped
rotational dynamics of :doc:`fix brownian/sphere <fix_brownian>` and its
rotational friction coefficient :math:`\gamma_r`, the orientations obey
:math:`d\mathbf{e}_i/dt = (\mathbf{\tau}_i/\gamma_r) \times
\mathbf{e}_i` plus rotational noise, which for polar alignment in 2d
becomes

.. math::

   \frac{d\theta_i}{dt} = \frac{K}{\gamma_r} \sum_j \sin(\theta_j - \theta_i) + \mathrm{noise}

with :math:`\theta_i` the angle of the orientation of atom *i*.  The
alignment rate :math:`K/\gamma_r` and the rotational diffusion
coefficient (the noise strength) then play the roles of the alignment
strength and the noise of the Vicsek model, respectively.  Since the
alignment enters as a rate, the resulting dynamics do not depend on the
size of the timestep, as long as it is small compared to the inverse of
the alignment rate.  With inertial rotational dynamics, as in :doc:`fix
nve/sphere <fix_nve_sphere>` with the *update dipole* keyword, the
torque changes the angular velocity of the particles, and rotational
damping, e.g. from :doc:`fix langevin <fix_langevin>` with the *omega*
keyword, is required for the particles to actually align.

For mode *dipole*, the orientation :math:`\mathbf{e}_i` of a particle is
its dipole vector :math:`\mathbf{\mu}_i` as defined by :doc:`atom_style
dipole <atom_style>`, and the torque is added to the per-atom torque
defined by :doc:`atom_style sphere <atom_style>`.  Both are available
with :doc:`atom_style hybrid sphere dipole <atom_style>`.  If the dipole
vectors are not unit vectors, the torque scales with the product of
their lengths.  The dipole moments can be set with the :doc:`set <set>`
command, e.g. to random unit vectors with the *dipole/random* keyword.

This fix builds its own neighbor list with the given *cutoff* for all
pairs of atom types, which is independent of the cutoff of the pair
style.  The cutoff plus the neighbor list skin must not exceed the
communication cutoff, which by default is determined by the largest pair
style cutoff.  If the alignment cutoff is larger than the pair style
cutoff (or there is no pair style), the communication cutoff has to be
increased accordingly with the :doc:`comm_modify cutoff <comm_modify>`
command, otherwise LAMMPS stops with an error.

Only atoms in the fix group receive alignment torques, and only atoms in
the fix group are included in the sums.  Atoms outside the group neither
exert nor receive alignment torques.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix.  This allows to set at which level of the :doc:`r-RESPA
<run_style>` integrator the fix is adding its torques.  Default is the
outermost level.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the BROWNIAN package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Mode *dipole* requires that atoms store torque as defined by the
:doc:`atom_style sphere <atom_style>` command, as well as a dipole
moment as defined by the :doc:`atom_style dipole <atom_style>` command,
which is part of the DIPOLE package.

Related commands
""""""""""""""""

:doc:`fix align/self <fix_align_self>`, :doc:`fix propel/self <fix_propel_self>`,
:doc:`fix brownian/sphere <fix_brownian>`, :doc:`fix nve/sphere <fix_nve_sphere>`,
:doc:`fix tumble <fix_tumble>`, :doc:`comm_modify <comm_modify>`

Default
"""""""

The default is *symmetry* = *polar*.

----------

.. _Vicsek1995:

**(Vicsek1995)** T. Vicsek, A. Czirok, E. Ben-Jacob, I. Cohen, and O. Shochet, Novel Type of Phase Transition in a System of Self-Driven Particles, Phys. Rev. Lett. 75, 1226 (1995).

.. _Ginelli2010:

**(Ginelli2010)** F. Ginelli, F. Peruani, M. Baer, and H. Chate, Large-Scale Collective Properties of Self-Propelled Rods, Phys. Rev. Lett. 104, 184502 (2010).

.. _Peruani2008:

**(Peruani2008)** F. Peruani, A. Deutsch, and M. Baer, A mean-field theory for self-propelled particles interacting by velocity alignment mechanisms, Eur. Phys. J. Special Topics 157, 111 (2008).

.. _Chate2020:

**(Chate2020)** H. Chate, Dry Aligning Dilute Active Matter, Annu. Rev. Condens. Matter Phys. 11, 189 (2020).

.. _Lebwohl1972:

**(Lebwohl1972)** P. A. Lebwohl and G. Lasher, Nematic-Liquid-Crystal Order - A Monte Carlo Calculation, Phys. Rev. A 6, 426 (1972).
