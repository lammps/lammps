.. index:: bond_style bpm/rotational

bond_style bpm/rotational command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/rotational keyword value attribute1 attribute2 ...

* optional keyword = *overlay/pair* or *store/local*

  .. parsed-literal::

       *store/local* values = ID of associated fix store/local followed by one or more attributes

          *id1, id2* = IDs of 2 atoms in the bond
          *time* = the timestep the bond broke
          *x, y, z* = the center of mass position of the 2 atoms when the bond broke (distance units)
          *x/ref, y/ref, z/ref* = the initial center of mass position of the 2 atoms (distance units)

       *overlay/pair* value = none
          bonded particles will still interact with pair forces

Examples
""""""""

.. code-block:: LAMMPS

   bond_style bpm/rotational
   bond_coeff 1 1.0 0.2 0.02 0.02 0.20 0.04 0.04 0.04 0.1 0.02 0.002 0.002

   bond_style bpm/rotational myfix time id1 id2
   fix myfix all store/local 1000 3
   dump 1 all local 1000 dump.broken f_myfix[1] f_myfix[2] f_myfix[3]
   dump_modify 1 write_header no

Description
"""""""""""

The *bpm/rotational* bond style computes forces and torques based on
deviations from the initial reference state of the two atoms.  The
reference state is stored by each bond when it is first computed in
the setup of a run. Data is then preserved across run commands and is
written to :doc:`binary restart files <restart>` such that restarting
the system will not reset the reference state of a bond.

Forces include a normal and tangential component. The base normal force
has a magnitude of

.. math::

   f_r = k_r (r - r_0)

where :math:`k_r` is a stiffness and :math:`r` is the current distance and
:math:`r_0` is the initial distance between the two particles.

A tangential force is applied perpendicular to the normal direction
which is proportional to the tangential shear displacement with a
stiffness of :math:`k_s`. This tangential force also induces a torque.
In addition, bending and twisting torques are also applied to
particles which are proportional to angular bending and twisting
displacements with stiffnesses of :math`k_b` and :math:`k_t',
respectively.  Details on the calculations of shear displacements and
angular displacements can be found in :ref:`(Wang) <Wang2009>` and
:ref:`(Wang and Mora) <WangMora2009b>`.

Bonds will break under sufficient stress. A breaking criteria is calculated

.. math::

   B = \mathrm{max}\{0, \frac{f_r}{f_{r,c}} + \frac{|f_s|}{f_{s,c}} +
       \frac{|\tau_b|}{\tau_{b,c}} + \frac{|\tau_t|}{\tau_{t,c}} \}

where :math:`|f_s|` is the magnitude of the shear force and
:math:`|\tau_b|` and :math:`|\tau_t|` are the magnitudes of the
bending and twisting forces, respectively. The corresponding variables
:math:`f_{r,c}` :math:`f_{s,c}`, :math:`\tau_{b,c}`, and
:math:`\tau_{t,c}` are critical limits to each force or torque.  If
:math:`B` is ever equal to or exceeds one, the bond will break.  This
is done by setting by setting its type to 0 such that forces and
torques are no longer computed.

After computing the base magnitudes of the forces and torques, they
are all multiplied by an extra factor :math:`w` to smoothly
interpolate forces and torques to zero as the bond breaks. This term
is calculated as :math:`w = (1.0 - B^4)`.

Finally, additional damping forces and torques are applied to the two
particles. A force is applied proportional to the difference in the
normal velocity of particles using a similar construction as
dissipative particle dynamics (:ref:`(Groot) <Groot1>`):

.. math::

   F_D = - \gamma_n w (\hat{r} \bullet \vec{v})

where :math:`\gamma_n` is the damping strength, :math:`\hat{r}` is the
radial normal vector, and :math:`\vec{v}` is the velocity difference
between the two particles. Similarly, tangential forces are applied to
each atom proportional to the relative differences in sliding
velocities with a constant prefactor :math:`\gamma_s` (:ref:`(Wang et
al.) <Wang2015>) along with their associated torques. The rolling and
twisting components of the relative angular velocities of the two
atoms are also damped by applying torques with prefactors of
:math:`\gamma_r` and :math:`\gamma_t`, respectively.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`k_r`           (force/distance units)
* :math:`k_s`           (force/distance units)
* :math:`k_t`           (force units)
* :math:`k_b`           (force units)
* :math:`f_{r,c}`       (force units)
* :math:`f_{s,c}`       (force units)
* :math:`\tau_{b,c}`    (force*distance units)
* :math:`\tau_{t,c}`    (force*distance units)
* :math:`\gamma_n`      (force/velocity units)
* :math:`\gamma_s`      (force/velocity units)
* :math:`\gamma_r`      (distance*force/seconds/radians units)
* :math:`\gamma_t`      (distance*force/seconds/radians units)

By default, pair forces are not calculated between bonded particles.
Pair forces can alternatively be overlaid on top of bond forces using
the *overlay/pair* keyword. These settings require specific
:doc:`special_bonds <special_bonds>` settings described in the
restrictions.  Further details can be found in the `:doc: how to
<Howto_BPM>` page on BPMs.

This bond style tracks broken bonds and can record them using an
instance of :doc:`fix store/local <fix_store_local>` if the
*store/local* keyword is used followed by the ID of the fix and then a
series of bond attributes.

Note that when bonds are dumped to a file via the :doc:`dump local <dump>`
command, bonds with type 0 (broken bonds) are not included.  The
:doc:`delete_bonds <delete_bonds>` command can also be used to query the
status of broken bonds or permanently delete them, e.g.:

.. code-block:: LAMMPS

   delete_bonds all stats
   delete_bonds all bond 0 remove


----------

Restart and other info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This bond style writes the reference state of each bond to
:doc:`binary restart files <restart>`. Loading a restart file will
properly resume bonds.

The single() function of these pair styles returns 0.0 for the energy
of a pairwise interaction, since energy is not conserved in these
dissipative potentials.  It also returns only the normal component of
the pairwise interaction force.

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the BPM
package. See the :doc:`Build package <Build_package>` doc page for
more info.

By default if pair interactions are to be disabled, this bond style
requires setting

.. code-block:: LAMMPS

   special_bonds lj 0 1 1 coul 1 1 1

and :doc:`newton <newton>` must be set to bond off.  If the
*overlay/pair* option is used, this bond style alternatively requires
setting

.. code-block:: LAMMPS

   special_bonds lj/coul 1 1 1

The *bpm/rotational* style requires :doc:`atom style sphere/bpm <atom_style>`.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`fix store/local <fix_store_local>`,
:doc:`fix nve/sphere/bpm <fix_nve_sphere_bpm>`

Default
"""""""

none

----------

.. _Wang2009:

**(Wang)** Wang, Acta Geotechnica, 4,
p 117-127 (2009).

.. _Wang2009b:

**(Wang and Mora)** Wang, Mora, Advances in Geocomputing,
119, p 183-228 (2009).

.. _Groot1:

**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).

.. _Wang2015:

**(Wang et al, 2015)** Wang, Y., Alonso-Marroquin, F., & Guo,
W. W. (2015).  Rolling and sliding in 3-D discrete element
models. Particuology, 23, 49-55.
