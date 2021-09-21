.. index:: bond_style bpm/spring

bond_style bpm/spring command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/spring keyword value attribute1 attribute2 ...

* optional keyword = *store/local*

  .. parsed-literal::

       *store/local* values = ID of associated fix store/local followed by one or more attributes

          *id1, id2* = IDs of 2 atoms in the bond
          *time* = the time the bond broke
          *x, y, z* = the center of mass position of the 2 atoms when the bond broke
          *x/ref, y/ref, z/ref* = the inintial center of mass position of the 2 atoms


Examples
""""""""

.. code-block:: LAMMPS

   bond_style bpm/spring
   bond_coeff 1 1.0 0.05 0.1

   bond_style bpm/spring myfix time id1 id2
   fix myfix all store/local 1000 3
   dump 1 all local 1000 dump.broken f_myfix[1] f_myfix[2] f_myfix[3]
   dump_modify 1 write_header no

Description
"""""""""""

The *bpm/spring* bond style computes forces and torques based
on deviations from the initial reference state of the two atoms.
The reference state is stored by each bond when it is first computed
in the setup of a run. Data is then preserved across run commands and
is written to :doc:`binary restart files <restart>` such that restarting
the system will not reset the reference state of a bond.

This bond style only applies central-body forces which conserve the translational
and rotational degrees of freedom of a bonded set of particles. The force
has a magnitude of

.. math::

   f = k (r - r_0) w

where :math:`k_r` is a stiffness, :math:`r` is the current distance and
:math:`r_0` is the initial distance between the two particles, and :math:`w`
is a smoothing factor.
Bonds will break at a strain of :math:`\epsilon_c`.
This is done by setting by setting its type to 0 such that forces are
no longer computed.
The smoothing factor :math:`w` is constructed such that forces smoothly
go to zero, avoiding discontinuities, as bonds approach the critical strain

.. math::

   w = 1.0 - \left( \frac{r - r_0}{r_0 \epsilon_c} \right^4 .

Finally, additional damping forces and torques are applied to the two
particles. A force is applied proportional to the difference in the
normal velocity of particles using a similar construction as
dissipative particle dynamics (:ref:`(Groot) <Groot1>`):

.. math::

   F_D = - \gamma_n w (\hat{r} \bullet \vec{v})

where :math:`\gamma_n` is the damping strength, :math:`\hat{r}` is the
radial normal vector, and :math:`\vec{v}` is the velocity difference
between the two particles.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`k`             (force/distance units)
* :math:`\epsilon_c`    (unit less)
* :math:`\gamma`        (force/velocity units)

As bonds can be broken between neighbor list builds, particular
:doc:`special_bonds <special_bonds>` are required. See the `:doc: how to <Howto_BPM>`
page on BPMs or `:doc: fix update/special/bonds <fix_update_special_bonds>`
for details.

This bond style tracks broken bonds and can record them using an instance of
:doc:`fix store/local <fix_store_local>` if the *store/local* keyword is
used followed by the ID of the fix and then a series of bond attributes.

Note that when bonds are dumped to a file via the :doc:`dump local <dump>`
command, bonds with type 0 (broken bonds) are not included.  The
:doc:`delete_bonds <delete_bonds>` command can also be used to query the
status of broken bonds or permanently delete them, e.g.:

.. code-block:: LAMMPS

   delete_bonds all stats
   delete_bonds all bond 0 remove


----------

Restart
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This bond style writes the reference state of each bond to
:doc:`binary restart files <restart>`. Loading a restart
file will properly resume bonds.

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the BPM
package. See the :doc:`Build package <Build_package>` doc page for more
info.

The *bpm/spring* style requires 1-3 and 1-4 :doc:`special_bonds <special_bonds>`
be turned off using the :doc:`special_bonds <special_bonds>` command.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`fix store/local <fix_store_local>`

Default
"""""""

none


.. _Groot1:

**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).
