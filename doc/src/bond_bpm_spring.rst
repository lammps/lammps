.. index:: bond_style bpm/spring

bond_style bpm/spring command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/spring keyword value attribute1 attribute2 ...

* optional keyword = *overlay/pair* or *store/local* or *smooth*

  .. parsed-literal::

       *store/local* values = fix_ID N attributes ...
          * fix_ID = ID of associated internal fix to store data
          * N = prepare data for output every this many timesteps
          * attributes = zero or more of the below attributes may be appended

            *id1, id2* = IDs of 2 atoms in the bond
            *time* = the timestep the bond broke
            *x, y, z* = the center of mass position of the 2 atoms when the bond broke (distance units)
            *x/ref, y/ref, z/ref* = the initial center of mass position of the 2 atoms (distance units)

       *overlay/pair* value = none
          bonded particles will still interact with pair forces

       *smooth* value = *yes* or *no*
          smooths bond forces near the breaking point

Examples
""""""""

.. code-block:: LAMMPS

   bond_style bpm/spring
   bond_coeff 1 1.0 0.05 0.1

   bond_style bpm/spring myfix 1000 time id1 id2
   dump 1 all local 1000 dump.broken f_myfix[1] f_myfix[2] f_myfix[3]
   dump_modify 1 write_header no

Description
"""""""""""

The *bpm/spring* bond style computes forces and torques based on
deviations from the initial reference state of the two atoms.  The
reference state is stored by each bond when it is first computed in
the setup of a run. Data is then preserved across run commands and is
written to :doc:`binary restart files <restart>` such that restarting
the system will not reset the reference state of a bond.

This bond style only applies central-body forces which conserve the
translational and rotational degrees of freedom of a bonded set of
particles. The force has a magnitude of

.. math::

   F = k (r - r_0) w

where :math:`k_r` is a stiffness, :math:`r` is the current distance
and :math:`r_0` is the initial distance between the two particles, and
:math:`w` is an optional smoothing factor discussed below. Bonds will
break at a strain of :math:`\epsilon_c`.  This is done by setting by
setting its type to 0 such that forces are no longer computed.

An additional damping force is applied to the bonded
particles.  This forces is proportional to the difference in the
normal velocity of particles using a similar construction as
dissipative particle dynamics (:ref:`(Groot) <Groot4>`):

.. math::

   F_D = - \gamma w (\hat{r} \bullet \vec{v})

where :math:`\gamma` is the damping strength, :math:`\hat{r}` is the
radial normal vector, and :math:`\vec{v}` is the velocity difference
between the two particles.

The smoothing factor :math:`w` can be added or removed using the
*smooth* keyword. It is constructed such that forces smoothly go
to zero, avoiding discontinuities, as bonds approach the critical strain

.. math::

   w = 1.0 - \left( \frac{r - r_0}{r_0 \epsilon_c} \right)^8 .

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands:

* :math:`k`             (force/distance units)
* :math:`\epsilon_c`    (unit less)
* :math:`\gamma`        (force/velocity units)

By default, pair forces are not calculated between bonded particles.
Pair forces can alternatively be overlaid on top of bond forces using
the *overlay/pair* keyword. These settings require specific
:doc:`special_bonds <special_bonds>` settings described in the
restrictions.  Further details can be found in the `:doc: how to
<Howto_BPM>` page on BPMs.

If the *store/local* keyword is used, this fix will track bonds that
break during the simulation. Whenever a bond breaks, data is processed
and transferred to an internal fix labeled *fix_ID*. This allows the
local data to be accessed by other LAMMPS commands.
Following any optional keyword/value arguments, a list of one or more
attributes is specified.  These include the IDs of the two atoms in
the bond. The other attributes for the two atoms include the timestep
during which the bond broke and the current/initial center of mass
position of the two atoms.

Data is continuously accumulated over intervals of *N*
timesteps. At the end of each interval, all of the saved accumulated
data is deleted to make room for new data. Individual datum may
therefore persist anywhere between *1* to *N* timesteps depending on
when they are saved. This data can be accessed using the *fix_ID* and a
:doc:`dump local <dump>` command. To ensure all data is output,
the dump frequency should correspond to the same interval of *N*
timesteps. A dump frequency of an integer multiple of *N* can be used
to regularly output a sample of the accumulated data.

Note that when unbroken bonds are dumped to a file via the
:doc:`dump local <dump>` command, bonds with type 0 (broken bonds)
are not included.
The :doc:`delete_bonds <delete_bonds>` command can also be used to
query the status of broken bonds or permanently delete them, e.g.:

.. code-block:: LAMMPS

   delete_bonds all stats
   delete_bonds all bond 0 remove

----------

Restart and other info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This bond style writes the reference state of each bond to
:doc:`binary restart files <restart>`. Loading a restart
file will properly resume bonds.

The single() function of these pair styles returns 0.0 for the energy
of a pairwise interaction, since energy is not conserved in these
dissipative potentials.

The accumulated data is not written to restart files and should be
output before a restart file is written to avoid missing data.

The internal fix calculates a local vector or local array depending on the
number of input values.  The length of the vector or number of rows in
the array is the number of recorded, lost interactions.  If a single
input is specified, a local vector is produced.  If two or more inputs
are specified, a local array is produced where the number of columns =
the number of inputs.  The vector or array can be accessed by any
command that uses local values from a compute as input.  See the
:doc:`Howto output <Howto_output>` page for an overview of LAMMPS
output options.

The vector or array will be floating point values that correspond to
the specified attribute.

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

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`pair bpm/spring <pair_bpm_spring>`

Default
"""""""

The option defaults are *smooth* = *yes*

----------

.. _Groot4:

**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).
