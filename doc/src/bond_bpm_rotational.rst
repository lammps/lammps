.. index:: bond_style bpm/rotational

bond_style bpm/rotational command
=================================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/rotational keyword value attribute1 attribute2 ...

* optional keyword = *overlay/pair* or *store/local* or *smooth* or *break/no*

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

       *break/no*
          indicates that bonds should not break during a run

Examples
""""""""

.. code-block:: LAMMPS

   bond_style bpm/rotational
   bond_coeff 1 1.0 0.2 0.02 0.02 0.20 0.04 0.04 0.04 0.1 0.02 0.002 0.002

   bond_style bpm/rotational store/local myfix 1000 time id1 id2
   dump 1 all local 1000 dump.broken f_myfix[1] f_myfix[2] f_myfix[3]
   dump_modify 1 write_header no

Description
"""""""""""

.. versionadded:: 4May2022

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
displacements with stiffnesses of :math:`k_b` and :math:`k_t`,
respectively.  Details on the calculations of shear displacements and
angular displacements can be found in :ref:`(Wang) <Wang2009>` and
:ref:`(Wang and Mora) <Wang2009b>`.

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
can be optionally multiplied by an extra factor :math:`w` to smoothly
interpolate forces and torques to zero as the bond breaks. This term
is calculated as :math:`w = (1.0 - B^4)`. This smoothing factor can be
added or removed using the *smooth* keyword.

Finally, additional damping forces and torques are applied to the two
particles. A force is applied proportional to the difference in the
normal velocity of particles using a similar construction as
dissipative particle dynamics (:ref:`(Groot) <Groot3>`):

.. math::

   F_D = - \gamma_n w (\hat{r} \bullet \vec{v})

where :math:`\gamma_n` is the damping strength, :math:`\hat{r}` is the
radial normal vector, and :math:`\vec{v}` is the velocity difference
between the two particles. Similarly, tangential forces are applied to
each atom proportional to the relative differences in sliding
velocities with a constant prefactor :math:`\gamma_s` (:ref:`(Wang et
al.) <Wang20152>`) along with their associated torques. The rolling and
twisting components of the relative angular velocities of the two
atoms are also damped by applying torques with prefactors of
:math:`\gamma_r` and :math:`\gamma_t`, respectively.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`k_r`           (force/distance units)
* :math:`k_s`           (force/distance units)
* :math:`k_t`           (force*distance/radians units)
* :math:`k_b`           (force*distance/radians units)
* :math:`f_{r,c}`       (force units)
* :math:`f_{s,c}`       (force units)
* :math:`\tau_{b,c}`    (force*distance units)
* :math:`\tau_{t,c}`    (force*distance units)
* :math:`\gamma_n`      (force/velocity units)
* :math:`\gamma_s`      (force/velocity units)
* :math:`\gamma_r`      (force*distance/velocity units)
* :math:`\gamma_t`      (force*distance/velocity units)

By default, pair forces are not calculated between bonded particles.
Pair forces can alternatively be overlaid on top of bond forces using
the *overlay/pair* keyword. These settings require specific
:doc:`special_bonds <special_bonds>` settings described in the
restrictions.  Further details can be found in the `:doc: how to
<Howto_BPM>` page on BPMs.

.. versionadded:: TBD

If the *break/no* keyword is used, then LAMMPS assumes bonds should not break
during a simulation run. This will prevent some unnecessary calculation.
However, if a bond does break, it will trigger an error.

If the *store/local* keyword is used, an internal fix will track bonds that
break during the simulation. Whenever a bond breaks, data is processed
and transferred to an internal fix labeled *fix_ID*. This allows the
local data to be accessed by other LAMMPS commands. Following this optional
keyword, a list of one or more attributes is specified.  These include the
IDs of the two atoms in the bond. The other attributes for the two atoms
include the timestep during which the bond broke and the current/initial
center of mass position of the two atoms.

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
:doc:`binary restart files <restart>`. Loading a restart file will
properly resume bonds. However, the reference state is NOT
written to data files. Therefore reading a data file will not
restore bonds and will cause their reference states to be redefined.

If the *store/local* option is used, an internal fix will calculate
a local vector or local array depending on the number of input values.
The length of the vector or number of rows in the array is the number
of recorded, broken bonds.  If a single input is specified, a local
vector is produced. If two or more inputs are specified, a local array
is produced where the number of columns = the number of inputs.  The
vector or array can be accessed by any command that uses local values
from a compute as input. See the :doc:`Howto output <Howto_output>` page
for an overview of LAMMPS output options.

The vector or array will be floating point values that correspond to
the specified attribute.

The single() function of this bond style returns 0.0 for the energy
of a bonded interaction, since energy is not conserved in these
dissipative potentials.  It also returns only the normal component of
the bonded interaction force.  However, the single() function also
calculates 7 extra bond quantities.  The first 4 are data from the
reference state of the bond including the initial distance between particles
:math:`r_0` followed by the :math:`x`, :math:`y`, and :math:`z` components
of the initial unit vector pointing to particle I from particle J. The next 3
quantities (5-7) are the  :math:`x`, :math:`y`, and :math:`z` components
of the total force, including normal and tangential contributions, acting
on particle I.

These extra quantities can be accessed by the :doc:`compute bond/local <compute_bond_local>`
command, as *b1*, *b2*, ..., *b7*\ .

Restrictions
""""""""""""

This bond style is part of the BPM package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

By default if pair interactions are to be disabled, this bond style
requires setting

.. code-block:: LAMMPS

   special_bonds lj 0 1 1 coul 1 1 1

and :doc:`newton <newton>` must be set to bond off.  If the
*overlay/pair* option is used, this bond style alternatively requires
setting

.. code-block:: LAMMPS

   special_bonds lj/coul 1 1 1

The *bpm/rotational* style requires :doc:`atom style bpm/sphere <atom_style>`.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`fix nve/bpm/sphere <fix_nve_bpm_sphere>`

Default
"""""""

The option defaults are *smooth* = *yes*

----------

.. _Wang2009:

**(Wang)** Wang, Acta Geotechnica, 4,
p 117-127 (2009).

.. _Wang2009b:

**(Wang and Mora)** Wang, Mora, Advances in Geocomputing,
119, p 183-228 (2009).

.. _Groot3:

**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).

.. _Wang20152:

**(Wang et al, 2015)** Wang, Y., Alonso-Marroquin, F., & Guo,
W. W. (2015).  Rolling and sliding in 3-D discrete element
models. Particuology, 23, 49-55.
