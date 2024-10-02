.. index:: bond_style rheo/shell

bond_style rheo/shell command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style rheo/shell keyword value attribute1 attribute2 ...

* required keyword = *t/form*
* optional keyword = *store/local*

  .. parsed-literal::

       *t/form* value = formation time for a bond (time units)

       *store/local* values = fix_ID N attributes ...
          * fix_ID = ID of associated internal fix to store data
          * N = prepare data for output every this many timesteps
          * attributes = zero or more of the below attributes may be appended

            *id1, id2* = IDs of two atoms in the bond
            *time* = the timestep the bond broke
            *x, y, z* = the center of mass position of the two atoms when the bond broke (distance units)
            *x/ref, y/ref, z/ref* = the initial center of mass position of the two atoms (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   bond_style rheo/shell t/form 10.0
   bond_coeff 1 1.0 0.05 0.1

Description
"""""""""""

.. versionadded:: 29Aug2024

The *rheo/shell* bond style is designed to work with
:doc:`fix rheo/oxidation <fix_rheo_oxidation>` which creates candidate
bonds between eligible surface or near-surface particles. When a bond
is first created, it computes no forces and starts a timer. Forces are
not computed until the timer reaches the specified bond formation time,
*t/form*, and the bond is enabled and applies forces. If the two particles
move outside of the maximum bond distance or move into the bulk before
the timer reaches *t/form*, the bond automatically deletes itself. This
deletion is not recorded as a broken bond in the optional *store/local* fix.

Before bonds are enabled, they are still treated as regular bonds by
all other parts of LAMMPS. This means they are written to data files
and counted in computes such as :doc:`nbond/atom <compute_nbond_atom>`.
To only count enabled bonds, use the *nbond/shell* attribute in
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`.

When enabled, the bond then computes forces based on deviations from
the initial reference state of the two atoms much like a BPM style
bond (as further discussed in the :doc:`BPM howto page <Howto_bpm>`).
The reference state is stored by each bond when it is first enabled.
Data is then preserved across run commands and is written to
:doc:`binary restart files <restart>` such that restarting the system
will not reset the reference state of a bond or the timer.

This bond style is based on a model described in
:ref:`(Clemmer) <rheo_clemmer>`. The force has a magnitude of

.. math::

   F = 2 k (r - r_0) + \frac{2 k}{r_0^2 \epsilon_c^2} (r - r_0)^3

where :math:`k` is a stiffness, :math:`r` is the current distance
and :math:`r_0` is the initial distance between the two particles, and
:math:`\epsilon_c` is maximum strain beyond which a bond breaks. This
is done by setting the bond type to 0 such that forces are no longer
computed.

A damping force proportional to the difference in the normal velocity
of particles is also applied to bonded particles:

.. math::

   F_D = - \gamma w (\hat{r} \bullet \vec{v})

where :math:`\gamma` is the damping strength, :math:`\hat{r}` is the
displacement normal vector, and :math:`\vec{v}` is the velocity difference
between the two particles.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands:

* :math:`k`             (force/distance units)
* :math:`\epsilon_c`    (unit less)
* :math:`\gamma`        (force/velocity units)

Unlike other BPM-style bonds, this bond style does not update special
bond settings when bonds are created or deleted. This bond style also
does not enforce specific :doc:`special_bonds <special_bonds>` settings.
This behavior is purposeful such :doc:`RHEO pair <pair_rheo>` forces
and heat flows are still calculated.

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
:doc:`binary restart files <restart>`. Loading a restart
file will properly restore bonds. However, the reference state is NOT
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
dissipative potentials.  The single() function also calculates two
extra bond quantities, the initial distance :math:`r_0` and a time.
These extra quantities can be accessed by the
:doc:`compute bond/local <compute_bond_local>` command as *b1* and *b2*\ .

Restrictions
""""""""""""

This bond style is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`fix rheo/oxidation <fix_rheo_oxidation>`

Default
"""""""

NA

----------

.. _rheo_clemmer:

**(Clemmer)** Clemmer, Pierce, O'Connor, Nevins, Jones, Lechman, Tencer, Appl. Math. Model., 130, 310-326 (2024).
