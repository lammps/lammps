.. index:: bond_style bpm/spring/plastic

bond_style bpm/spring/plastic command
=====================================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/spring/plastic keyword value attribute1 attribute2 ...

* optional keyword = *overlay/pair* or *store/local* or *write/history* or *read/history* or *smooth* or *normalize* or *break*

  .. parsed-literal::

       *store/local* values = fix_ID N attributes ...
          * fix_ID = ID of associated internal fix to store data
          * N = prepare data for output every this many timesteps
          * attributes = zero or more of the below attributes may be appended

            *id1, id2* = IDs of two atoms in the bond
            *time* = the timestep the bond broke
            *x, y, z* = the center of mass position of the two atoms when the bond broke (distance units)
            *x/ref, y/ref, z/ref* = the initial center of mass position of the two atoms (distance units)

       *write/history* values = fix_ID N
          * fix_ID = ID of associated internal fix to write data
          * N = prepare data for output every this many timesteps

       *read/history* values = filename
          * filename = name of history file to read data from

       *overlay/pair* value = *yes* or *no*
          bonded particles will still interact with pair forces

       *smooth* value = *yes* or *no*
          smooths bond forces near the breaking point

       *normalize* value = *yes* or *no*
          normalizes bond forces by the reference length

       *break* value = *yes* or *no*
          indicates whether bonds break during a run

Examples
""""""""

.. code-block:: LAMMPS

   bond_style bpm/spring/plastic
   bond_coeff 1 1.0 0.05 0.1 0.02

   bond_style bpm/spring/plastic myfix 1000 time id1 id2
   dump 1 all local 1000 dump.broken f_myfix[1] f_myfix[2] f_myfix[3]
   dump_modify 1 write_header no

   bond_style bpm/spring/plastic write/history myfix 1000 read/history bond.ref
   dump 1 all local 1000 bond*.ref f_myfix[*]

Description
"""""""""""

.. versionadded:: 2Apr2025

The *bpm/spring/plastic* bond style computes forces based on
deviations from the initial reference state of the two atoms and the
strain history.  The reference length of the bond :math:`r_0` is stored
by each bond when it is first computed in the setup of a run. Initially,
the equilibrium length of each bond :math:`r_\mathrm{eq}` is set equal
to :math:`r_0` but can evolve. data is then preserved across run commands
and is written to :doc:`binary restart files <restart>` such that restarting
the system will not modify either of these quantities.

This bond style only applies central-body forces which conserve the
translational and rotational degrees of freedom of a bonded set of
particles. The force has a magnitude of

.. math::

   F = -k (r_\mathrm{eq} - r) w

where :math:`k` is a stiffness, :math:`r` is the current distance between
the two particles, and :math:`w` is an optional smoothing factor discussed
below. If the bond stretches beyond a strain of :math:`\epsilon_p` in compression
or extension, it will plastically activate and :math:`r_\mathrm{eq}` will evolve
to ensure :math:`|(r-r_\mathrm{eq})/r_\mathrm{eq}|` never exceeds :math:`\epsilon_p`.
Therefore, if a bond is continually loaded in either tension or compression, the
force will initially grow elastically before plateauing. See
:ref:`(Clemmer4) <plastic-Clemmer>` for more details on these mechanics.

Bonds will break at a strain of :math:`\epsilon_c`.  This is done by setting
the bond type to 0 such that forces are no longer computed.

An additional damping force is applied to the bonded
particles.  This forces is proportional to the difference in the
normal velocity of particles:

.. math::

   F_D = - \gamma w (\hat{r} \bullet \vec{v})

where :math:`\gamma` is the damping strength, :math:`\hat{r}` is the
radial normal vector, and :math:`\vec{v}` is the velocity difference
between the two particles.

The smoothing factor :math:`w`  is constructed such that forces smoothly
go to zero, avoiding discontinuities, as bonds approach the critical
breaking strain

.. math::

   w = 1.0 - \left( \frac{r - r_0}{r_0 \epsilon_c} \right)^8 .

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands:

* :math:`k`             (force/distance units)
* :math:`\epsilon_c`    (unitless)
* :math:`\gamma`        (force/velocity units)
* :math:`\epsilon_p`    (unitless)

See the :doc:`bpm/spring doc page <bond_bpm_spring>` for information on
the *smooth*, *normalize*, *break*, *overlay/pair*, *store/local*,
*write/history*, and *read/history* keywords.

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

This bond style writes the history data (e.g. reference state) of each
bond and its per-type coefficients to :doc:`binary restart files <restart>`.
Loading a restart file restores bonds and their history data.  The history
data is NOT written to data files.  Reading a data file will therefore not
restore bond history such that bond reference states will be redefined.
Alternatively, bond history data can be saved and restored using the
*write/history* and *read/history* options.

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

Any settings with the *store/local* option are not saved to a restart
file and must be redefined.

If the *write/history* keyword is used, an internal fix will process
and transfer the internal bond history (e.g. :math:`r_0`) to an
internal fix labeled *fix_ID*.  This allows the internal bond history data,
as well as the IDs of the two atoms in the bond,  to be accessed by other
LAMMPS commands, in particular, :doc:`dump local <dump>`.

If the *read/history* keyword is used, history data is read from the file
labeled *filename*.  This is expected to be in the format of a LAMMPS
local dump file. The first two columns of the history file must contain
the IDs of the two atoms in the bond, with the remaining columns corresponding
 to the internal bond data.  For more details on formatting the remaining file
 see :doc:`Howto bpm <Howto_bpm>`.

The potential energy and the single() function of this bond style
returns zero.  The single() function also calculates two extra bond
quantities, the initial distance :math:`r_0` and the current equilibrium
length :math:`r_{eq}`. These extra quantities can be accessed by the
:doc:`compute bond/local <compute_bond_local>` command as *b1* and *b2*,
respectively.

Restrictions
""""""""""""

This bond style is part of the BPM package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

To handle breaking bonds, BPM bond styles have extra requirements for
special bonds. If bonds cannot break (*break no*), then one can use any
special bond weights. Otherwise, restrictions depend on whether pair
forces are overlaid (*pair/overlay yes*). If so, then all weights must
be one:

.. code-block:: LAMMPS

   special_bonds lj/coul 1 1 1

If pair forces are disabled (*pair/overlay no*), the default, then the
weights must be

.. code-block:: LAMMPS

   special_bonds lj 0 1 1 coul 1 1 1

and :doc:`newton <newton>` must be set to bond off.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`bond bpm/spring <bond_bpm_spring>`

Default
"""""""

The option defaults are *overlay/pair* = *no*, *smooth* = *yes*, *normalize* = *no*, and *break* = *yes*

----------

.. _plastic-Clemmer:

**(Clemmer4)** Clemmer and Lechman, Powder Technology (2025).

