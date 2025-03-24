.. index:: bond_style bpm/spring/plastic

bond_style bpm/spring/plastic command
=====================================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/spring/plastic keyword value attribute1 attribute2 ...

* optional keyword = *overlay/pair* or *store/local* or *smooth* or *normalize* or *break*

  .. parsed-literal::

       *store/local* values = fix_ID N attributes ...
          * fix_ID = ID of associated internal fix to store data
          * N = prepare data for output every this many timesteps
          * attributes = zero or more of the below attributes may be appended

            *id1, id2* = IDs of two atoms in the bond
            *time* = the timestep the bond broke
            *x, y, z* = the center of mass position of the two atoms when the bond broke (distance units)
            *x/ref, y/ref, z/ref* = the initial center of mass position of the two atoms (distance units)

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

Description
"""""""""""

.. versionadded:: TBD

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
:ref:`(Clemmer) <plastic-Clemmer>` for more details on these mechanics.

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
the *smooth*, *normalize*, *break*, *overlay/pair*, and *store/local*
keywords.

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

This bond style writes the reference state and plastic history of each
bond to :doc:`binary restart files <restart>`. Loading a restart file
will properly restore bonds. However, the reference state is NOT written
to data files.  Therefore reading a data file will not restore bonds and
will cause their reference states to be redefined.

The potential energy and the single() function of this bond style
returns zero.  The single() function also calculates two extra bond
quantities, the initial distance :math:`r_0` and the current equilibrium
length :math:`r_eq`. These extra quantities can be accessed by the
:doc:`compute bond/local <compute_bond_local>` command as *b1* and *b2*,
respectively.

Restrictions
""""""""""""

This bond style is part of the BPM package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

By default if pair interactions between bonded atoms are to be disabled,
this bond style requires setting

.. code-block:: LAMMPS

   special_bonds lj 0 1 1 coul 1 1 1

and :doc:`newton <newton>` must be set to bond off.  If the *overlay/pair*
keyword is set to *yes*, this bond style alternatively requires setting

.. code-block:: LAMMPS

   special_bonds lj/coul 1 1 1

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`bond bpm/spring <bond_bpm_spring>`

Default
"""""""

The option defaults are *overlay/pair* = *no*, *smooth* = *yes*, *normalize* = *no*, and *break* = *yes*

----------

.. _plastic-Clemmer:

**(Clemmer)** Clemmer and Lechman, Powder Technology (2025).

