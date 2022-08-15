.. index:: angle_style mesocnt

angle_style mesocnt command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   angle_style mesocnt

Examples
""""""""

.. code-block:: LAMMPS

   angle_style mesocnt
   angle_coeff 1 buckling C 10 10 20.0
   angle_coeff 4 harmonic C 8 4 10.0
   angle_coeff 2 buckling custom 400.0 50.0 5.0
   angle_coeff 1 harmonic custom 300.0

Description
"""""""""""

The *mesocnt* angle style uses the potential

.. math::

   E = K_\text{H} \Delta \theta^2 \qquad |\Delta \theta| < \Delta \theta_\text{B} \\ 
   E = K_\text{H} \theta_\text{B}^2 + K_\text{B} (\Delta \theta - \theta_\text{B}) \qquad |\Delta \theta| \geq \Delta \theta_\text{B}

where :math:`\Delta \theta = \theta - \pi` is the bending angle of the nanotube, :math:`K_\text{H}` and :math:`K_\text{B}` are prefactors for the harmonic and linear regime respectively and :math:`\theta_\text{B}` is the buckling angle. Note that the usual 1/2 factor for the harmonic potential is included in :math:`K_\text{H}`.

The style implements parametrisation presets of :math:`K_\text{H}`, :math:`K_\text{B}` and :math:`\theta_\text{B}` for mesoscopic simulations of 
carbon nanotubes based on the atomistic simulations of :ref:`(Zhigilei) <Zhigilei>` and buckling considerations of :ref:`(Volkov) <Volkov>`.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the examples above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* mode = *buckling* or *harmonic*
* preset = *C* or *custom*
* additional parameters depending on preset

If mode *harmonic* is chosen, the potential is simply harmonic and does not switch to the linear term when the buckling angle is reached. In *buckling* mode, the full piecewise potential is used.

Preset *C* is for carbon nanotubes, and the additional parameters are:

* chiral index :math:`n` (unitless)
* chiral index :math:`m` (unitless)
* :math:`r_0` (distance)

Here, :math:`r_0` is the equilibrium distance of the bonds included in the angle, see :doc:`bond_style mesocnt <bond_mesocnt>`.

In harmonic mode with preset *custom*, the additional paramter is:

* :math:`K_\text{H}` (energy)

Hence, this setting is simply a wrapper for :doc:`bond_style harmonic <bond_harmonic>` with an equilibrium angle of 180 degrees.

In harmonic mode with preset *custom*, the additional paramters are:

* :math:`K_\text{H}` (energy)
* :math:`K_\text{B}` (energy)
* :math:`\theta_\text{B}` (degrees)

:math:`\theta_\text{B}` is specified in degrees, but LAMMPS converts it to
radians internally; hence :math:`K_\text{H}` is effectively energy per
radian\^2 and :math:`K_\text{B}` is energy per radian.

----------

In *buckling* mode, this angle style adds the *buckled* property to all atoms in the simulation, which is an integer flag indicating whether the bending angle at a given atom has exceeded :math:`\theta_\text{B}`. It can be accessed as an atomic variable, e.g. for custom dump commands, as *i_buckled*. 

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
MOLECULE and MESONT packages.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none

----------

.. _Zhigilei:

**(Zhigilei)** Zhigilei, Wei and Srivastava, Phys. Rev. B 71, 165417 (2005).

.. _Volkov:

**(Volkov)** Volkov and Zhigilei, ACS Nano 4, 10, 6187–6195 (2010).
