.. index:: fix rhok

fix rhok command
================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID rhok nx ny nz K a

* ID, group-ID are documented in :doc:`fix <fix>` command
* nx, ny, nz = k-vector of collective density field
* K = spring constant of bias potential
* a = anchor point of bias potential

Examples
""""""""

.. code-block:: LAMMPS

   fix bias all rhok 16 0 0 4.0 16.0
   fix 1 all npt temp 0.8 0.8 4.0 z 2.2 2.2 8.0
   # output of 4 values from fix rhok: U_bias rho_k_RE  rho_k_IM  \|rho_k\|
   thermo_style custom step temp pzz lz f_bias f_bias[1] f_bias[2] f_bias[3]

Description
"""""""""""

The fix applies a force to atoms given by the potential

.. math::

   U  = &  \frac{1}{2} K (|\rho_{\vec{k}}| - a)^2 \\
   \rho_{\vec{k}}  = & \sum_j^N \exp(-i\vec{k} \cdot \vec{r}_j )/\sqrt{N} \\
   \vec{k}  = & (2\pi n_x /L_x , 2\pi n_y  /L_y , 2\pi n_z/L_z )

as described in :ref:`(Pedersen) <Pedersen>`.

This field, which biases configurations with long-range order, can be
used to study crystal-liquid interfaces and determine melting
temperatures :ref:`(Pedersen) <Pedersen>`.

An example of using the interface pinning method is located in the
*examples/PACKAGES/rhok* directory.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by
this fix to add the potential energy calculated by the fix to the
global potential energy of the system as part of :doc:`thermodynamic
output <thermo_style>`.  The default setting for this fix is
:doc:`fix_modify energy no <fix_modify>`.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the potential
energy discussed in the preceding paragraph.  The scalar stored by
this fix is "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the EXTRA-FIX package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`thermo_style <thermo_style>`

Default
"""""""

none

----------

.. _Pedersen:

**(Pedersen)** Pedersen, J. Chem. Phys., 139, 104102 (2013).
