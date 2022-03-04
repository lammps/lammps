.. index:: improper_style inversion/harmonic

improper_style inversion/harmonic command
=========================================

Syntax
""""""

.. code-block:: LAMMPS

   improper_style inversion/harmonic

Examples
""""""""

.. code-block:: LAMMPS

   improper_style inversion/harmonic
   improper_coeff 1 18.776340 0.000000

Description
"""""""""""

The *inversion/harmonic* improper style follows the Wilson-Decius
out-of-plane angle definition and uses an harmonic potential:

.. math::

   E = K \left(\omega - \omega_0\right)^2

where :math:`K` is the force constant and :math:`\omega` is the angle
evaluated for all three axis-plane combinations centered around the atom I.
For the IL axis and the IJK plane :math:`\omega` looks as follows:

.. image:: JPG/umbrella.jpg
   :align: center

Note that the *inversion/harmonic* angle term evaluation differs to
the :doc:`improper_umbrella <improper_umbrella>` due to the cyclic
evaluation of all possible angles :math:`\omega`.

The following coefficients must be defined for each improper type via
the :doc:`improper_coeff <improper_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`K` (energy)
* :math:`\omega_0` (degrees)

If :math:`\omega_0 = 0` the potential term has a single minimum for
the planar structure.  Otherwise it has two minima at +/- :math:`\omega_0`,
with a barrier in between.

----------

Restrictions
""""""""""""

This improper style can only be used if LAMMPS was built with the
MOFFF package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

Default
"""""""

none
