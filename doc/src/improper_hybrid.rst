.. index:: improper_style hybrid

improper_style hybrid command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   improper_style hybrid style1 style2 ...

* style1,style2 = list of one or more improper styles

Examples
""""""""

.. code-block:: LAMMPS

   improper_style hybrid harmonic helix
   improper_coeff 1 harmonic 120.0 30
   improper_coeff 2 cvff 20.0 -1 2

Description
"""""""""""

The *hybrid* style enables the use of multiple improper styles in one
simulation.  An improper style is assigned to each improper type.  For
example, impropers in a polymer flow (of improper type 1) could be
computed with a *harmonic* potential and impropers in the wall
boundary (of improper type 2) could be computed with a *cvff*
potential.  The assignment of improper type to style is made via the
:doc:`improper_coeff <improper_coeff>` command or in the data file.

In the improper_coeff command, the first coefficient sets the improper
style and the remaining coefficients are those appropriate to that
style.  In the example above, the 2 improper_coeff commands would set
impropers of improper type 1 to be computed with a *harmonic*
potential with coefficients 120.0, 30 for :math:`K`, :math:`\chi_0`.
Improper type 2 would be computed with a *cvff* potential with coefficients
20.0, -1, 2 for K, d, and n, respectively.

If improper coefficients are specified in the data file read via the
:doc:`read_data <read_data>` command, then the same rule applies.
E.g. "harmonic" or "cvff", must be added after the improper type, for
each line in the "Improper Coeffs" section, e.g.

.. parsed-literal::

   Improper Coeffs

   1 harmonic 120.0 30
   2 cvff 20.0 -1 2
   ...

If *class2* is one of the improper hybrid styles, the same rule holds
for specifying additional AngleAngle coefficients either via the input
script or in the data file. I.e. *class2* must be added to each line
after the improper type.  For
lines in the AngleAngle Coeffs section of the data
file for dihedral types that are not *class2*, you must use an
improper style of *skip* as a placeholder, e.g.

.. parsed-literal::

   AngleAngle Coeffs

   1 skip
   2 class2 0.0 0.0 0.0 115.06 130.01 115.06
   ...

Note that it is not necessary to use the improper style *skip* in the
input script, since AngleAngle coefficients
need not be specified at all for improper types that are not *class2*.

An improper style of *none* can be specified as the second argument to
the improper_coeff command, if you desire to turn off certain improper
types.

----------

Restrictions
""""""""""""

This improper style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Unlike other improper styles, the hybrid improper style does not store
improper coefficient info for individual sub-styles in a :doc:`binary restart files <restart>`.
Thus when restarting a simulation from a
restart file, you need to re-specify improper_coeff commands.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

Default
"""""""

none
