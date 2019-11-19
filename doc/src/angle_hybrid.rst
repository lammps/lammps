.. index:: angle_style hybrid

angle_style hybrid command
==========================

Syntax
""""""


.. code-block:: LAMMPS

   angle_style hybrid style1 style2 ...

* style1,style2 = list of one or more angle styles

Examples
""""""""


.. code-block:: LAMMPS

   angle_style hybrid harmonic cosine
   angle_coeff 1 harmonic 80.0 30.0
   angle_coeff 2* cosine 50.0

Description
"""""""""""

The *hybrid* style enables the use of multiple angle styles in one
simulation.  An angle style is assigned to each angle type.  For
example, angles in a polymer flow (of angle type 1) could be computed
with a *harmonic* potential and angles in the wall boundary (of angle
type 2) could be computed with a *cosine* potential.  The assignment
of angle type to style is made via the :doc:`angle_coeff <angle_coeff>`
command or in the data file.

In the :doc:`angle_coeff <angle_coeff>` commands, the name of an angle style must be added
after the angle type, with the remaining coefficients being those
appropriate to that style.  In the example above, the 2 angle\_coeff
commands set angles of angle type 1 to be computed with a *harmonic*
potential with coefficients 80.0, 30.0 for :math:`K`, :math:`\theta_0`.  All other angle
types :math:`(2 - N)` are computed with a *cosine* potential with coefficient
50.0 for :math:`K`.

If angle coefficients are specified in the data file read via the
:doc:`read_data <read_data>` command, then the same rule applies.
E.g. "harmonic" or "cosine", must be added after the angle type, for each
line in the "Angle Coeffs" section, e.g.


.. parsed-literal::

   Angle Coeffs

   1 harmonic 80.0 30.0
   2 cosine 50.0
   ...

If *class2* is one of the angle hybrid styles, the same rule holds for
specifying additional BondBond (and BondAngle) coefficients either via
the input script or in the data file.  I.e. *class2* must be added to
each line after the angle type.  For lines in the BondBond (or
BondAngle) section of the data file for angle types that are not
*class2*\ , you must use an angle style of *skip* as a placeholder, e.g.


.. parsed-literal::

   BondBond Coeffs

   1 skip
   2 class2 3.6512 1.0119 1.0119
   ...

Note that it is not necessary to use the angle style *skip* in the
input script, since BondBond (or BondAngle) coefficients need not be
specified at all for angle types that are not *class2*\ .

An angle style of *none* with no additional coefficients can be used
in place of an angle style, either in a input script :doc:`angle_coeff <angle_coeff>`
command or in the data file, if you desire to turn off interactions
for specific angle types.


----------


Restrictions
""""""""""""


This angle style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Unlike other angle styles, the hybrid angle style does not store angle
coefficient info for individual sub-styles in a :doc:`binary restart files <restart>`.  Thus when restarting a simulation from a restart
file, you need to re-specify :doc:`angle_coeff <angle_coeff>` commands.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

**Default:** none
