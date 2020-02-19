.. index:: improper\_style hybrid

improper\_style hybrid command
==============================

Syntax
""""""


.. parsed-literal::

   improper_style hybrid style1 style2 ...

* style1,style2 = list of one or more improper styles

Examples
""""""""


.. parsed-literal::

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

In the improper\_coeff command, the first coefficient sets the improper
style and the remaining coefficients are those appropriate to that
style.  In the example above, the 2 improper\_coeff commands would set
impropers of improper type 1 to be computed with a *harmonic*
potential with coefficients 120.0, 30 for :math:`K`, :math:`\chi_0`.
Improper type 2 would be computed with a *cvff* potential with coefficients
20.0, -1, 2 for K, d, and n, respectively.

If the improper *class2* potential is one of the hybrid styles, it
requires additional AngleAngle coefficients be specified in the data
file.  These lines must also have an additional "class2" argument
added after the improper type.  For improper types which are assigned
to other hybrid styles, use the style name (e.g. "harmonic")
appropriate to that style.  The AngleAngle coeffs for that improper
type will then be ignored.

An improper style of *none* can be specified as the 2nd argument to
the improper\_coeff command, if you desire to turn off certain improper
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
restart file, you need to re-specify improper\_coeff commands.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

**Default:** none
