.. index:: pair\_style lj/mdf

pair\_style lj/mdf command
==========================

pair\_style buck/mdf command
============================

pair\_style lennard/mdf command
===============================

Syntax
""""""


.. parsed-literal::

   pair_style style args

* style = *lj/mdf* or *buck/mdf* or *lennard/mdf*
* args = list of arguments for a particular style
  
  .. parsed-literal::
  
       *lj/mdf* args = cutoff1 cutoff2
         cutoff1 = inner cutoff for the start of the tapering function
         cutoff1 = out cutoff for the end of the tapering function
       *buck/mdf* args = cutoff1 cutoff2
         cutoff1 = inner cutoff for the start of the tapering function
         cutoff1 = out cutoff for the end of the tapering function
       *lennard/mdf* args = cutoff1 cutoff2
         cutoff1 = inner cutoff for the start of the tapering function
         cutoff1 = out cutoff for the end of the tapering function



Examples
""""""""


.. parsed-literal::

   pair_style lj/mdf 2.5 3.0
   pair_coeff \* \* 1.0 1.0
   pair_coeff 1 1 1.1 2.8 3.0 3.2

   pair_style buck 2.5 3.0
   pair_coeff \* \* 100.0 1.5 200.0
   pair_coeff \* \* 100.0 1.5 200.0 3.0 3.5

   pair_style lennard/mdf 2.5 3.0
   pair_coeff \* \* 1.0 1.0
   pair_coeff 1 1 1021760.3664 2120.317338 3.0 3.2

Description
"""""""""""

The *lj/mdf*\ , *buck/mdf* and *lennard/mdf* compute the standard 12-6
Lennard-Jones and Buckingham potential with the addition of a taper
function that ramps the energy and force smoothly to zero between an
inner and outer cutoff.

.. math::

   E_{smooth}(r) = E(r)*f(r)


The tapering, *f(r)*\ , is done by using the Mei, Davenport, Fernando
function :ref:`(Mei) <Mei>`.

.. math::

   f(r) & = 1.0  \qquad \qquad \mathrm{for} \qquad r < r_m \\
   f(r) & = (1 - x)^3*(1+3x+6x^2) \quad \mathrm{for} \qquad r_m < r < r_{cut} \\
   f(r) & = 0.0  \qquad \qquad \mathrm{for} \qquad  r >= r_{cut} \\

where

.. math::

   x = \frac{(r-r_m)}{(r_{cut}-r_m)}


Here :math:`r_m` is the inner cutoff radius and :math:`r_{cut}` is the
outer cutoff radius.


----------


For the *lj/mdf* pair\_style, the potential energy, *E(r)*\ , is the
standard 12-6 Lennard-Jones written in the epsilon/sigma form:

.. math::

   E(r) = 4\epsilon\biggl[\bigl(\frac{\sigma}{r}\bigr)^{12} - \bigl(\frac{\sigma}{r}\bigr)^6\biggr]


Either the first two or all of the following coefficients must be
defined for each pair of atoms types via the pair\_coeff command as in
the examples above, or in the data file read by the :doc:`read_data
<read_data>`. The two cutoffs default to the global values and
:math:`\epsilon` and :math:`\sigma` can also be determined by mixing as
described below:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* :math:`r_m` (distance units)
* :math:`r_{cut}` (distance units)

----------

For the *buck/mdf* pair\_style, the potential energy, *E(r)*\ , is the
standard Buckingham potential with three required coefficients.
The two cutoffs can be omitted and default to the corresponding
global values:

.. math::

   E(r) = A e^{(-r/\rho)} -\frac{C}{r^6}


* *A* (energy units)
* :math:`\rho` (distance units)
* *C* (energy-distance\^6 units)
* :math:`r_m` (distance units)
* :math:`r_{cut}` (distance units)

----------

For the *lennard/mdf* pair\_style, the potential energy, *E(r)*\ , is the
standard 12-6 Lennard-Jones written in the A/B form:

.. math::

   E(r) = \frac{A}{r^{12}} - \frac{B}{r^{6}}


The following coefficients must be defined for each pair of atoms
types via the pair\_coeff command as in the examples above, or in the
data file read by the read\_data commands, or by mixing as described below.
The two cutoffs default to their global values and must be either both
given or both left out:

* *A* (energy-distance\^12 units)
* *B* (energy-distance\^6 units)
* :math:`r_m` (distance units)
* :math:`r_{cut}` (distance units)

----------

**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, the :math:`\epsilon` and
:math:`sigma` coefficients and cutoff distances for the lj/mdf pair
style can be mixed.  The default mix value is *geometric*\ .  See the
"pair\_modify" command for details. The other two pair styles buck/mdf
and lennard/mdf do not support mixing, so all I,J pairs of coefficients
must be specified explicitly.

None of the lj/mdf, buck/mdf, or lennard/mdf pair styles supports
the :doc:`pair_modify <pair_modify>` shift option or long-range
tail corrections to pressure and energy.

These styles write their information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

These styles can only be used via the *pair* keyword of the :doc:`run_style respa <run_style>` command.  They do not support the *inner*\ ,
*middle*\ , *outer* keywords.

----------

Restrictions
""""""""""""

These pair styles can only be used if LAMMPS was built with the
USER-MISC package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none


----------


.. _Mei:



**(Mei)** Mei, Davenport, Fernando, Phys Rev B, 43 4653 (1991)
