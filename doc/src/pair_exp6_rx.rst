.. index:: pair_style exp6/rx

pair_style exp6/rx command
==========================

pair_style exp6/rx/kk command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style exp6/rx cutoff ...

* cutoff = global cutoff for DPD interactions (distance units)
* weighting = fractional or molecular (optional)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style exp6/rx 10.0
   pair_style exp6/rx 10.0 fractional
   pair_style exp6/rx 10.0 molecular
   pair_coeff * * exp6.params h2o h2o exponent 1.0 1.0 10.0
   pair_coeff * * exp6.params h2o 1fluid exponent 1.0 1.0 10.0
   pair_coeff * * exp6.params 1fluid 1fluid exponent 1.0 1.0 10.0
   pair_coeff * * exp6.params 1fluid 1fluid none 10.0
   pair_coeff * * exp6.params 1fluid 1fluid polynomial filename 10.0

Description
"""""""""""

Style *exp6/rx* is used in reaction DPD simulations, where the
coarse-grained (CG) particles are composed of *m* species whose
reaction rate kinetics are determined from a set of *n* reaction rate
equations through the :doc:`fix rx <fix_rx>` command.  The species of
one CG particle can interact with a species in a neighboring CG
particle through a site-site interaction potential model.  The
*exp6/rx* style computes an exponential-6 potential given by

.. math::

   U_{ij}(r) = \frac{\epsilon}{\alpha-6}\{6\exp[\alpha(1-\frac{r_{ij}}{R_{m}})]-\alpha(\frac{R_{m}}{r_{ij}})^6\}

where the :math:`\epsilon` parameter determines the depth of the
potential minimum located at :math:`R_m`, and :math:`\alpha` determines
the softness of the repulsion.

The coefficients must be defined for each species in a given particle
type via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, where the first argument is the filename that includes the
exponential-6 parameters for each species.  The file includes the
species tag followed by the :math:`\alpha, \epsilon` and :math:`R_m`
parameters. The format of the file is described below.

The second and third arguments specify the site-site interaction
potential between two species contained within two different
particles.  The species tags must either correspond to the species
defined in the reaction kinetics files specified with the :doc:`fix rx <fix_rx>` command or they must correspond to the tag "1fluid",
signifying interaction with a product species mixture determined
through a one-fluid approximation.  The interaction potential is
weighted by the geometric average of either the mole fraction concentrations
or the number of molecules associated with the interacting coarse-grained
particles (see the *fractional* or *molecular* weighting pair style options).
The coarse-grained potential is stored before and after the
reaction kinetics solver is applied, where the difference is defined
to be the internal chemical energy (uChem).

The fourth argument specifies the type of scaling that will be used
to scale the EXP-6 parameters as reactions occur.  Currently, there
are three scaling options:  *exponent*\ , *polynomial* and *none*\ .

Exponent scaling requires two additional arguments for scaling
the :math:`R_m` and :math:`\epsilon` parameters, respectively.  The scaling factor
is computed by phi\^exponent, where phi is the number of molecules
represented by the coarse-grain particle and exponent is specified
as a pair coefficient argument for :math:`R_m` and :math:`\epsilon`, respectively.
The :math:`R_m` and :math:`\epsilon` parameters are multiplied by the scaling
factor to give the scaled interaction parameters for the CG particle.

Polynomial scaling requires a filename to be specified as a pair
coeff argument.  The file contains the coefficients to a fifth order
polynomial for the :math:`\alpha`, :math:`\epsilon` and :math:`R_m` parameters that depend
upon phi (the number of molecules represented by the CG particle).
The format of a polynomial file is provided below.

The *none* option to the scaling does not have any additional pair coeff
arguments.  This is equivalent to specifying the *exponent* option with
:math:`R_m` and :math:`\epsilon` exponents of 0.0 and 0.0, respectively.

The final argument specifies the interaction cutoff (optional).

----------

The format of a tabulated file is as follows (without the
parenthesized comments):

.. parsed-literal::

   # exponential-6 parameters for various species      (one or more comment or blank lines)

   h2o  exp6  11.00 0.02 3.50                          (species, exp6, alpha, Rm, epsilon)
   no2  exp6  13.60 0.01 3.70
   ...
   co2  exp6  13.00 0.03 3.20

The format of the polynomial scaling file as follows (without the
parenthesized comments):

.. parsed-literal::

   # POLYNOMIAL FILE          (one or more comment or blank lines)

   #  General Functional Form:
   #  A\*phi\^5 + B\*phi\^4 + C\*phi\^3 + D\*phi\^2 + E\*phi + F
   #
   #  Parameter  A        B         C        D         E        F
                              (blank)
   alpha        0.0000   0.00000   0.00008  0.04955  -0.73804  13.63201
   epsilon      0.0000   0.00478  -0.06283  0.24486  -0.33737   2.60097
   rm           0.0001  -0.00118  -0.00253  0.05812  -0.00509   1.50106

A section begins with a non-blank line whose 1st character is not a
"#"; blank lines or lines starting with "#" can be used as comments
between sections.

Following a blank line, the next N lines list the species and their
corresponding parameters.  The first argument is the species tag, the
second argument is the exp6 tag, the 3rd argument is the :math:`\alpha`
parameter (energy units), the 4th argument is the :math:`\epsilon` parameter
(energy-distance\^6 units), and the 5th argument is the :math:`R_m` parameter
(distance units).  If a species tag of "1fluid" is listed as a pair
coefficient, a one-fluid approximation is specified where a
concentration-dependent combination of the parameters is computed
through the following equations:

.. math::

   R_{m}^{3} = & \sum_{a}\sum_{b} x_{a}x_{b}R_{m,ab}^{3} \\
   \epsilon  = & \frac{1}{R_{m}^{3}}\sum_{a}\sum_{b} x_{a}x_{b}\epsilon_{ab}R_{m,ab}^{3} \\
   \alpha    = & \frac{1}{\epsilon R_{m}^{3}}\sum_{a}\sum_{b} x_{a}x_{b}\alpha_{ab}\epsilon_{ab}R_{m,ab}^{3}

where

.. math::

   \epsilon_{ab} = & \sqrt{\epsilon_{a}\epsilon_{b}} \\
   R_{m,ab}      = & \frac{R_{m,a}+R_{m,b}}{2} \\
   \alpha_{ab}   = & \sqrt{\alpha_{a}\alpha_{b}}

and :math:`x_a` and :math:`x_b` are the mole fractions of a and b, respectively, which
comprise the gas mixture.

----------

**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair style does not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

This style does not support the :doc:`pair_modify <pair_modify>` shift option
for the energy of the exp() and 1/r\^6 portion of the pair interaction.

This style does not support the pair_modify tail option for adding long-range
tail corrections to energy and pressure for the A,C terms in the
pair interaction.

----------

Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.

----------

Restrictions
""""""""""""

This command is part of the USER-DPD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** fractional weighting
