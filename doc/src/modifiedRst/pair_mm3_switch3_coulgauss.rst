.. index:: pair\_style mm3/switch3/coulgauss/long

pair\_style mm3/switch3/coulgauss/long command
==============================================

Syntax
""""""


.. parsed-literal::

   pair_style style args

* style = *mm3/switch3/coulgauss/long*
* args = list of arguments for a particular style


.. parsed-literal::

     *mm3/switch3/coulgauss/long* args = cutoff (cutoff2) width
       cutoff  = global cutoff for MM3 (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
       width  = width parameter of the smoothing function (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style mm3/switch3/coulgauss/long    12.0 3.0
   pair_coeff 1  0.2 2.5 1.2

   pair_style mm3/switch3/coulgauss/long   12.0 10.0 3.0
   pair_coeff 1  0.2 2.5 1.2

Description
"""""""""""

The *mm3/switch3/coulgauss* style evaluates the MM3
vdW potential :ref:`(Allinger) <mm3-allinger1989>`

.. math::

  E &=& \epsilon_{ij} \left[ -2.25 \left(\frac{r_{v,ij}}{r_{ij}}\right)^6 + 1.84(10)^5 \exp\left[-12.0 r_{ij}/r_{v,ij}\right] \right] S_3(r_{ij}) \\
  r_{v,ij} &=& r_{v,i} + r_{v,j} \\ 
  \epsilon_{ij} &=& \sqrt{\epsilon_i \epsilon_j}


, which goes smoothly to zero at the cutoff r\_c as defined
by the switching function

.. math::

 S_3(r) = \left\lbrace \begin{array}{ll}
                     1 & \quad\mathrm{if}\quad r < r_\mathrm{c} - w \\
                     3x^2 - 2x^3 & \quad\mathrm{if}\quad r < r_\mathrm{c} \quad\mathrm{with\quad} x=\frac{r_\mathrm{c} - r}{w} \\
                     0 & \quad\mathrm{if}\quad r >= r_\mathrm{c}
                 \end{array} \right.


where w is the width defined in the arguments. This potential
is combined with Coulomb interaction between Gaussian charge densities:

.. math::

  E &=& \frac{q_i q_j \mathrm{erf}\left( r/\sqrt{\gamma_1^2+\gamma_2^2} \right) }{\epsilon r_{ij}}


where qi and qj are the
charges on the 2 atoms, epsilon is the dielectric constant which
can be set by the :doc:`dielectric <dielectric>` command, gamma\_i and gamma\_j
are the widths of the Gaussian charge distribution and erf() is the error-function.
This style has to be used in conjunction with the :doc:`kspace\_style <kspace_style>` command

If one cutoff is specified it is used for both the vdW and Coulomb
terms.  If two cutoffs are specified, the first is used as the cutoff
for the vdW terms, and the second is the cutoff for the Coulombic term.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair\_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands:

* epsilon (energy)
* r\_v (distance)
* gamma (distance)


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

Mixing rules are fixed for this style as defined above.

Shifting the potential energy is not necessary because the switching
function ensures that the potential is zero at the cut-off.

Restrictions
""""""""""""


These styles are part of the USER-YAFF package.  They are only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
