.. index:: pair_style hbond/dreiding/lj/angleoffset
.. index:: pair_style hbond/dreiding/lj/angleoffset/omp
.. index:: pair_style hbond/dreiding/morse/angleoffset
.. index:: pair_style hbond/dreiding/morse/angleoffest/omp

pair_style hbond/dreiding/lj/angleoffset command
====================================

Accelerator Variants: *hbond/dreiding/lj/angleoffset/omp*

pair_style hbond/dreiding/morse/angleoffset command
=======================================

Accelerator Variants: *hbond/dreiding/morse/angleoffset/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style N inner_distance_cutoff outer_distance_cutoff angle_cutoff equilibrium_angle

* style = *hbond/dreiding/lj/angleoffset* or *hbond/dreiding/morse/angleoffset*
* N = power of cosine of sum of angle theta and angle offset (integer)
* inner_distance_cutoff = global inner cutoff for Donor-Acceptor interactions (distance units)
* outer_distance_cutoff = global cutoff for Donor-Acceptor interactions (distance units)
* angle_cutoff = global angle cutoff for Acceptor-Hydrogen-Donor interactions (degrees)
* equilibrium_angle = global equilibrium angle for Acceptor-Hydrogen-Donor interactions (degrees)

(Tips: angle offset is the supplementary angle of equilibrium angle. It means if equilibrium angle is 166.6, the angle offset is 13.4)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style hybrid/overlay lj/cut 10.0 hbond/dreiding/lj/angleoffset 4 9.0 11.0 90.0 170.0
   pair_coeff 1 2 hbond/dreiding/lj 3 i 9.5 2.75 4 9.0 11.0 90.0 160.0

   pair_style hybrid/overlay lj/cut 10.0 hbond/dreiding/morse/angleoffset 2 9.0 11.0 90.0 170.0
   pair_coeff 1 2 hbond/dreiding/morse 3 i 3.88 1.7241379 2.9 2 9.0 11.0 90.0 160.0

   labelmap atom 1 C 2 O 3 H
   pair_coeff C O hbond/dreiding/morse H i 3.88 1.7241379 2.9 2 9.0 11.0 90.0 160.0

Description
"""""""""""

The *hbond/dreiding/\*/angleoffset* styles are modified version of *hbond/dreiding* styles.

In some cases, the angle of acceptor-hydrogen-donor in the equilibrium state could not achieve 180 degree especially in some coarse grained models.
In these cases, an angle offset is required to ensure that equilibrium state could be the minimum energy state.

The *hbond/dreiding/\*/angleoffset* styles compute the Acceptor-Hydrogen-Donor (AHD)
3-body hydrogen bond interaction for the :doc:`DREIDING <Howto_bioFF>`
force field with an angle offset, given by:

.. math::

   E  = & \left[LJ(r) | Morse(r) \right] \qquad \qquad \qquad r < r_{\rm in} and  \theta + \theta_{offset} < \theta_c \\
      = & S(r) * \left[LJ(r) | Morse(r) \right] \qquad \qquad r_{\rm in} < r < r_{\rm out} and \theta + \theta_{offset} < \theta_c \\
      = & 0 \qquad \qquad \qquad \qquad \qquad \qquad \qquad r > r_{\rm out} and \theta + \theta_{offset} < \theta_c \\
   LJ(r)  = & AR^{-12}-BR^{-10}cos^n(\theta + \theta_{offset})=
         \epsilon\left\lbrace 5\left[ \frac{\sigma}{r}\right]^{12}-
         6\left[ \frac{\sigma}{r}\right]^{10}  \right\rbrace cos^n(\theta + \theta_{offset})\\
   Morse(r)  = & D_0\left\lbrace \chi^2 - 2\chi\right\rbrace cos^n(\theta + \theta_{offset})=
         D_{0}\left\lbrace e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)}
         \right\rbrace cos^n(\theta + \theta_{offset})\
   S(r)  = & \frac{ \left[r_{\rm out}^2 - r^2\right]^2
   \left[r_{\rm out}^2 + 2r^2 - 3{r_{\rm in}^2}\right]}
   { \left[r_{\rm out}^2 - {r_{\rm in}}^2\right]^3 }

where :math:`r_{\rm in}` is the inner spline distance cutoff,
:math:`r_{\rm out}` is the outer distance cutoff, :math:`\theta_c` is
the angle cutoff, :math:`\theta_offset` is the angle offset, and :math:`n` is the power of the cosine of the sum of the angle :math:`\theta` and the angle offset :math:`\theta_offset`

Here, *r* is the radial distance between the donor (D) and acceptor
(A) atoms and :math:`\theta` is the bond angle between the acceptor, the
hydrogen (H) and the donor atoms:

.. image:: JPG/dreiding_hbond.jpg
   :align: center

For the *hbond/dreiding/lj/angleoffset* style the list of coefficients is as
follows:

* K = hydrogen atom type = 1 to Ntypes, or type label
* donor flag = *i* or *j*
* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* *n* = exponent in formula above
* distance cutoff :math:`r_{\rm in}` (distance units)
* distance cutoff :math:`r_{\rm out}` (distance units)
* angle cutoff (degrees)
* equilibrium angle (degrees)

(Tips: angle offset is the supplementary angle of equilibrium angle)

For the *hbond/dreiding/morse/angleoffset* style the list of coefficients is as
follows:

* K = hydrogen atom type = 1 to Ntypes, or type label
* donor flag = *i* or *j*
* :math:`D_0` (energy units)
* :math:`\alpha` (1/distance units)
* :math:`r_0` (distance units)
* *n* = exponent in formula above
* distance cutoff :math:`r_{\rm in}` (distance units)
* distance cutoff :math:`r_{out}` (distance units)
* angle cutoff (degrees)
* equilibrium angle (degrees)

(Tips: angle offset is the supplementary angle of equilibrium angle)

----------

Additional Information
""""""""""""

For more information about DREIDING force field and other notes, please refer to the documentation of *hbond/dreiding* styles.

----------

Restrictions
""""""""""""

This pair style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

