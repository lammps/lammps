.. index:: pair_style lj/cut/coul/long/gauss

pair_style lj/cut/coul/long/gauss command
=========================================

.. index:: pair_style lj/cut/coul/wolf/gauss

pair_style lj/cut/coul/wolf/gauss command
=========================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/cut/coul/long/gauss* or *lj/cut/coul/wolf/gauss*
* args = list of arguments for a particular style

.. parsed-literal::

     *lj/cut/coul/long/gauss* args = cutoff (cutoff2)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *lj/cut/coul/wolf/gauss* args = alpha cutoff (cutoff2)
       alpha = damping parameter (inverse distance units)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/long/gauss 10.0
   pair_style lj/cut/coul/long/gauss 10.0 8.0
   pair_style lj/cut/coul/wolf/gauss 0.2 10.0
   pair_coeff 1 1 1.13  5.06 NULL
   pair_coeff 2 2 0.055 3.37 1.979

Description
"""""""""""

The style computes the standard 12/6 Lennard-Jones potential, given by

.. math::

   E = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
       \left(\frac{\sigma}{r}\right)^6 \right]
                       \qquad r < r_c

:math:`r_c` is the cutoff.

The Coulomb interaction of atoms with a Gaussian charge width

.. math::

   \rho_i(\mathbf r) = \left(\frac{\eta_i}{\sqrt{\pi}}\right)^3 q_i
      \mbox{exp} \left( -\eta_i^2 | \mathbf r - \mathbf r_i | ^2 \right)

or point charges is given by

.. math::

   E = q_i q_j \frac{\mbox{erf} (\eta_{ij} r_{ij})}{r_{ij}}

and the self-interaction of Gaussian charges is given by

.. math::

   E = \frac{\eta_i}{\sqrt{2\pi}} q_i^2.

The pair style *lj/cut/coul/long/gauss* calculates the short-range term of the
Coulomb energy as the energy of point charges with a correction for the Gaussian
charge width as derived by :ref:`Gingrich and Wilson <GingrichWilson>`. The
short-range term has to be used in combination with a :doc:`kspace_style
<kspace_style>`.

The pair style *lj/cut/coul/wolf/gauss* computes the Coulomb energy via the
:ref:`Wolf <Wolf6>` summation method with a correction for Gaussian charges.

.. note::

   The inequality :math:`\eta_i > \sqrt{2} \alpha` needs to hold (:ref:`Gingrich
   and Wilson <GingrichWilson>`). :math:`\alpha` is the damping parameter and is
   either gewald (cf. :doc:`kspace_modify <kspace_modify>`) or alpha in the Wolf
   summation. The damping parameter needs to be sufficiently small relative to
   the reciprocal width of Gaussian charges. In particluar, the inequality is
   required to ensure the matrix in :doc:`fix electrode <fix_electrode>` is
   positive-definite.

Coefficients
""""""""""""

For every atom type :math:`i`, the :math:`\eta_i` parameter in 1/distance units
must be defined via the :doc:`pair_coeff <pair_coeff>` command as in the
examples above or by the :doc:`read_data <read_data>` or :doc:`read_restart
<read_restart>` commands.  For the point charges the value has to be defined as
"NULL" instead of a real value. The mixed parameters are

.. math::

   \eta_{ij} = \frac{\eta_i \eta_j}{\sqrt{\eta_i^2 + \eta_j^2}}

for two Gaussian charges and :math:`\eta_{ij} = \eta_i` if  atoms of type
:math:`i` have Gaussian charges and atoms of type :math:`j` are point charges
and vice versa.

The following coefficients must be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data
file or restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, or by mixing as described for
:doc:`pair_lj_cut_coul <pair_lj_cut_coul>`.

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* :math:`\eta_i` (1/distance units)
* cutoff1 (distance units)
* cutoff2 (distance units)

Restrictions
""""""""""""

*lj/cut/coul/long/gauss* and *lj/cut/coul/wolf/gauss* do not support the
run_style respa.

These pair styles are part of the ELECTRODE package. They are only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix electrode <fix_electrode>`

Default
"""""""

none

----------

.. _GingrichWilson:

**(Gingrich and Wilson)** Gingrich and Wilson, Chem. Phys. Lett., 500, 178-183
(2010).

.. _Wolf6:

**(Wolf)** D. Wolf, P. Keblinski, S. R. Phillpot, J. Eggebrecht, J Chem
Phys, 110, 8254 (1999).
