.. index:: pair_style lj/switch3/coulgauss/long
.. index:: pair_style mm3/switch3/coulgauss/long

pair_style lj/switch3/coulgauss/long command
============================================

pair_style mm3/switch3/coulgauss/long command
=============================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/switch3/coulgauss/long* or *mm3/switch3/coulgauss/long*
* args = list of arguments for a particular style

.. parsed-literal::

     *lj/switch3/coulgauss/long* args = cutoff (cutoff2) width
       cutoff  = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
       width  = width parameter of the smoothing function (distance units)

     *mm3/switch3/coulgauss/long* args = cutoff (cutoff2) width
       cutoff  = global cutoff for MM3 (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
       width  = width parameter of the smoothing function (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/switch3/coulgauss/long    12.0 3.0
   pair_coeff 1  0.2 2.5 1.2

   pair_style lj/switch3/coulgauss/long   12.0 10.0 3.0
   pair_coeff 1  0.2 2.5 1.2

   pair_style mm3/switch3/coulgauss/long    12.0 3.0
   pair_coeff 1  0.2 2.5 1.2

   pair_style mm3/switch3/coulgauss/long   12.0 10.0 3.0
   pair_coeff 1  0.2 2.5 1.2

Description
"""""""""""

The *lj/switch3/coulgauss* style evaluates the LJ
vdW potential

.. math::

   E = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^{6} \right]

The *mm3/switch3/coulgauss/long* style evaluates the MM3
vdW potential :ref:`(Allinger) <mm3-allinger1989>`

.. math::

   E & = \epsilon_{ij} \left[ -2.25 \left(\frac{r_{v,ij}}{r_{ij}}\right)^6 + 1.84(10)^5 \exp\left[-12.0 r_{ij}/r_{v,ij}\right] \right] S_3(r_{ij}) \\
   r_{v,ij} & =  r_{v,i} + r_{v,j} \\
   \epsilon_{ij} & = \sqrt{\epsilon_i \epsilon_j}

Both potentials go smoothly to zero at the cutoff r_c as defined by the
switching function

.. math::

 S_3(r) = \left\lbrace \begin{array}{ll}
                     1 & \quad\mathrm{if}\quad r < r_\mathrm{c} - w \\
                     3x^2 - 2x^3 & \quad\mathrm{if}\quad r < r_\mathrm{c} \quad\mathrm{with\quad} x=\frac{r_\mathrm{c} - r}{w} \\
                     0 & \quad\mathrm{if}\quad r >= r_\mathrm{c}
                 \end{array} \right.

where w is the width defined in the arguments. This potential
is combined with Coulomb interaction between Gaussian charge densities:

.. math::

   E = \frac{q_i q_j \mathrm{erf}\left( r/\sqrt{\gamma_1^2+\gamma_2^2} \right) }{\epsilon r_{ij}}

where :math:`q_i` and :math:`q_j` are the charges on the 2 atoms,
:math:`\epsilon` is the dielectric constant which can be set by the
:doc:`dielectric <dielectric>` command, :math:`\gamma_i` and
:math:`\gamma_j` are the widths of the Gaussian charge distribution and
erf() is the error-function.  This style has to be used in conjunction
with the :doc:`kspace_style <kspace_style>` command

If one cutoff is specified it is used for both the vdW and Coulomb
terms.  If two cutoffs are specified, the first is used as the cutoff
for the vdW terms, and the second is the cutoff for the Coulombic term.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`\epsilon` (energy)
* :math:`\sigma` (distance)
* :math:`\gamma` (distance)

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/long pair styles can be mixed.
The default mix value is *geometric*\ .  See the "pair_modify" command
for details.

Shifting the potential energy is not necessary because the switching
function ensures that the potential is zero at the cut-off.

These pair styles support the :doc:`pair_modify <pair_modify>` table and
options since they can tabulate the short-range portion of the
long-range Coulombic interactions.

Thes pair styles do not support the :doc:`pair_modify <pair_modify>`
tail option for adding a long-range tail correction to the
Lennard-Jones portion of the energy and pressure.

These pair styles write their information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*, *middle*, *outer* keywords.

Restrictions
""""""""""""

These styles are part of the YAFF package.  They are only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none
