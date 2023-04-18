.. index:: pair_style lj/cut
.. index:: pair_style lj/cut/gpu
.. index:: pair_style lj/cut/intel
.. index:: pair_style lj/cut/kk
.. index:: pair_style lj/cut/opt
.. index:: pair_style lj/cut/omp

pair_style lj/cut command
=========================

Accelerator Variants: *lj/cut/gpu*, *lj/cut/intel*, *lj/cut/kk*, *lj/cut/opt*, *lj/cut/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/cut*
* args = list of arguments for a particular style

.. parsed-literal::

     *lj/cut* args = cutoff
       cutoff = global cutoff for Lennard Jones interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut 2.5
   pair_coeff * * 1 1
   pair_coeff 1 1 1 1.1 2.8

Description
"""""""""""

The *lj/cut* styles compute the standard 12/6 Lennard-Jones potential,
given by

.. math::

   E = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
       \left(\frac{\sigma}{r}\right)^6 \right]
                       \qquad r < r_c

:math:`r_c` is the cutoff.

See the :doc:`lj/cut/coul <pair_lj_cut_coul>` styles to add a Coulombic
pairwise interaction and the :doc:`lj/cut/tip4p <pair_lj_cut_tip4p>` styles to
add the TIP4P water model.

Coefficients
""""""""""""

The following coefficients must be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data
file or restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, or by mixing as described below:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* LJ cutoff (distance units)

Note that :math:`\sigma` is defined in the LJ formula as the zero-crossing
distance for the potential, not as the energy minimum at :math:`2^{\frac{1}{6}} \sigma`.

The last coefficient is optional.  If not specified, the global
LJ cutoff specified in the pair_style command is used.

----------

A version of these styles with a soft core, *lj/cut/soft*, suitable
for use in free energy calculations, is part of the FEP package and
is documented with the :doc:`pair_style */soft <pair_fep_soft>`
styles.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/cut pair styles can be mixed.
The default mix value is *geometric*.  See the "pair_modify" command
for details.

All of the *lj/cut* pair styles support the
:doc:`pair_modify <pair_modify>` shift option for the energy of the
Lennard-Jones portion of the pair interaction.

All of the *lj/cut* pair styles support the
:doc:`pair_modify <pair_modify>` tail option for adding a long-range
tail correction to the energy and pressure for the Lennard-Jones
portion of the pair interaction.

All of the *lj/cut* pair styles write their information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do
not need to be specified in an input script that reads a restart file.

The *lj/cut* pair styles support the use of the *inner*, *middle*, and
*outer* keywords of the :doc:`run_style respa <run_style>` command,
meaning the pairwise forces can be partitioned by distance at different
levels of the rRESPA hierarchy.  The other styles only support the
*pair* keyword of run_style respa.  See the :doc:`run_style <run_style>`
command for details.

----------

Related commands
""""""""""""""""

* :doc:`pair_coeff <pair_coeff>`
* :doc:`pair_style lj/cut/coul/cut <pair_lj_cut_coul>`
* :doc:`pair_style lj/cut/coul/debye <pair_lj_cut_coul>`
* :doc:`pair_style lj/cut/coul/dsf <pair_lj_cut_coul>`
* :doc:`pair_style lj/cut/coul/long <pair_lj_cut_coul>`
* :doc:`pair_style lj/cut/coul/msm <pair_lj_cut_coul>`
* :doc:`pair_style lj/cut/coul/wolf <pair_lj_cut_coul>`
* :doc:`pair_style lj/cut/tip4p/cut <pair_lj_cut_tip4p>`
* :doc:`pair_style lj/cut/tip4p/long <pair_lj_cut_tip4p>`

Default
"""""""

none
