.. index:: pair_style lj/sphere
.. index:: pair_style lj/sphere/omp

pair_style lj/sphere command
============================

Accelerator Variant: *lj/sphere/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/sphere*
* args = list of arguments for a particular style

.. parsed-literal::

     *lj/sphere* args = cutoff
       cutoff = global cutoff for Lennard Jones interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/sphere 2.5
   pair_coeff * * 1.0
   pair_coeff 1 1 1.1 2.8

Description
"""""""""""

The *lj/sphere* style compute the standard 12/6 Lennard-Jones potential,
given by

.. math::

   E = 4 \epsilon \left[ \left(\frac{\sigma_{ij}}{r}\right)^{12} -
       \left(\frac{\sigma_{ij}}{r}\right)^6 \right]
                       \qquad r < r_c

:math:`r_c` is the cutoff.

This is the same potential function as used by the :doc:`lj/cut
<pair_lj>` pair style, but the :math:`\sigma_{ij}` parameter is not set
as a per-type parameter via the :doc:`pair_coeff command <pair_coeff>`,
but taken from the per-atom diameter attribute of :doc:`atom_style
sphere <atom_style>`.  The individual value of :math:`\sigma_{ij}` is
computed from the diameters of the two atoms using the mixing rule for
pair coefficients as set by the :doc:`pair_modify mix <pair_modify>`
command (defaults to geometric mixing).

Note that :math:`\sigma_{ij}` is defined in the LJ formula above as the
zero-crossing distance for the potential, *not* as the energy minimum which
is at :math:`2^{\frac{1}{6}} \sigma_{ij}`.

Coefficients
""""""""""""

The following coefficients must be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data
file or restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, or by mixing as described below:

* :math:`\epsilon` (energy units)
* LJ cutoff (distance units) (optional)

The last coefficient is optional.  If not specified, the global
LJ cutoff specified in the pair_style command is used.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon coefficients and cutoff
distance for the *lj/sphere* pair style can be mixed.  The default mix value
is *geometric*.  See the "pair_modify" command for details.

The *lj/sphere* pair style supports the :doc:`pair_modify <pair_modify>`
shift option for the energy of the Lennard-Jones portion of the pair
interaction.

The *lj/sphere* pair style does *not* support the :doc:`pair_modify
<pair_modify>` tail option for adding a long-range tail corrections to
the energy and pressure.

The *lj/sphere* pair style writes its information to :doc:`binary
restart files <restart>`, so pair_style and pair_coeff commands do not
need to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does *not* support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

The *lj/sphere* pair style is only enabled if LAMMPS was built with the
EXTRA-PAIR package.  See the :doc:`Build package <Build_package>` page
for more info.

The *lj/sphere* pair style does not support the *sixthpower* mixing rule.

----------

Related commands
""""""""""""""""

* :doc:`pair_coeff <pair_coeff>`
* :doc:`pair_style lj/cut <pair_lj>`

Default
"""""""

none
