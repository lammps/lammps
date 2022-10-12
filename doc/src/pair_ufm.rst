.. index:: pair_style ufm
.. index:: pair_style ufm/gpu
.. index:: pair_style ufm/omp
.. index:: pair_style ufm/opt

pair_style ufm command
======================

Accelerator Variants: *ufm/gpu*, *ufm/omp*, *ufm/opt*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style ufm cutoff

* cutoff = global cutoff for *ufm* interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style ufm 4.0
   pair_coeff 1 1 100.0 1.0 2.5
   pair_coeff * * 100.0 1.0

   pair_style ufm 4.0
   pair_coeff * * 10.0 1.0
   variable prefactor equal ramp(10,100)
   fix 1 all adapt 1 pair ufm epsilon * * v_prefactor

Description
"""""""""""

Style *ufm* computes pairwise interactions using the Uhlenbeck-Ford model (UFM) potential :ref:`(Paula Leite2016) <PL2>` which is given by

.. math::

   E & = -\varepsilon\, \ln{\left[1-\exp{\left(-r^{2}/\sigma^{2}\right)}\right]} \qquad  r < r_c \\
   \varepsilon & = p\,k_B\,T

where :math:`r_c` is the cutoff, :math:`\sigma` is a distance-scale and
:math:`\epsilon` is an energy-scale, i.e., a product of Boltzmann constant
:math:`k_B`, temperature :math:`T` and the Uhlenbeck-Ford p-parameter which
is responsible
to control the softness of the interactions :ref:`(Paula Leite2017) <PL1>`.
This model is useful as a reference system for fluid-phase free-energy calculations :ref:`(Paula Leite2016) <PL2>`.

The following coefficients must be defined for each pair of atom types
via the :doc:`pair_coeff <pair_coeff>` command as in the examples above,
or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global *ufm*
cutoff is used.

The :doc:`fix adapt <fix_adapt>` command can be used to vary epsilon and sigma for this pair style over the course of a simulation, in which case
pair_coeff settings for epsilon and sigma must still be specified, but will be
overridden.  For example these commands will vary the prefactor epsilon for
all pairwise interactions from 10.0 at the beginning to 100.0 at the end
of a run:

.. code-block:: LAMMPS

   variable prefactor equal ramp(10,100)
   fix 1 all adapt 1 pair ufm epsilon * * v_prefactor

.. note::

   The thermodynamic integration procedure can be performed with this
   potential using :doc:`fix adapt <fix_adapt>`. This command will
   rescale the force on each atom by varying a scale variable, which
   always starts with value 1.0. The syntax is the same described above,
   however, changing epsilon to scale. A detailed explanation of how to
   use this command and perform nonequilibrium thermodynamic integration
   in LAMMPS is given in the paper by :ref:`(Freitas) <Freitas2>`.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients and cutoff
distance for this pair style can be mixed. The default mix value is *geometric*\ .  See the
"pair_modify" command for details.

This pair style support the :doc:`pair_modify <pair_modify>` shift option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table and tail are not relevant for this
pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>` tail option for adding long-range tail corrections to energy and pressure.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the EXTRA-PAIR package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`fix adapt <fix_adapt>`

Default
"""""""

none

.. _PL1:

**(Paula Leite2017)** Paula Leite, Santos-Florez, and de Koning, Phys Rev E, 96,
32115 (2017).

.. _PL2:

**(Paula Leite2016)** Paula Leite , Freitas, Azevedo, and de Koning, J Chem Phys, 126,
044509 (2016).

.. _Freitas2:

**(Freitas)** Freitas, Asta, and de Koning, Computational Materials Science, 112, 333 (2016).
