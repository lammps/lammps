.. index:: pair_style edip

pair_style edip command
=======================

pair_style edip/omp command
===========================

pair_style edip/multi command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style

* style = *edip* or *edip/multi*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style edip
   pair_coeff * * Si.edip Si

Description
"""""""""""

The *edip* and *edip/multi* styles compute a 3-body :ref:`EDIP <EDIP>`
potential which is popular for modeling silicon materials where
it can have advantages over other models such as the
:doc:`Stillinger-Weber <pair_sw>` or :doc:`Tersoff <pair_tersoff>`
potentials. The *edip* style has been programmed for single element
potentials, while *edip/multi* supports multi-element EDIP runs.

In EDIP, the energy E of a system of atoms is

.. math::

   E  = & \sum_{j \ne i} \phi_{2}(R_{ij}, Z_{i}) + \sum_{j \ne i} \sum_{k \ne i,k > j} \phi_{3}(R_{ij}, R_{ik}, Z_{i}) \\
   \phi_{2}(r, Z)  = & A\left[\left(\frac{B}{r}\right)^{\rho} - e^{-\beta Z^2}\right]exp{\left(\frac{\sigma}{r-a}\right)} \\
   \phi_{3}(R_{ij}, R_{ik}, Z_i)  = & exp{\left(\frac{\gamma}{R_{ij}-a}\right)}exp{\left(\frac{\gamma}{R_{ik}-a}\right)}h(cos\theta_{ijk},Z_i) \\
   Z_i  = & \sum_{m \ne i} f(R_{im}) \qquad
   f(r) = \begin{cases}
          1 & \quad r<c \\
          \exp\left(\frac{\alpha}{1-x^{-3}}\right) & \quad c<r<a \\
          0 & \quad r>a
          \end{cases} \\
   h(l,Z)  = & \lambda [(1-e^{-Q(Z)(l+\tau(Z))^2}) + \eta Q(Z)(l+\tau(Z))^2 ] \\
   Q(Z)  = & Q_0 e^{-\mu Z} \qquad \tau(Z) = u_1 + u_2 (u_3 e^{-u_4 Z} - e^{-2u_4 Z})

where :math:`\phi_2` is a two-body term and :math:`\phi_3` is a
three-body term.  The summations in the formula are over all neighbors J
and K of atom I within a cutoff distance = a.  Both terms depend on the
local environment of atom I through its effective coordination number
defined by Z, which is unity for a cutoff distance < c and gently goes
to 0 at distance = a.

Only a single pair_coeff command is used with the *edip* style which
specifies a EDIP potential file with parameters for all
needed elements.  These are mapped to LAMMPS atom types by specifying
N additional arguments after the filename in the pair_coeff command,
where N is the number of LAMMPS atom types:

* filename
* N element names = mapping of EDIP elements to atom types

See the :doc:`pair_coeff <pair_coeff>` doc page for alternate ways
to specify the path for the potential file.

As an example, imagine a file Si.edip has EDIP values for Si.

EDIP files in the *potentials* directory of the LAMMPS
distribution have a ".edip" suffix.  Lines that are not blank or
comments (starting with #) define parameters for a triplet of
elements.  The parameters in a single entry correspond to the two-body
and three-body coefficients in the formula above:

* element 1 (the center atom in a 3-body interaction)
* element 2
* element 3
* A (energy units)
* B (distance units)
* cutoffA (distance units)
* cutoffC (distance units)
* :math:`\alpha`
* :math:`\beta`
* :math:`\eta`
* :math:`\gamma` (distance units)
* :math:`lambda` (energy units)
* :math:`\mu`
* :math:`\tau`
* :math:`\sigma` (distance units)
* Q0
* u1
* u2
* u3
* u4

The A, B, beta, sigma parameters are used only for two-body interactions.
The eta, gamma, lambda, mu, Q0 and all u1 to u4 parameters are used only
for three-body interactions. The alpha and cutoffC parameters are used
for the coordination environment function only.

The EDIP potential file must contain entries for all the
elements listed in the pair_coeff command.  It can also contain
entries for additional elements not being used in a particular
simulation; LAMMPS ignores those entries.

For a single-element simulation, only a single entry is required
(e.g. SiSiSi).  For a two-element simulation, the file must contain 8
entries (for SiSiSi, SiSiC, SiCSi, SiCC, CSiSi, CSiC, CCSi, CCC), that
specify EDIP parameters for all permutations of the two elements
interacting in three-body configurations.  Thus for 3 elements, 27
entries would be required, etc.

At the moment, only a single element parameterization is
implemented. However, the author is not aware of other
multi-element EDIP parameterization. If you know any and
you are interest in that, please contact the author of
the EDIP package.

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

**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.

----------

Restrictions
""""""""""""

This pair style can only be used if LAMMPS was built with the
USER-MISC package.  See the :doc:`Build package <Build_package>` doc
page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The EDIP potential files provided with LAMMPS (see the potentials directory)
are parameterized for metal :doc:`units <units>`.
You can use the EDIP potential with any LAMMPS units, but you would need
to create your own EDIP potential file with coefficients listed in the
appropriate units if your simulation does not use "metal" units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none

----------

.. _EDIP:

**(EDIP)** J F Justo et al, Phys Rev B 58, 2539 (1998).
