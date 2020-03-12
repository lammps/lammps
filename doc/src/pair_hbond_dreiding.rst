.. index:: pair_style hbond/dreiding/lj

pair_style hbond/dreiding/lj command
====================================

pair_style hbond/dreiding/lj/omp command
========================================

pair_style hbond/dreiding/morse command
=======================================

pair_style hbond/dreiding/morse/omp command
===========================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style N inner_distance_cutoff outer_distance_cutoff angle_cutof

* style = *hbond/dreiding/lj* or *hbond/dreiding/morse*
* n = cosine angle periodicity
* inner\_distance\_cutoff = global inner cutoff for Donor-Acceptor interactions (distance units)
* outer\_distance\_cutoff = global cutoff for Donor-Acceptor interactions (distance units)
* angle\_cutoff = global angle cutoff for Acceptor-Hydrogen-Donor
* interactions (degrees)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style hybrid/overlay lj/cut 10.0 hbond/dreiding/lj 4 9.0 11.0 90
   pair_coeff 1 2 hbond/dreiding/lj 3 i 9.5 2.75 4 9.0 11.0 90.0

   pair_style hybrid/overlay lj/cut 10.0 hbond/dreiding/morse 2 9.0 11.0 90
   pair_coeff 1 2 hbond/dreiding/morse 3 i 3.88 1.7241379 2.9 2 9 11 90

Description
"""""""""""

The *hbond/dreiding* styles compute the Acceptor-Hydrogen-Donor (AHD)
3-body hydrogen bond interaction for the :doc:`DREIDING <Howto_bioFF>`
force field, given by:

.. math::

   E  = & \left[LJ(r) | Morse(r) \right] \qquad \qquad \qquad r < r_{\rm in} \\
      = & S(r) * \left[LJ(r) | Morse(r) \right] \qquad \qquad r_{\rm in} < r < r_{\rm out} \\
      = & 0 \qquad \qquad \qquad \qquad \qquad \qquad \qquad r > r_{\rm out} \\
   LJ(r)  = & AR^{-12}-BR^{-10}cos^n\theta=
         \epsilon\left\lbrace 5\left[ \frac{\sigma}{r}\right]^{12}-
         6\left[ \frac{\sigma}{r}\right]^{10}  \right\rbrace cos^n\theta\\
   Morse(r)  = & D_0\left\lbrace \chi^2 - 2\chi\right\rbrace cos^n\theta=
         D_{0}\left\lbrace e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)}
         \right\rbrace cos^n\theta \\
   S(r)  = & \frac{ \left[r_{\rm out}^2 - r^2\right]^2
   \left[r_{\rm out}^2 + 2r^2 - 3{r_{\rm in}^2}\right]}
   { \left[r_{\rm out}^2 - {r_{\rm in}}^2\right]^3 }

where :math:`r_{\rm in}` is the inner spline distance cutoff,
:math:`r_{\rm out}` is the outer distance cutoff, :math:`\theta_c` is
the angle cutoff, and *n* is the cosine periodicity.

Here, *r* is the radial distance between the donor (D) and acceptor
(A) atoms and :math:`\theta` is the bond angle between the acceptor, the
hydrogen (H) and the donor atoms:

.. image:: JPG/dreiding_hbond.jpg
   :align: center

These 3-body interactions can be defined for pairs of acceptor and
donor atoms, based on atom types.  For each donor/acceptor atom pair,
the 3rd atom in the interaction is a hydrogen permanently bonded to
the donor atom, e.g. in a bond list read in from a data file via the
:doc:`read_data <read_data>` command.  The atom types of possible
hydrogen atoms for each donor/acceptor type pair are specified by the
:doc:`pair_coeff <pair_coeff>` command (see below).

Style *hbond/dreiding/lj* is the original DREIDING potential of
:ref:`(Mayo) <pair-Mayo>`.  It uses a LJ 12/10 functional for the Donor-Acceptor
interactions. To match the results in the original paper, use n = 4.

Style *hbond/dreiding/morse* is an improved version using a Morse
potential for the Donor-Acceptor interactions. :ref:`(Liu) <Liu>` showed
that the Morse form gives improved results for Dendrimer simulations,
when n = 2.

See the :doc:`Howto bioFF <Howto_bioFF>` doc page for more information
on the DREIDING force field.

.. note::

   Because the Dreiding hydrogen bond potential is only one portion
   of an overall force field which typically includes other pairwise
   interactions, it is common to use it as a sub-style in a :doc:`pair_style hybrid/overlay <pair_hybrid>` command, where another pair style
   provides the repulsive core interaction between pairs of atoms, e.g. a
   1/r\^12 Lennard-Jones repulsion.

.. note::

   When using the hbond/dreiding pair styles with :doc:`pair_style hybrid/overlay <pair_hybrid>`, you should explicitly define pair
   interactions between the donor atom and acceptor atoms, (as well as
   between these atoms and ALL other atoms in your system).  Whenever
   :doc:`pair_style hybrid/overlay <pair_hybrid>` is used, ordinary mixing
   rules are not applied to atoms like the donor and acceptor atoms
   because they are typically referenced in multiple pair styles.
   Neglecting to do this can cause difficult-to-detect physics problems.

.. note::

   In the original Dreiding force field paper 1-4 non-bonded
   interactions ARE allowed.  If this is desired for your model, use the
   special\_bonds command (e.g. "special\_bonds lj 0.0 0.0 1.0") to turn
   these interactions on.

----------

The following coefficients must be defined for pairs of eligible
donor/acceptor types via the :doc:`pair_coeff <pair_coeff>` command as
in the examples above.

.. note::

   Unlike other pair styles and their associated
   :doc:`pair_coeff <pair_coeff>` commands, you do not need to specify
   pair\_coeff settings for all possible I,J type pairs.  Only I,J type
   pairs for atoms which act as joint donors/acceptors need to be
   specified; all other type pairs are assumed to be inactive.

.. note::

   A :doc:`pair_coeff <pair_coeff>` command can be specified multiple
   times for the same donor/acceptor type pair.  This enables multiple
   hydrogen types to be assigned to the same donor/acceptor type pair.
   For other pair\_styles, if the pair\_coeff command is re-used for the
   same I.J type pair, the settings for that type pair are overwritten.
   For the hydrogen bond potentials this is not the case; the settings
   are cumulative.  This means the only way to turn off a previous
   setting, is to re-use the pair\_style command and start over.

For the *hbond/dreiding/lj* style the list of coefficients is as
follows:

* K = hydrogen atom type = 1 to Ntypes
* donor flag = *i* or *j*
* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* *n* = exponent in formula above
* distance cutoff :math:`r_{\rm in}` (distance units)
* distance cutoff :math:`r_{\rm out}` (distance units)
* angle cutoff (degrees)

For the *hbond/dreiding/morse* style the list of coefficients is as
follows:

* K = hydrogen atom type = 1 to Ntypes
* donor flag = *i* or *j*
* :math:`D_0` (energy units)
* :math:`\alpha` (1/distance units)
* :math:`r_0` (distance units)
* *n* = exponent in formula above
* distance cutoff :math:`r_{\rm in}` (distance units)
* distance cutoff :math:`r_{out}` (distance units)
* angle cutoff (degrees)

A single hydrogen atom type K can be specified, or a wild-card asterisk
can be used in place of or in conjunction with the K arguments to
select multiple types as hydrogen atoms.  This takes the form
"\*" or "\*n" or "n\*" or "m\*n".  See the :doc:`pair_coeff <pair_coeff>`
command doc page for details.

If the donor flag is *i*\ , then the atom of type I in the pair\_coeff
command is treated as the donor, and J is the acceptor.  If the donor
flag is *j*\ , then the atom of type J in the pair\_coeff command is
treated as the donor and I is the donor.  This option is required
because the :doc:`pair_coeff <pair_coeff>` command requires that I <= J.

:math:`\epsilon` and :math:`\sigma` are settings for the hydrogen bond
potential based on a Lennard-Jones functional form.  Note that sigma is
defined as the zero-crossing distance for the potential, not as the
energy minimum at :math:`2^{1/6} \sigma`.

:math:`D_0` and :math:`\alpha` and :math:`r_0` are settings for the
hydrogen bond potential based on a Morse functional form.

The last 3 coefficients for both styles are optional.  If not
specified, the global n, distance cutoff, and angle cutoff specified
in the pair\_style command are used.  If you wish to only override the
2nd or 3rd optional parameter, you must also specify the preceding
optional parameters.

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

These pair styles do not support mixing. You must explicitly identify
each donor/acceptor type pair.

These styles do not support the :doc:`pair_modify <pair_modify>` shift
option for the energy of the interactions.

The :doc:`pair_modify <pair_modify>` table option is not relevant for
these pair styles.

These pair styles do not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

These pair styles do not write their information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands need to be
re-specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*\ , *middle*\ , *outer* keywords.

These pair styles tally a count of how many hydrogen bonding
interactions they calculate each timestep and the hbond energy.  These
quantities can be accessed via the :doc:`compute pair <compute_pair>`
command as a vector of values of length 2.

To print these quantities to the log file (with a descriptive column
heading) the following commands could be included in an input script:

.. code-block:: LAMMPS

   compute hb all pair hbond/dreiding/lj
   variable n_hbond equal c_hb[1] #number hbonds
   variable E_hbond equal c_hb[2] #hbond energy
   thermo_style custom step temp epair v_E_hbond

----------

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none

----------

.. _pair-Mayo:

**(Mayo)** Mayo, Olfason, Goddard III, J Phys Chem, 94, 8897-8909
(1990).

.. _Liu:

**(Liu)** Liu, Bryantsev, Diallo, Goddard III, J. Am. Chem. Soc 131 (8)
2798 (2009)
