.. index:: pair_style lj/relres
.. index:: pair_style lj/relres/omp

pair_style lj/relres command
============================

Accelerator Variants: *lj/relres/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/relres Rsi Rso Rci Rco

* Rsi = inner switching cutoff between the fine-grained and coarse-grained potentials (distance units)
* Rso = outer switching cutoff between the fine-grained and coarse-grained potentials (distance units)
* Rci = inner cutoff beyond which the force smoothing for all interactions is applied (distance units)
* Rco = outer cutoff for all interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/relres 4.0 5.0 8.0 10.0
   pair_coeff 1 1 0.5 1.0 1.5 1.1
   pair_coeff 2 2 0.5 1.0 0.0 0.0 3.0 3.5 6.0 7.0

Description
"""""""""""

Pair style *lj/relres* computes a LJ interaction using the Relative
Resolution (RelRes) framework which applies a fine-grained (FG)
potential between near neighbors and a coarse-grained (CG) potential
between far neighbors :ref:`(Chaimovich1) <Chaimovich1>`. This approach
can improve the computational efficiency by almost an order of
magnitude, while maintaining the correct static and dynamic behavior of
a reference system :ref:`(Chaimovich2) <Chaimovich2>`.

.. math::

   E = \left\{\begin{array}{lr}
        4 \epsilon^{\scriptscriptstyle FG} \left[ \left(\frac{\sigma^{FG}}{r}\right)^{12} - \left(\frac{\sigma^{FG}}{r}\right)^6 \right]-\Gamma_{si}, & \quad\mathrm{if}\quad  r< r_{si}, \\
        \sum_{m=0}^{4} \gamma_{sm}\left(r-r_{si}\right)^m-\Gamma_{so} ,   & \quad\mathrm{if}\quad  r_{si}\leq r< r_{so}, \\
        4 \epsilon^{\scriptscriptstyle CG} \left[ \left(\frac{\sigma^{CG}}{r}\right)^{12} -     \left(\frac{\sigma^{CG}}{r}\right)^6 \right]-\Gamma_c, &  \quad\mathrm{if}\quad  r_{so}\leq r<r_{ci}, \\
        \sum_{m=0}^{4} \gamma_{cm}\left(r-r_{ci}\right)^m -\Gamma_c, & \quad\mathrm{if}\quad  r_{ci}\leq r< r_{co}, \\
        0, & \quad\mathrm{if}\quad  r\geq r_{co}.\end{array}\right.

The FG parameters of the LJ potential (:math:`\epsilon^{FG}` and
:math:`\sigma^{FG}`) are applied up to the inner switching cutoff,
:math:`r_{si}`, while the CG parameters of the LJ potential
(:math:`\epsilon^{CG}` and :math:`\sigma^{CG}`) are applied beyond the
outer switching cutoff, :math:`r_{so}`. Between :math:`r_{si}` and
:math:`r_{so}` a polynomial smoothing function is applied so that the
force and its derivative are continuous between the FG and CG
potentials.  An analogous smoothing function is applied between the
inner and outer cutoffs (:math:`r_{ci}` and :math:`r_{co}`).
The offsets :math:`\Gamma_{si}`, :math:`\Gamma_{so}` and
:math:`\Gamma_{c}` ensure the continuity of the energy over the entire
domain.  The corresponding polynomial coefficients :math:`\gamma_{sm}`
and :math:`\gamma_{cm}`, as well as the offsets are automatically
computed by LAMMPS.

.. note::

   Energy and force resulting from this methodology can be plotted via the
   :doc:`pair_write <pair_write>` command.

The following coefficients must be defined for each pair of atom types
via the :doc:`pair_coeff <pair_coeff>` command as in the examples above,
or in the data file or restart files read by the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands, or by
mixing as will be described below:

* :math:`\epsilon^{FG}` (energy units)
* :math:`\sigma^{FG}` (distance units)
* :math:`\epsilon^{CG}` (energy units)
* :math:`\sigma^{CG}` (distance units)

Additional parameters can be defined to specify different
:math:`r_{si}`, :math:`r_{so}`, :math:`r_{ci}`, :math:`r_{co}` for
a particular set of atom types:

* :math:`r_{si}` (distance units)
* :math:`r_{so}` (distance units)
* :math:`r_{ci}` (distance units)
* :math:`r_{co}` (distance units)

These parameters are optional, and they are used to override the global
cutoffs as defined in the pair_style command.  If not specified, the
global values for :math:`r_{si}`, :math:`r_{so}`, :math:`r_{ci}`, and
:math:`r_{co}` are used.  If this override option is employed, all four
arguments must be specified.

----------

Here are some guidelines for using the pair_style *lj/relres* command.

In general, RelRes focuses on the speedup of pairwise interactions between
all LJ sites. Importantly, it works with any settings and flags (e.g.,
:doc:`special_bonds <special_bonds>` settings and :doc:`newton <newton>`
flags) that can be used in a molecular simulation with the
conventional LJ potential. In particular, all intramolecular topology
with its energetics (i.e., bonds, angles, etc.) remains unaltered.

At the most basic level in the RelRes framework, all sites are mapped into
clusters. Each cluster is just a collection of sites bonded together (the
bonds themselves are not part of the cluster). In general, a molecule may
be comprised of several clusters, and preferably, no two sites in a cluster
are separated by more than two bonds. There are two categories of sites in
RelRes: "hybrid" sites embody both FG and CG models, while "ordinary" sites
embody just FG characteristics with no CG features. A given cluster has
a single hybrid site (typically its central site) and several ordinary sites
(typically its peripheral sites). Notice that while clusters are necessary
for the RelRes parameterization (discussed below), they are not actually
defined in LAMMPS. Besides, the total number of sites in the cluster are
called the "mapping ratio", and this substantially impacts the computational
efficiency of RelRes: For a mapping ratio of 3, the efficiency factor is
around 4, and for a mapping ratio of 5, the efficiency factor is around 5
:ref:`(Chaimovich2) <Chaimovich2>`.

The flexibility of LAMMPS allows placing any values for the LJ
parameters in the input script. However, here are the optimal
recommendations for the RelRes parameters, which yield the correct
structural and thermal behavior in a system of interest
:ref:`(Chaimovich1) <Chaimovich1>`.  One must first assign a complete set of
parameters for the FG interactions that are applicable to all atom types.
Regarding the parameters for the CG interactions, the rules rely on the
site category (if it is a hybrid or an ordinary site). For atom types of
ordinary sites, :math:`\epsilon^{CG}` must be set to 0 (zero) while the
specific value of :math:`\sigma^{CG}` is irrelevant. For atom types of
hybrid sites, the CG parameters should be generally calculated using the
following equations:

.. math::

   \sigma_I^{CG}=\frac{\left((\sum_{\alpha\in A}\sqrt{\epsilon_\alpha^{FG}\left(\sigma_\alpha^{FG}\right)^{12}}\right)^{1/2}}{\left((\sum_{\alpha\in A}\sqrt{\epsilon_\alpha^{FG}\left(\sigma_\alpha^{FG}\right)^6}\right)^{1/3}}
   \quad\mathrm{and}\quad
   \epsilon_I^{CG}=\frac{\left((\sum_{\alpha\in A}\sqrt{\epsilon_\alpha^{FG}\left(\sigma_\alpha^{FG}\right)^6}\right)^4}{\left((\sum_{\alpha\in A}\sqrt{\epsilon_\alpha^{FG}\left(\sigma_\alpha^{FG}\right)^{12}}\right)^2}

where :math:`I` is an atom type of a hybrid site of a particular cluster
:math:`A`, and corresponding with this cluster, the summation proceeds over
all of its sites :math:`\alpha`.  These equations are derived from the
monopole term in the underlying Taylor series, and they are indeed relevant
only if geometric mixing is applicable for the FG model; if this is not the
case, Ref. :ref:`(Chaimovich2) <Chaimovich2>` discusses the alternative
formula, and in such a situation, the pair_coeff command should be explicitly
used for all combinations of atom types :math:`I\;!=J`.

The switching distance (the midpoint between inner and outer switching
cutoffs) is another crucial factor in RelRes: decreasing it improves the
computational efficiency, yet if it is too small, the molecular simulations
may not capture the system behavior correctly.  As a rule of thumb,
the switching distance should be approximately :math:`\,\sim\! 1.5\sigma`
:ref:`(Chaimovich1) <Chaimovich1>`; recommendations can be found in
Ref. :ref:`(Chaimovich2) <Chaimovich2>`.
Regarding the switching smoothing zone, :math:`\,\sim\!0.1\sigma` is
recommended; if desired, smoothing can be eliminated by setting
the inner switching cutoff, :math:`r_{si}`, equal to the outer
switching cutoff, :math:`r_{so}` (the same is true for the other cutoffs
:math:`r_{ci}` and :math:`r_{co}`).

----------

As an example, imagine that in your system, a molecule is comprised just
of one cluster such that one atom type (#1) is associated with
its hybrid site, and another atom type (#2) is associated with its ordinary
sites (in total, there are 2 atom types). If geometric mixing is applicable,
the following commands should be used:

.. code-block:: LAMMPS

   pair_style lj/relres Rsi Rso Rci Rco
   pair_coeff 1 1 epsilon_FG1 sigma_FG1 epsilon_CG1 sigma_CG1
   pair_coeff 2 2 epsilon_FG2 sigma_FG2 0.0         0.0
   pair_modify shift yes

In a more complex situation, there may be two distinct clusters in a system
(these two clusters may be on same molecule or on different molecules),
each with its own switching cutoffs.  If there are still two atom types
in each cluster as in the earlier example, the commands should be:

.. code-block:: LAMMPS

   pair_style lj/relres Rsi Rso Rci Rco
   pair_coeff 1 1 epsilon_FG1 sigma_FG1 epsilon_CG1 sigma_CG1 Rsi1 Rso1 Rci Rco
   pair_coeff 2 2 epsilon_FG2 sigma_FG2 0.0         0.0       Rsi1 Rso1 Rci Rco
   pair_coeff 3 3 epsilon_FG3 sigma_FG3 epsilon_CG3 sigma_CG3
   pair_coeff 4 4 epsilon_FG4 sigma_FG4 0.0         0.0
   pair_modify shift yes

In this example, the switching cutoffs for the first cluster (atom types 1
and 2) is defined explicitly in the pair_coeff command which overrides the
global values, while the second cluster (atom types 3 and 4) uses the global
definition from the pair_style command. The emphasis here is that the atom
types that belong to a specific cluster should have the same switching/cutoff
arguments.

In the case that geometric mixing is not applicable, for simulating the
system from the previous example, we recommend using the following commands:

.. code-block:: LAMMPS

   pair_style lj/relres Rsi Rso Rci Rco
   pair_coeff 1 1 epsilon_FG1  sigma_FG1  epsilon_CG1  sigma_CG1  Rsi1  Rso1  Rci Rco
   pair_coeff 1 2 epsilon_FG12 sigma_FG12 0.0          0.0        Rsi1  Rso1  Rci Rco
   pair_coeff 1 3 epsilon_FG13 sigma_FG13 epsilon_CG13 sigma_CG13 Rsi13 Rso13 Rci Rco
   pair_coeff 1 4 epsilon_FG14 sigma_FG14 0.0          0.0        Rsi13 Rso13 Rci Rco
   pair_coeff 2 2 epsilon_FG2  sigma_FG2  0.0          0.0        Rsi1  Rso1  Rci Rco
   pair_coeff 2 3 epsilon_FG23 sigma_FG23 0.0          0.0        Rsi13 Rso13 Rci Rco
   pair_coeff 2 4 epsilon_FG24 sigma_FG24 0.0          0.0        Rsi13 Rso13 Rci Rco
   pair_coeff 3 3 epsilon_FG3  sigma_FG3  epsilon_CG3  sigma_CG3
   pair_coeff 3 4 epsilon_FG34 sigma_FG34 0.0          0.0
   pair_coeff 4 4 epsilon_FG4  sigma_FG4  0.0          0.0
   pair_modify shift yes

Notice that the CG parameters are mixed only for interactions between atom
types associated with hybrid sites, and that the cutoffs are
mixed on the cluster basis.

More examples can be found in the *examples/relres* folder.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs :math:`I,\:J` with :math:`I\;!=J`, the
:math:`\epsilon^{FG}`, :math:`\sigma^{FG}`, :math:`\epsilon^{CG}`,
:math:`\sigma^{CG}`, :math:`r_{si}`, :math:`r_{so}`, :math:`r_{ci}`,
and :math:`r_{co}` parameters for this pair style can be mixed, if
not defined explicitly. All parameters are mixed according to the
pair_modify mix option.  The default mix value is *geometric*,
and it is recommended to use with this *lj/relres* style.  See the
"pair_modify" command for details.

This pair style supports the :doc:`pair_modify <pair_modify>` shift
option for the energy of the pair interaction. It is recommended to set
this option to *yes*\ .  Otherwise, the offset :math:`\Gamma_{c}`
is set to zero. Constants :math:`\Gamma_{si}` and :math:`\Gamma_{so}` are
not impacted by this option.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure, since the energy of the pair interaction is smoothed to 0.0
at the cutoff.

This pair style writes its information to :doc:`binary restart files
<restart>`, so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

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

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Chaimovich1:

**(Chaimovich1)** A. Chaimovich, C. Peter and K. Kremer, J. Chem. Phys. 143,
243107 (2015).

.. _Chaimovich2:

**(Chaimovich2)** M. Chaimovich and A. Chaimovich, J. Chem. Theory Comput. 17,
1045-1059 (2021).

