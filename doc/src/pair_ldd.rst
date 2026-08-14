.. index:: pair_style ldd

pair_style ldd command
======================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style ldd
   pair_coeff * * file.ldd S1 S2 ...

* ldd = name of this pair style (takes no arguments)
* file.ldd = name of the local-density potential file
* S1, S2, ... = mapping of N LAMMPS atom types to species names from the file (N = # of atom types)

Examples
""""""""

.. code-block:: LAMMPS

   # two atom types mapped to species A and B
   pair_style ldd
   pair_coeff * * LDINDSET1.ldd A B

   # combine local-density interactions with other pair styles
   pair_style hybrid/overlay ldd lj/cut 14.0
   pair_coeff 1 1 lj/cut 1.0 1.0
   pair_coeff * * ldd ld_input_file.ldd A

A potential file entry has the form (one ordered species pair per line):

.. parsed-literal::

   *Si Sj* indicator *wtype r0 rc* self *yes/no* potential *ptype args* [gradient *gtype args*] [ignore]

Description
"""""""""""

.. versionadded:: 4Jul2026

Style *ldd* implements the local density potential as first described by
Pagonabarraga and Frenkel :ref:`(Pagonabarraga)<Pagonabarraga>` and
additionally the square gradient of local densities first introduced by
:ref:`(DeLyser)<DeLyser>`.  The pair_style *ldd* is compatible with a
variety of molecular and atomic topologies (see the :doc:`Howto ldd
<Howto_ldd>` page for details) and offers a variety of options for how to
define the local density.

Following the manybody potential convention (as for :doc:`pair_style sw
<pair_sw>` or :doc:`tersoff <pair_tersoff>`), all interactions are read
from a potential file with a single ``pair_coeff * *`` command.  The file
defines a set of *species*; the arguments after the file name map each
LAMMPS atom type to one of these species (one species name per atom type,
in order).  Several atom types may map to the same species.  The set and
number of species are inferred from the file.  The cutoff of each
interaction is the ``rc`` of its indicator function; there is no separate
global cutoff.

Each line of the potential file specifies the interaction for one
*ordered* species pair ``Si Sj`` (the local density of species ``Sj``
around a central atom of species ``Si``), so the ``Si Sj`` and ``Sj Si``
interactions can differ.  All :math:`N_{\text{species}}^2` ordered pairs
must be listed (use the *ignore* keyword for pairs with no interaction).
Comment lines (starting with ``#``) and the customary
``# DATE: ... UNITS: ... CONTRIBUTOR: ...`` header are allowed.

Theory
""""""

Here, for notational simplicity, we outline the theory for style *ldd*
potentials with just one species of particle.  The :doc:`Howto ldd
<Howto_ldd>` page explains the more general
:math:`n_{\text{species}}` case.

Consider an indicator function of pair distance, :math:`w(r)`, where
:math:`w(0)=1`, :math:`w(r_{c})=0`, and :math:`w(r)` is
continuous/differentiable for :math:`0 \leq r \leq r_{c}`.  We define
:math:`[w]` as the spatial integral of :math:`w(r)` and the normalized
indicator function as :math:`\bar{w}(r) = w(r)/[w]`.  The local density
around a given particle is then defined as :math:`\rho_{I} = \sum_{J}
\bar{w}(r_{IJ})`, where the sum over :math:`J` can either include or
exclude :math:`I`, depending on whether the argument following the
*self* keyword is yes or no.  Practically, the only difference is a
horizontal shift in the argument of the potential by :math:`1/[w]`. The
local density potential and corresponding pair forces are given by:

.. math::
   U_{LD}(\mathbf{R}) &= \sum_{I} U_{\rho}(\rho_{I}) \\
   F_{IJ;\rho}(\mathbf{R}) &= (F_{\rho}(\rho_{I}) + F_{\rho}(\rho_{J}))\bar{w}'(r_{IJ}) \mathbf{e}_{IJ}

where :math:`F_{\rho}(x) = -dU_{\rho}(x)/dx`.

The optional *gradient* keyword implements a potential that is a
function of the square gradient of the local density as described in
:ref:`(DeLyser)<DeLyser>`.  This keyword adds an extra square gradient
(SG) term to the overall potential:

.. math::
   U_{SG}(\mathbf{R}) = \sum_{I} U_{\nabla}(\rho_{I}) \| \nabla_{I} \rho_{I} \|^{2}

The *self* argument indicates whether the particle :math:`I=J` term is
included in the local densities and gradients calculated.  Note that if
site :math:`I` is of species :math:`\alpha`, then it cannot contribute to
the local density of :math:`\beta` sites around :math:`I` when
:math:`\beta \neq \alpha`.  Thus in that case the self term is
automatically excluded and a warning is generated.

The *ignore* keyword is used in simulations where only some of the
species pairs have local density potentials acting between them.  Once a
pair is marked *ignore*, no local density or gradient is computed for it.

The *noforce* potential is an alternative way to turn an interaction off:
it computes the local densities and gradients (so they can be inspected,
see below) without applying any force.

.. _ldd_indicator:

Indicator functions
"""""""""""""""""""""

The (required) *indicator* keyword defines the form for :math:`w(r)` and
takes three arguments, ``wtype r0 rc``: the indicator style, the inner
range :math:`r_0`, and the outer cutoff :math:`r_c`.  In every case
:math:`\theta(x)` is the Heaviside function (:math:`1` for
:math:`x \geq 0`, :math:`0` otherwise).

.. index:: pair_coeff ldd indicator lucy

*lucy* employs the Lucy function popularized from smoothed particle
hydrodynamics :ref:`(Lucy)<Lucy>` (requires :math:`r_0 = 0`):

.. math::
   w(r) &= (1-r/r_{c})^{3}(1+3r/r_{c}) \\
   [w] &= 16 \pi r_{c}^{3} / 105 \\
   \bar{w}(r) &= 105 (1-r/r_{c})^{3}(1+3r/r_{c}) \theta(r_c - r)/(16 \pi r_{c}^{3})

.. index:: pair_coeff ldd indicator dpd

*dpd* employs a simple quadratic indicator function (requires
:math:`r_0 = 0`):

.. math::
   w(r) &= (1-r/r_{c})^{2} \\
   [w] &= 2 \pi r_{c}^{3} / 15 \\
   \bar{w}(r) &= 15 (1-r/r_{c})^{2} \theta(r_{c}-r) /(2 \pi r_{c}^{3})

.. index:: pair_coeff ldd indicator sphere

*sphere* employs a smoothed Heaviside indicator function.  If you
consider a large sphere of radius :math:`X` and a small sphere of radius
:math:`x`, with :math:`r_{0} = X - x` and :math:`r_{c} = X + x`, then
:math:`w(r)` is the fraction of the small sphere contained in the large
sphere when their centers are separated by :math:`r`:

.. math::
   w(r) =
   \begin{cases}
   1 & r \leq r_{0} \\
   c_{3}r^{3} + c_{1} r + c_{0} + c_{-1}/r  &r_{0} \leq r \leq r_{c} \\
   0 & r \geq r_{c}
   \end{cases}

.. index:: pair_coeff ldd indicator shell

*shell* employs a smoothed Heaviside indicator function originally
parameterized in :ref:`(Sanyal)<Sanyal>`.  It is one for
:math:`r \leq r_0`, zero for :math:`r \geq r_c`, and smoothly
interpolates in between with a continuous first derivative (but a
discontinuous second derivative at :math:`r_0` and :math:`r_c`):

.. math::
   w(r) =
   \begin{cases}
   1 & r \leq r_{0} \\
   c_{0} + c_{2} r^{2} + c_{4} r^{4} + c_{6}r^{6} & r_{0} \leq r \leq r_{c} \\
   0 & r \geq r_{c}
   \end{cases}

.. index:: pair_coeff ldd indicator smooth

*smooth* is similar to *shell* but has a continuous second derivative:

.. math::
   w(r) =
   \begin{cases}
   1 & r \leq r_{0} \\
   \sum_{i=0}^{5} c_{i}r^{i} &r_{0} \leq r \leq r_{c} \\
   0 & r \geq r_{c}
   \end{cases}

.. _ldd_potential:

Potential and gradient functions
""""""""""""""""""""""""""""""""

The (required) *potential* keyword defines the form for
:math:`U_{\rho}`, and the (optional) *gradient* keyword defines the form
for :math:`U_{\nabla}`.  Both take a style name followed by
style-specific arguments.  Below, :math:`X` is a placeholder for either
:math:`\rho` or :math:`\nabla`, and :math:`F_{X} = -dU_{X}/d\rho`.

.. index:: pair_coeff ldd potential noforce

*noforce* (potential only) computes the densities/gradients without
applying a force: :math:`U_{\rho}(\rho) = 0`, :math:`F_{\rho}(\rho) = 0`.

.. index:: pair_coeff ldd potential constant

*constant* takes one argument :math:`c`.  A constant *potential* does not
change the dynamics, but a constant *gradient* does:

.. math::
   U_{X}(\rho) = c, \qquad F_{X}(\rho) = 0.

.. index:: pair_coeff ldd potential linear

*linear* takes two arguments, the slope :math:`m` and intercept
:math:`b`:

.. math::
   U_{X}(\rho) = m\rho + b, \qquad F_{X}(\rho) = -m.

.. index:: pair_coeff ldd potential quadratic

*quadratic* takes three arguments :math:`a`, :math:`b`, :math:`c`:

.. math::
   U_{X}(\rho) = a\rho^{2} + b\rho + c, \qquad F_{X}(\rho) = -2a\rho - b.

.. index:: pair_coeff ldd potential mdpd

*mdpd* takes one argument :math:`B` and applies the quadratic local
density potential used in many MDPD studies (e.g. :ref:`(Warren)<Warren>`
or :ref:`(Ghoufi)<Ghoufi>`).  Combined with the *dpd* indicator and
:doc:`pair_style dpd <pair_dpd>` it reproduces the traditional
many-body DPD conservative force:

.. math::
   U_{X}(\rho) = B \pi r_{c}^{4} \rho^{2} / 30, \qquad
   F_{X}(\rho) = - B \pi r_{c}^{4} \rho / 15.

.. index:: pair_coeff ldd potential table/lin
.. index:: pair_coeff ldd potential table/spline

*table/lin*, *table/spline*, *table/gradlin*, *table/gradspline* each
take one argument, the name of a tabulated potential file.  The
``table/lin`` and ``table/spline`` forms interpolate :math:`U_{X}` from
the table linearly or with a cubic spline; the ``table/gradlin`` and
``table/gradspline`` forms are the corresponding interpolating styles
used after the *gradient* keyword.  Each line of a table file holds three numbers,
:math:`\rho`, :math:`U(\rho)`, and :math:`-dU/d\rho`, on a uniform
:math:`\rho` grid.  Example tables are provided under
``examples/PACKAGES/bocs``.

Accessing the local-density data
""""""""""""""""""""""""""""""""

The per-atom local densities, density gradients, and energies are
recomputed every step but are not stored in data or restart files.  They
can be made available to the rest of LAMMPS with the :doc:`fix pair
<fix_pair>` command, which exposes the following per-atom fields (one
column per species, in the order of the type map; *grad_density* has
three columns per species):

* *local_density* -- local density of each species around an atom
* *grad_density* -- gradient of the local density (x, y, z per species)
* *energy* -- :math:`U_{\rho}` contribution from each species
* *grad_energy* -- :math:`U_{\nabla}` contribution from each species
* *total_energy* -- sum of all local-density energies on an atom (scalar)

For example, ``fix f all pair 1 ldd local_density 0 total_energy 0``
creates the per-atom array ``f_f`` whose first columns are the
per-species local densities followed by the total energy; these can then
be output with :doc:`dump custom <dump>` or read by other per-atom
computes and fixes.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support automatic mixing.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart
files <restart>`.  Therefore, you must re-specify the pair_style and
pair_coeff commands in an input script that reads a restart file.

Restrictions
""""""""""""

This pair style is part of the BOCS package. It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This pair style requires ``newton on`` for pair interactions.

The ``pair_style ldd`` command must be issued *after* the simulation box
has been created (e.g. with :doc:`read_data <read_data>`,
:doc:`read_restart <read_restart>`, or :doc:`create_box <create_box>`),
because it sizes its per-type species map and inter-processor communication
buffers from the number of atom types.  Since ``ldd`` is a manybody-style
potential that does not read per-type ``Pair Coeffs`` from a data file,
there is no need to define it earlier.  When continuing from a binary
restart, re-specify ``pair_style ldd`` and its ``pair_coeff`` settings after
the :doc:`read_restart <read_restart>` command.

The *indicator*, *self*, and *potential* keywords are mandatory for each
species pair unless the *ignore* keyword is provided; the *gradient*
keyword is optional.  Every ordered species pair must appear exactly once
in the potential file.

Indicator styles with a non-zero :math:`r_0` (``sphere``, ``shell``,
``smooth``) are non-zero inside :math:`r_0`; *lucy* and *dpd* require
:math:`r_0 = 0`.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`fix pair <fix_pair>`,
:doc:`Howto ldd <Howto_ldd>`

----------

.. _Pagonabarraga:

**(Pagonabarraga)** I. Pagonabarraga, D. Frenkel. "Dissipative particle dynamics for interacting systems.", J. Chem. Phys., 115, 5015 (2001).

.. _DeLyser:

**(DeLyser)** M.R. DeLyser, W.G. Noid. "Coarse-grained models for local density gradients." J. Chem. Phys., 156, 034106 (2021).

.. _Lucy:

**(Lucy)** L. B. Lucy, "A numerical approach to the testing of the fission hypothesis.", Astronomical Journal, 82, 1013-1024 (1977).

.. _Sanyal:

**(Sanyal)** T. Sanyal, M. S. Shell, "Coarse-grained models using local-density potentials optimized with the relative entropy: Application to implicit solvation.", J. Chem. Phys., 145, 034109 (2016).

.. _Warren:

**(Warren)** P. B. Warren, "Vapor-liquid coexistence in many-body dissipative particle dynamics.", Phys. Rev. E, 68, 066702 (2003).

.. _Ghoufi:

**(Ghoufi)** A. Ghoufi, P. Malfreyt, "Mesoscale modeling of the water liquid-vapor interface: A surface tension calculation.", Phys. Rev. E, 83, 051601 (2011).
