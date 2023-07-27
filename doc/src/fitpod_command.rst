.. index:: fitpod

fitpod command
======================

Syntax
""""""

.. code-block:: LAMMPS

   fitpod Ta_param.pod Ta_data.pod

* fitpod = style name of this command
* Ta_param.pod = an input file that describes proper orthogonal descriptors (PODs)
* Ta_data.pod = an input file that specifies DFT data used to fit a POD potential

Examples
""""""""

.. code-block:: LAMMPS

   fitpod Ta_param.pod Ta_data.pod

Description
"""""""""""
.. versionadded:: 22Dec2022

Fit a machine-learning interatomic potential (ML-IAP) based on proper
orthogonal descriptors (POD).  Two input files are required for this
command. The first input file describes a POD potential parameter
settings, while the second input file specifies the DFT data used for
the fitting procedure.

The table below has one-line descriptions of all the keywords that can
be used in the first input file (i.e. ``Ta_param.pod`` in the example
above):

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Keyword
     - Default
     - Type
     - Description
   * - species
     - (none)
     - STRING
     - Chemical symbols for all elements in the system and have to match XYZ training files.
   * - pbc
     - 1 1 1
     - INT
     - three integer constants specify boundary conditions
   * - rin
     - 1.0
     - REAL
     - a real number specifies the inner cut-off radius
   * - rcut
     - 5.0
     - REAL
     - a real number specifies the outer cut-off radius
   * - bessel_polynomial_degree
     - 3
     - INT
     - the maximum degree of Bessel polynomials
   * - inverse_polynomial_degree
     - 6
     - INT
     - the maximum degree of inverse radial basis functions
   * - onebody
     - 1
     - BOOL
     - turns on/off one-body potential
   * - twobody_number_radial_basis_functions
     - 6
     - INT
     - number of radial basis functions for two-body potential
   * - threebody_number_radial_basis_functions
     - 5
     - INT
     - number of radial basis functions for three-body potential
   * - threebody_number_angular_basis_functions
     - 5
     - INT
     - number of angular basis functions for three-body potential
   * - fourbody_snap_twojmax
     - 0
     - INT
     - band limit for SNAP bispectrum components (0,2,4,6,8... allowed)
   * - fourbody_snap_chemflag
     - 0
     - BOOL
     - turns on/off the explicit multi-element variant of the SNAP bispectrum components
   * - quadratic_pod_potential
     - 0
     - BOOL
     - turns on/off quadratic POD potential

All keywords except *species* have default values. If a keyword is not
set in the input file, its default value is used.  The next table
describes all keywords that can be used in the second input file
(i.e. ``Ta_data.pod`` in the example above):

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Keyword
     - Default
     - Type
     - Description
   * - file_format
     - extxyz
     - STRING
     - only the extended xyz format (extxyz) is currently supported
   * - file_extension
     - xyz
     - STRING
     - extension of the data files
   * - path_to_training_data_set
     - (none)
     - STRING
     - specifies the path to training data files in double quotes
   * - path_to_test_data_set
     - ""
     - STRING
     - specifies the path to test data files in double quotes
   * - fraction_training_data_set
     - 1.0
     - REAL
     - a real number (<= 1.0) specifies the fraction of the training set used to fit POD
   * - randomize_training_data_set
     - 0
     - BOOL
     - turns on/off randomization of the training set
   * - fitting_weight_energy
     - 100.0
     - REAL
     - a real constant specifies the weight for energy in the least-squares fit
   * - fitting_weight_force
     - 1.0
     - REAL
     - a real constant specifies the weight for force in the least-squares fit
   * - fitting_regularization_parameter
     - 1.0e-10
     - REAL
     - a real constant specifies the regularization parameter in the least-squares fit
   * - error_analysis_for_training_data_set
     - 0
     - BOOL
     - turns on/off error analysis for the training data set
   * - error_analysis_for_test_data_set
     - 0
     - BOOL
     - turns on/off error analysis for the test data set
   * - basename_for_output_files
     - pod
     - STRING
     - a basename string added to the output files
   * - precision_for_pod_coefficients
     - 8
     - INT
     - number of digits after the decimal points for numbers in the coefficient file

All keywords except *path_to_training_data_set* have default values. If
a keyword is not set in the input file, its default value is used.  After
successful training, a number of output files are produced, if enabled:

* ``<basename>_training_errors.pod``  reports the errors in energy and forces for the training data set
* ``<basename>_training_analysis.pod`` reports detailed errors for all training configurations
* ``<basename>_test_errors.pod`` reports errors for the test data set
* ``<basename>_test_analysis.pod`` reports detailed errors for all test configurations
* ``<basename>_coefficients.pod`` contains the coefficients of the POD potential

After training the POD potential, ``Ta_param.pod`` and ``<basename>_coefficients.pod``
are the two files needed to use the POD potential in LAMMPS. See
:doc:`pair_style pod <pair_pod>` for using the POD potential. Examples
about training and using POD potentials are found in the directory
lammps/examples/PACKAGES/pod.

Parameterized Potential Energy Surface
""""""""""""""""""""""""""""""""""""""

We consider a multi-element system of *N* atoms with :math:`N_{\rm e}`
unique elements.  We denote by :math:`\boldsymbol r_n` and :math:`Z_n`
position vector and type of an atom *n* in the system,
respectively. Note that we have :math:`Z_n \in \{1, \ldots, N_{\rm e}
\}`, :math:`\boldsymbol R = (\boldsymbol r_1, \boldsymbol r_2, \ldots,
\boldsymbol r_N) \in \mathbb{R}^{3N}`, and :math:`\boldsymbol Z = (Z_1,
Z_2, \ldots, Z_N) \in \mathbb{N}^{N}`. The potential energy surface
(PES) of the system can be expressed as a many-body expansion of the
form

.. math::

    E(\boldsymbol R, \boldsymbol Z, \boldsymbol{\eta}, \boldsymbol{\mu}) \ = \ & \sum_{i} V^{(1)}(\boldsymbol r_i, Z_i, \boldsymbol \mu^{(1)} ) + \frac12 \sum_{i,j} V^{(2)}(\boldsymbol r_i, \boldsymbol r_j, Z_i, Z_j, \boldsymbol \eta, \boldsymbol \mu^{(2)})  \\
    & + \frac16 \sum_{i,j,k} V^{(3)}(\boldsymbol r_i, \boldsymbol r_j, \boldsymbol r_k, Z_i, Z_j, Z_k, \boldsymbol \eta, \boldsymbol \mu^{(3)}) + \ldots

where :math:`V^{(1)}` is the one-body potential often used for
representing external field or energy of isolated elements, and the
higher-body potentials :math:`V^{(2)}, V^{(3)}, \ldots` are symmetric,
uniquely defined, and zero if two or more indices take identical values.
The superscript on each potential denotes its body order. Each *q*-body
potential :math:`V^{(q)}` depends on :math:`\boldsymbol \mu^{(q)}` which
are sets of parameters to fit the PES. Note that :math:`\boldsymbol \mu`
is a collection of all potential parameters :math:`\boldsymbol
\mu^{(1)}`, :math:`\boldsymbol \mu^{(2)}`, :math:`\boldsymbol
\mu^{(3)}`, etc, and that :math:`\boldsymbol \eta` is a set of
hyper-parameters such as inner cut-off radius :math:`r_{\rm in}` and
outer cut-off radius :math:`r_{\rm cut}`.

Interatomic potentials rely on parameters to learn relationship between
atomic environments and interactions.  Since interatomic potentials are
approximations by nature, their parameters need to be set to some
reference values or fitted against data by necessity.  Typically,
potential fitting finds optimal parameters, :math:`\boldsymbol \mu^*`,
to minimize a certain loss function of the predicted quantities and
data. Since the fitted potential depends on the data set used to fit it,
different data sets will yield different optimal parameters and thus
different fitted potentials. When fitting the same functional form on
*Q* different data sets, we would obtain *Q* different optimized
potentials, :math:`E(\boldsymbol R,\boldsymbol Z, \boldsymbol \eta,
\boldsymbol \mu_q^*), 1 \le q \le Q`.  Consequently, there exist many
different sets of optimized parameters for empirical interatomic
potentials.

Instead of optimizing the potential parameters, inspired by the reduced
basis method :ref:`(Grepl) <Grepl20072>` for parameterized partial
differential equations, we view the parameterized PES as a parametric
manifold of potential energies

.. math::

    \mathcal{M} = \{E(\boldsymbol R, \boldsymbol Z, \boldsymbol \eta, \boldsymbol \mu) \ | \  \boldsymbol \mu \in \Omega^{\boldsymbol \mu} \}

where :math:`\Omega^{\boldsymbol \mu}` is a parameter domain in which
:math:`\boldsymbol \mu` resides.  The parametric manifold
:math:`\mathcal{M}` contains potential energy surfaces for all values of
:math:`\boldsymbol \mu \in \Omega^{\boldsymbol \mu}`.  Therefore, the
parametric manifold yields a much richer and more transferable atomic
representation than any particular individual PES :math:`E(\boldsymbol
R, \boldsymbol Z, \boldsymbol \eta, \boldsymbol \mu^*)`.

We propose specific forms of the parameterized potentials for one-body,
two-body, and three-body interactions. We apply the Karhunen-Loeve
expansion to snapshots of the parameterized potentials to obtain sets of
orthogonal basis functions. These basis functions are aggregated
according to the chemical elements of atoms, thus leading to
multi-element proper orthogonal descriptors.

Proper Orthogonal Descriptors
"""""""""""""""""""""""""""""

Proper orthogonal descriptors are finger prints characterizing the
radial and angular distribution of a system of atoms. The detailed
mathematical definition is given in the paper by Nguyen and Rohskopf
:ref:`(Nguyen) <Nguyen20222>`.

The descriptors for the one-body interaction are used to capture energy
of isolated elements and defined as follows

.. math::

    D_{ip}^{(1)} =  \left\{
        \begin{array}{ll}
        1, & \mbox{if } Z_i = p \\
        0, & \mbox{if } Z_i \neq p
        \end{array}
    \right.

for :math:`1 \le i \le N, 1 \le p \le N_{\rm e}`. The number of one-body
descriptors per atom is equal to the number of elements. The one-body
descriptors are independent of atom positions, but dependent on atom
types. The one-body descriptors are active only when the keyword
*onebody* is set to 1.

We adopt the usual assumption that the direct interaction between two
atoms vanishes smoothly when their distance is greater than the outer
cutoff distance :math:`r_{\rm cut}`. Furthermore, we assume that two
atoms can not get closer than the inner cutoff distance :math:`r_{\rm
in}` due to the Pauli repulsion principle. Let :math:`r \in (r_{\rm in},
r_{\rm cut})`, we introduce the following parameterized radial functions

.. math::

    \phi(r, r_{\rm in}, r_{\rm cut}, \alpha, \beta)  = \frac{\sin (\alpha \pi x) }{r - r_{\rm in}}, \qquad  \varphi(r, \gamma)  = \frac{1}{r^\gamma} ,

where the scaled distance function :math:`x` is defined below to enrich the two-body manifold

.. math::

    x(r, r_{\rm in}, r_{\rm cut}, \beta) = \frac{e^{-\beta(r - r_{\rm in})/(r_{\rm cut} - r_{\rm in})} - 1}{e^{-\beta} - 1} .

We introduce the following function as a convex combination of the two functions

.. math::

    \psi(r, r_{\rm in}, r_{\rm cut}, \alpha, \beta, \gamma, \kappa)  = \kappa \phi(r, r_{\rm in}, r_{\rm cut}, \alpha, \beta) + (1- \kappa)  \varphi(r, \gamma) .

We see that :math:`\psi` is a function of distance :math:`r`, cut-off
distances :math:`r_{\rm in}` and :math:`r_{\rm cut}`, and parameters
:math:`\alpha, \beta, \gamma, \kappa`. Together these parameters allow
the function :math:`\psi` to characterize a diverse spectrum of two-body
interactions within the cut-off interval :math:`(r_{\rm in}, r_{\rm
cut})`.

Next, we introduce the following parameterized potential

.. math::

    W^{(2)}(r_{ij}, \boldsymbol \eta, \boldsymbol \mu^{(2)})  = f_{\rm c}(r_{ij}, \boldsymbol \eta) \psi(r_{ij}, \boldsymbol \eta, \boldsymbol \mu^{(2)})

where :math:`\eta_1 = r_{\rm in}, \eta_2 = r_{\rm cut}, \mu_1^{(2)} =
\alpha, \mu_2^{(2)} = \beta, \mu_3^{(2)} = \gamma`, and
:math:`\mu_4^{(2)} = \kappa`. Here the cut-off function :math:`f_{\rm
c}(r_{ij}, \boldsymbol \eta)` proposed in [refs] is used to ensure the
smooth vanishing of the potential and its derivative for :math:`r_{ij}
\ge r_{\rm cut}`:

.. math::

    f_{\rm c}(r_{ij},  r_{\rm in}, r_{\rm cut})  =  \exp \left(1 -\frac{1}{\sqrt{\left(1 - \frac{(r-r_{\rm in})^3}{(r_{\rm cut} - r_{\rm in})^3} \right)^2 + 10^{-6}}} \right)

Based on the parameterized potential, we form a set of snapshots as
follows.  We assume that we are given :math:`N_{\rm s}` parameter tuples
:math:`\boldsymbol \mu^{(2)}_\ell, 1 \le \ell \le N_{\rm s}`. We
introduce the following set of snapshots on :math:`(r_{\rm in}, r_{\rm
cut})`:

.. math::

    \xi_\ell(r_{ij}, \boldsymbol \eta) =  W^{(2)}(r_{ij}, \boldsymbol \eta, \boldsymbol \mu^{(2)}_\ell),  \quad \ell = 1, \ldots, N_{\rm s} .

To ensure adequate sampling of the PES for different parameters, we
choose :math:`N_{\rm s}` parameter points :math:`\boldsymbol
\mu^{(2)}_\ell = (\alpha_\ell, \beta_\ell, \gamma_\ell, \kappa_\ell), 1
\le \ell \le N_{\rm s}` as follows. The parameters :math:`\alpha \in [1,
N_\alpha]` and :math:`\gamma \in [1, N_\gamma]` are integers, where
:math:`N_\alpha` and :math:`N_\gamma` are the highest degrees for
:math:`\alpha` and :math:`\gamma`, respectively. We next choose
:math:`N_\beta` different values of :math:`\beta` in the interval
:math:`[\beta_{\min}, \beta_{\max}]`, where :math:`\beta_{\min} = 0` and
:math:`\beta_{\max} = 4`. The parameter :math:`\kappa` can be set either
0 or 1.  Hence, the total number of parameter points is :math:`N_{\rm s}
= N_\alpha N_\beta + N_\gamma`.  Although :math:`N_\alpha, N_\beta,
N_\gamma` can be chosen conservatively large, we find that
:math:`N_\alpha = 6, N_\beta = 3, N_\gamma = 8` are adequate for most
problems.  Note that :math:`N_\alpha` and :math:`N_\gamma` correspond to
*bessel_polynomial_degree* and *inverse_polynomial_degree*,
respectively.

We employ the Karhunen-Loeve (KL) expansion to generate an orthogonal
basis set which is known to be optimal for representation of the
snapshot family :math:`\{\xi_\ell\}_{\ell=1}^{N_{\rm s}}`. The two-body
orthogonal basis functions are computed as follows

.. math::

    U^{(2)}_m(r_{ij}, \boldsymbol \eta) = \sum_{\ell = 1}^{N_{\rm s}} A_{\ell m}(\boldsymbol \eta) \,  \xi_\ell(r_{ij}, \boldsymbol \eta), \qquad m = 1, \ldots, N_{\rm 2b} ,

where the matrix :math:`\boldsymbol A \in \mathbb{R}^{N_{\rm s} \times
N_{\rm s}}` consists of eigenvectors of the eigenvalue problem

.. math::

    \boldsymbol C \boldsymbol a = \lambda \boldsymbol a

with the entries of :math:`\boldsymbol C \in \mathbb{R}^{N_{\rm s} \times N_{\rm s}}` being given by

.. math::

    C_{ij}  = \frac{1}{N_{\rm s}} \int_{r_{\rm in}}^{r_{\rm cut}} \xi_i(x, \boldsymbol \eta) \xi_j(x, \boldsymbol \eta) dx, \quad 1 \le i, j \le N_{\rm s}

Note that the eigenvalues :math:`\lambda_\ell, 1 \le \ell \le N_{\rm
s}`, are ordered such that :math:`\lambda_1 \ge \lambda_2 \ge \ldots \ge
\lambda_{N_{\rm s}}`, and that the matrix :math:`\boldsymbol A` is
pe-computed and stored for any given :math:`\boldsymbol \eta`.  Owing to
the rapid convergence of the KL expansion, only a small number of
orthogonal basis functions is needed to obtain accurate
approximation. The value of :math:`N_{\rm 2b}` corresponds to
*twobody_number_radial_basis_functions*.

The two-body proper orthogonal descriptors at each atom *i* are computed
by summing the orthogonal basis functions over the neighbors of atom *i*
and numerating on the atom types as follows

.. math::

    D^{(2)}_{im l(p, q) }(\boldsymbol \eta)  = \left\{
    \begin{array}{ll}
    \displaystyle \sum_{\{j | Z_j = q\}} U^{(2)}_m(r_{ij},  \boldsymbol \eta), & \mbox{if } Z_i = p \\
    0, & \mbox{if } Z_i \neq p
    \end{array}
    \right.

for :math:`1 \le i \le N, 1 \le m \le N_{\rm 2b}, 1 \le q, p \le N_{\rm
e}`. Here :math:`l(p,q)` is a symmetric index mapping such that

.. math::

    l(p,q)  = \left\{
    \begin{array}{ll}
    q + (p-1) N_{\rm e} - p(p-1)/2, & \mbox{if } q \ge p \\
    p + (q-1) N_{\rm e} - q(q-1)/2, & \mbox{if } q < p .
    \end{array}
    \right.

The number of two-body descriptors per atom is thus :math:`N_{\rm 2b}
N_{\rm e}(N_{\rm e}+1)/2`.

It is important to note that the orthogonal basis functions do not
depend on the atomic numbers :math:`Z_i` and :math:`Z_j`. Therefore, the
cost of evaluating the basis functions and their derivatives with
respect to :math:`r_{ij}` is independent of the number of elements
:math:`N_{\rm e}`. Consequently, even though the two-body proper
orthogonal descriptors depend on :math:`\boldsymbol Z`, their
computational complexity is independent of :math:`N_{\rm e}`.

In order to provide proper orthogonal descriptors for three-body
interactions, we need to introduce a three-body parameterized
potential. In particular, the three-body potential is defined as a
product of radial and angular functions as follows

.. math::

    W^{(3)}(r_{ij}, r_{ik}, \theta_{ijk}, \boldsymbol \eta, \boldsymbol \mu^{(3)})  =  \psi(r_{ij}, r_{\rm min}, r_{\rm max}, \alpha, \beta, \gamma, \kappa) f_{\rm c}(r_{ij}, r_{\rm min}, r_{\rm max}) \\
    \psi(r_{ik}, r_{\rm min}, r_{\rm max}, \alpha, \beta, \gamma, \kappa) f_{\rm c}(r_{ik}, r_{\rm min}, r_{\rm max}) \\
    \cos (\sigma \theta_{ijk} + \zeta)

where :math:`\sigma` is the periodic multiplicity, :math:`\zeta` is the
equilibrium angle, :math:`\boldsymbol \mu^{(3)} = (\alpha, \beta,
\gamma, \kappa, \sigma, \zeta)`. The three-body potential provides an
angular fingerprint of the atomic environment through the bond angles
:math:`\theta_{ijk}` formed with each pair of neighbors :math:`j` and
:math:`k`.  Compared to the two-body potential, the three-body potential
has two extra parameters :math:`(\sigma, \zeta)` associated with the
angular component.

Let :math:`\boldsymbol \varrho = (\alpha, \beta, \gamma, \kappa)`. We
assume that we are given :math:`L_{\rm r}` parameter tuples
:math:`\boldsymbol \varrho_\ell, 1 \le \ell \le L_{\rm r}`.  We
introduce the following set of snapshots on :math:`(r_{\min},
r_{\max})`:

.. math::

    \zeta_\ell(r_{ij}, r_{\rm min}, r_{\rm max} ) =  \psi(r_{ij}, r_{\rm min}, r_{\rm max}, \boldsymbol \varrho_\ell) f_{\rm c}(r_{ij}, r_{\rm min},  r_{\rm max}), \quad 1 \le \ell \le L_{\rm r} .

We apply the Karhunen-Loeve (KL) expansion to this set of snapshots to
obtain orthogonal basis functions as follows

.. math::

    U^{r}_m(r_{ij}, r_{\rm min}, r_{\rm max} ) = \sum_{\ell = 1}^{L_{\rm r}} A_{\ell m} \,  \zeta_\ell(r_{ij}, r_{\rm min}, r_{\rm max} ), \qquad m = 1, \ldots, N_{\rm r} ,

where the matrix :math:`\boldsymbol A \in \mathbb{R}^{L_{\rm r} \times L_{\rm r}}` consists
of eigenvectors of the eigenvalue problem. For the parameterized angular function,
we consider angular basis functions

.. math::

    U^{a}_n(\theta_{ijk}) = \cos ((n-1) \theta_{ijk}), \qquad  n = 1,\ldots, N_{\rm a},

where :math:`N_{\rm a}` is the number of angular basis functions. The orthogonal
basis functions for the parameterized potential are computed as follows

.. math::

    U^{(3)}_{mn}(r_{ij}, r_{ik}, \theta_{ijk}, \boldsymbol \eta) = U^{r}_m(r_{ij}, \boldsymbol \eta) U^{r}_m(r_{ik}, \boldsymbol \eta) U^{a}_n(\theta_{ijk}),

for :math:`1 \le m \le N_{\rm r}, 1 \le n \le N_{\rm a}`. The number of three-body
orthogonal basis functions is equal to :math:`N_{\rm 3b} = N_{\rm r} N_{\rm a}` and
independent of the number of elements. The value of :math:`N_{\rm r}` corresponds to
*threebody_number_radial_basis_functions*, while that of :math:`N_{\rm a}` to
*threebody_number_angular_basis_functions*.

The three-body proper orthogonal descriptors at each atom *i*
are obtained by summing over the neighbors *j* and *k* of atom *i* as

.. math::

    D^{(3)}_{imn \ell(p, q, s)}(\boldsymbol \eta)  = \left\{
    \begin{array}{ll}
    \displaystyle \sum_{\{j | Z_j = q\}} \sum_{\{k | Z_k = s\}} U^{(3)}_{mn}(r_{ij}, r_{ik}, \theta_{ijk}, \boldsymbol \eta), & \mbox{if } Z_i = p \\
    0, & \mbox{if } Z_i \neq p
    \end{array}
    \right.

for :math:`1 \le i \le N, 1 \le m \le N_{\rm r}, 1 \le n \le N_{\rm a}, 1 \le q, p, s \le N_{\rm e}`,
where

.. math::

    \ell(p,q,s)  = \left\{
    \begin{array}{ll}
    s + (q-1) N_{\rm e} - q(q-1)/2 + (p-1)N_{\rm e}(1+N_{\rm e})/2 , & \mbox{if } s \ge q \\
    q + (s-1) N_{\rm e} - s(s-1)/2 + (p-1)N_{\rm e}(1+N_{\rm e})/2, & \mbox{if } s < q .
    \end{array}
    \right.

The number of three-body descriptors per atom is thus :math:`N_{\rm 3b} N_{\rm e}^2(N_{\rm e}+1)/2`.
While the number of three-body PODs is cubic function of the number of elements,
the computational complexity of the three-body PODs is independent of the number of elements.

Four-Body SNAP Descriptors
""""""""""""""""""""""""""

In addition to the proper orthogonal descriptors described above, we also employ
the spectral neighbor analysis potential (SNAP) descriptors. SNAP uses bispectrum components
to characterize the local neighborhood of each atom in a very general way. The mathematical definition
of the bispectrum calculation and its derivatives w.r.t. atom positions is described in
:doc:`compute snap <compute_sna_atom>`. In SNAP, the
total energy is decomposed into a sum over atom energies. The energy of
atom *i* is expressed as a weighted sum over bispectrum components.

.. math::

   E_i^{\rm SNAP} = \sum_{k=1}^{N_{\rm 4b}} \sum_{p=1}^{N_{\rm e}} c_{kp}^{(4)} D_{ikp}^{(4)}


where the SNAP descriptors are related to the bispectrum components by

.. math::

    D^{(4)}_{ikp}  = \left\{
    \begin{array}{ll}
    \displaystyle B_{ik}, & \mbox{if } Z_i = p \\
    0, & \mbox{if } Z_i \neq p
    \end{array}
    \right.

Here :math:`B_{ik}` is the *k*\ -th bispectrum component of atom *i*. The number of
bispectrum components :math:`N_{\rm 4b}` depends on the value of *fourbody_snap_twojmax* :math:`= 2 J_{\rm max}`
and *fourbody_snap_chemflag*. If *fourbody_snap_chemflag* = 0
then :math:`N_{\rm 4b} = (J_{\rm max}+1)(J_{\rm max}+2)(J_{\rm max}+1.5)/3`.
If *fourbody_snap_chemflag* = 1 then :math:`N_{\rm 4b} = N_{\rm e}^3 (J_{\rm max}+1)(J_{\rm max}+2)(J_{\rm max}+1.5)/3`.
The bispectrum calculation is described in more detail in :doc:`compute sna/atom <compute_sna_atom>`.

Linear Proper Orthogonal Descriptor Potentials
""""""""""""""""""""""""""""""""""""""""""""""

The proper orthogonal descriptors and SNAP descriptors are used to define the atomic energies
in the following expansion

.. math::

    E_{i}(\boldsymbol \eta) = \sum_{p=1}^{N_{\rm e}} c^{(1)}_p D^{(1)}_{ip} + \sum_{m=1}^{N_{\rm 2b}}  \sum_{l=1}^{N_{\rm e}(N_{\rm e}+1)/2} c^{(2)}_{ml} D^{(2)}_{iml}(\boldsymbol \eta) + \sum_{m=1}^{N_{\rm r}} \sum_{n=1}^{N_{\rm a}}  \sum_{\ell=1}^{N_{\rm e}^2(N_{\rm e}+1)/2} c^{(3)}_{mn\ell} D^{(3)}_{imn\ell}(\boldsymbol \eta) + \sum_{k=1}^{N_{\rm 4b}} \sum_{p=1}^{N_{\rm e}} c_{kp}^{(4)} D_{ikp}^{(4)}(\boldsymbol \eta),

where :math:`D^{(1)}_{ip}, D^{(2)}_{iml}, D^{(3)}_{imn\ell}, D^{(4)}_{ikp}` are the  one-body, two-body, three-body, four-body descriptors,
respectively, and :math:`c^{(1)}_p, c^{(2)}_{ml}, c^{(3)}_{mn\ell}, c^{(4)}_{kp}` are their respective expansion
coefficients. In a more compact notation that implies summation over descriptor indices
the atomic energies can be written as

.. math::

    E_i(\boldsymbol \eta) =  \sum_{m=1}^{N_{\rm e}} c^{(1)}_m D^{(1)}_{im} +  \sum_{m=1}^{N_{\rm d}^{(2)}} c^{(2)}_k D^{(2)}_{im} + \sum_{m=1}^{N_{\rm d}^{(3)}} c^{(3)}_m D^{(3)}_{im} + \sum_{m=1}^{N_{\rm d}^{(4)}} c^{(4)}_m D^{(4)}_{im}

where :math:`N_{\rm d}^{(2)} = N_{\rm 2b} N_{\rm e} (N_{\rm e}+1)/2`,
:math:`N_{\rm d}^{(3)} = N_{\rm 3b} N_{\rm e}^2 (N_{\rm e}+1)/2`, and
:math:`N_{\rm d}^{(4)} = N_{\rm 4b} N_{\rm e}` are
the number of two-body, three-body, and four-body descriptors, respectively.

The potential energy is then obtained by summing local atomic energies :math:`E_i`
for all atoms :math:`i` in the system

.. math::

    E(\boldsymbol \eta) = \sum_{i}^N E_{i}(\boldsymbol \eta)

Because the descriptors are one-body, two-body, and three-body terms,
the resulting POD potential is a three-body PES. We can express the potential
energy as a linear combination of the global descriptors as follows

.. math::

    E(\boldsymbol \eta) = \sum_{m=1}^{N_{\rm e}} c^{(1)}_m d^{(1)}_{m} +  \sum_{m=1}^{N_{\rm d}^{(2)}} c^{(2)}_m d^{(2)}_{m} + \sum_{m=1}^{N_{\rm d}^{(3)}} c^{(3)}_m d^{(3)}_{m} + \sum_{m=1}^{N_{\rm d}^{(4)}} c^{(4)}_m d^{(4)}_{m}

where  the global descriptors are given by

.. math::

    d_{m}^{(1)}(\boldsymbol \eta) = \sum_{i=1}^N D_{im}^{(1)}(\boldsymbol \eta), \quad d_{m}^{(2)}(\boldsymbol \eta) = \sum_{i=1}^N D_{im}^{(2)}(\boldsymbol \eta), \quad d_{m}^{(3)}(\boldsymbol \eta) = \sum_{i=1}^N D_{im}^{(3)}(\boldsymbol \eta), \quad d_{m}^{(4)}(\boldsymbol \eta) = \sum_{i=1}^N D_{im}^{(4)}(\boldsymbol \eta)

Hence, we obtain the atomic forces as

.. math::

    \boldsymbol F = -\nabla E(\boldsymbol \eta) = - \sum_{m=1}^{N_{\rm d}^{(2)}}  c^{(2)}_m  \nabla d_m^{(2)} - \sum_{m=1}^{N_{\rm d}^{(3)}}  c^{(3)}_m \nabla d_m^{(3)} - \sum_{m=1}^{N_{\rm d}^{(4)}}  c^{(4)}_m \nabla d_m^{(4)}

where :math:`\nabla d_m^{(2)}`, :math:`\nabla d_m^{(3)}` and :math:`\nabla d_m^{(4)}` are derivatives of the two-body
three-body, and four-body global descriptors with respect to atom positions, respectively.
Note that since the first-body global descriptors are constant, their derivatives are zero.

Quadratic Proper Orthogonal Descriptor Potentials
"""""""""""""""""""""""""""""""""""""""""""""""""

We recall two-body PODs :math:`D^{(2)}_{ik}, 1 \le k \le N_{\rm d}^{(2)}`,
and three-body PODs :math:`D^{(3)}_{im}, 1 \le m \le N_{\rm d}^{(3)}`,
with :math:`N_{\rm d}^{(2)} = N_{\rm 2b} N_{\rm e} (N_{\rm e}+1)/2` and
:math:`N_{\rm d}^{(3)} = N_{\rm 3b} N_{\rm e}^2 (N_{\rm e}+1)/2` being
the number of descriptors per atom for the two-body PODs and three-body PODs,
respectively. We employ them to define a new set of atomic descriptors as follows

.. math::

    D^{(2*3)}_{ikm} = \frac{1}{2N}\left( D^{(2)}_{ik} \sum_{j=1}^N D^{(3)}_{jm} + D^{(3)}_{im} \sum_{j=1}^N D^{(2)}_{jk}  \right)

for :math:`1 \le i \le N, 1 \le k \le N_{\rm d}^{(2)}, 1 \le m \le N_{\rm d}^{(3)}`.
The new descriptors are four-body because they involve central atom :math:`i` together
with three neighbors :math:`j, k` and :math:`l`. The total number of new  descriptors per atom is equal to

.. math::

    N_{\rm d}^{(2*3)} = N_{\rm d}^{(2)} * N_{\rm d}^{(3)} = N_{\rm 2b} N_{\rm 3b} N_{\rm e}^3 (N_{\rm e}+1)^2/4 .

The new global descriptors are calculated as

.. math::

    d^{(2*3)}_{km} = \sum_{i=1}^N D^{(2*3)}_{ikm} = \left( \sum_{i=1}^N D^{(2)}_{ik} \right) \left( \sum_{i=1}^N D^{(3)}_{im} \right) = d^{(2)}_{k} d^{(3)}_m,

for :math:`1 \le k \le N_{\rm d}^{(2)}, 1 \le m \le N_{\rm d}^{(3)}`. Hence, the gradient
of the new global descriptors with respect to atom positions is calculated as

.. math::

    \nabla d^{(2*3)}_{km} = d^{(3)}_m \nabla d^{(2)}_{k}  +  d^{(2)}_{k} \nabla d^{(3)}_m, \quad 1 \le k \le N_{\rm d}^{(2)}, 1 \le m \le N_{\rm d}^{(3)} .

The quadratic  POD potential is defined as a linear combination of the
original and new global descriptors as follows

.. math::

    E^{(2*3)} = \sum_{k=1}^{N_{\rm 2d}^{(2*3)}} \sum_{m=1}^{N_{\rm 3d}^{(2*3)}} c^{(2*3)}_{km} d^{(2*3)}_{km} .

It thus follows that

.. math::

    E^{(2*3)} = 0.5 \sum_{k=1}^{N_{\rm 2d}^{(2*3)}} \left( \sum_{m=1}^{N_{\rm 3d}^{(2*3)}} c^{(2*3)}_{km} d_m^{(3)} \right) d_k^{(2)} + 0.5 \sum_{m=1}^{N_{\rm 3d}^{(2*3)}} \left( \sum_{k=1}^{N_{\rm 2d}^{(2*3)}} c^{(2*3)}_{km} d_k^{(2)} \right) d_m^{(3)}  ,

which is simplified to

.. math::

    E^{(2*3)} =  0.5 \sum_{k=1}^{N_{\rm 2d}^{(2*3)}} b_k^{(2)} d_k^{(2)} +  0.5 \sum_{m=1}^{N_{\rm 3d}^{(2*3)}}   b_m^{(3)} d_m^{(3)}

where

.. math::

    b_k^{(2)} & = \sum_{m=1}^{N_{\rm 3d}^{(2*3)}} c^{(2*3)}_{km} d_m^{(3)}, \quad k = 1,\ldots, N_{\rm 2d}^{(2*3)}, \\
    b_m^{(3)} & = \sum_{k=1}^{N_{\rm 2d}^{(2*3)}} c^{(2*3)}_{km} d_k^{(2)}, \quad m = 1,\ldots, N_{\rm 3d}^{(2*3)} .

The quadratic POD potential results in the following atomic forces

.. math::

    \boldsymbol F^{(2*3)} = - \sum_{k=1}^{N_{\rm 2d}^{(2*3)}} \sum_{m=1}^{N_{\rm 3d}^{(2*3)}} c^{(2*3)}_{km}  \nabla d^{(2*3)}_{km} .

It can be shown that

.. math::

    \boldsymbol F^{(2*3)} = - \sum_{k=1}^{N_{\rm 2d}^{(2*3)}}   b^{(2)}_k \nabla d_k^{(2)} - \sum_{m=1}^{N_{\rm 3d}^{(2*3)}}  b^{(3)}_m  \nabla d_m^{(3)} .

The calculation of the atomic forces for the quadratic POD  potential
only requires the extra calculation of :math:`b_k^{(2)}` and :math:`b_m^{(3)}` which can be negligible.
As a result, the quadratic  POD potential does not increase the computational complexity.


Training
""""""""

POD potentials are trained using the least-squares regression against
density functional theory (DFT) data.  Let :math:`J` be the number of
training configurations, with :math:`N_j` being the number of atoms in
the j-th configuration. Let :math:`\{E^{\star}_j\}_{j=1}^{J}` and
:math:`\{\boldsymbol F^{\star}_j\}_{j=1}^{J}` be the DFT energies and
forces for :math:`J` configurations. Next, we calculate the global
descriptors and their derivatives for all training configurations. Let
:math:`d_{jm}, 1 \le m \le M`, be the global descriptors associated with
the j-th configuration, where :math:`M` is the number of global
descriptors. We then form a matrix :math:`\boldsymbol A \in
\mathbb{R}^{J \times M}` with entries :math:`A_{jm} = d_{jm}/ N_j` for
:math:`j=1,\ldots,J` and :math:`m=1,\ldots,M`.  Moreover, we form a
matrix :math:`\boldsymbol B \in \mathbb{R}^{\mathcal{N} \times M}` by
stacking the derivatives of the global descriptors for all training
configurations from top to bottom, where :math:`\mathcal{N} =
3\sum_{j=1}^{J} N_j`.

The coefficient vector :math:`\boldsymbol c` of the POD potential is
found by solving the following least-squares problem

.. math::

    {\min}_{\boldsymbol c \in \mathbb{R}^{M}} \ w_E \|\boldsymbol A(\boldsymbol \eta) \boldsymbol c - \bar{\boldsymbol E}^{\star} \|^2 + w_F \|\boldsymbol B(\boldsymbol \eta) \boldsymbol c + \boldsymbol F^{\star} \|^2 + w_R \|\boldsymbol c \|^2,

where :math:`w_E` and :math:`w_F` are weights for the energy
(*fitting_weight_energy*) and force (*fitting_weight_force*),
respectively; and :math:`w_R` is the regularization parameter (*fitting_regularization_parameter*).  Here :math:`\bar{\boldsymbol E}^{\star} \in
\mathbb{R}^{J}` is a vector of with entries :math:`\bar{E}^{\star}_j =
E^{\star}_j/N_j` and :math:`\boldsymbol F^{\star}` is a vector of
:math:`\mathcal{N}` entries obtained by stacking :math:`\{\boldsymbol
F^{\star}_j\}_{j=1}^{J}` from top to bottom.

The training procedure is the same for both the linear and quadratic POD
potentials.  However, since the quadratic POD potential has a
significantly larger number of the global descriptors, it is more
expensive to train the linear POD potential. This is because the
training of the quadratic POD potential still requires us to calculate
and store the quadratic global descriptors and their
gradient. Furthermore, the quadratic POD potential may require more
training data in order to prevent over-fitting. In order to reduce the
computational cost of fitting the quadratic POD potential and avoid
over-fitting, we can use subsets of two-body and three-body PODs for
constructing the new descriptors.


Restrictions
""""""""""""

This command is part of the ML-POD package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_style pod <pair_pod>`

Default
"""""""

The keyword defaults are also given in the description of the input files.

----------

.. _Grepl20072:

**(Grepl)** Grepl, Maday, Nguyen, and Patera, ESAIM: Mathematical Modelling and Numerical Analysis 41(3), 575-605, (2007).

.. _Nguyen20222:

**(Nguyen)** Nguyen and Rohskopf, arXiv preprint arXiv:2209.02362 (2022).
