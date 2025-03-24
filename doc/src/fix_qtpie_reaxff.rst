.. index:: fix qtpie/reaxff

fix qtpie/reaxff command
========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID qtpie/reaxff Nevery cutlo cuthi tolerance params gfile args

* ID, group-ID are documented in :doc:`fix <fix>` command
* qtpie/reaxff = style name of this fix command
* Nevery = perform QTPIE every this many steps
* cutlo,cuthi = lo and hi cutoff for Taper radius
* tolerance = precision to which charges will be equilibrated
* params = reaxff or a filename
* gfile = the name of a file containing Gaussian orbital exponents
* one or more keywords or keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *maxiter*
       *maxiter* N = limit the number of iterations to *N*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all qtpie/reaxff 1 0.0 10.0 1.0e-6 reaxff exp.qtpie
   fix 1 all qtpie/reaxff 1 0.0 10.0 1.0e-6 params.qtpie exp.qtpie maxiter 500

Description
"""""""""""

.. versionadded:: 19Nov2024

The QTPIE charge equilibration method is an extension of the QEq charge
equilibration method. With QTPIE, the partial charges on individual atoms
are computed by minimizing the electrostatic energy of the system in the
same way as the QEq method but where the absolute electronegativity,
:math:`\chi_i`, of each atom in the QEq charge equilibration scheme
:ref:`(Rappe and Goddard) <Rappe3>` is replaced with an effective
electronegativity given by :ref:`(Chen) <qtpie-Chen>`

.. math::
   \chi_{\mathrm{eff},i} = \frac{\sum_{j=1}^{N} (\chi_i - \chi_j) S_{ij}}
                                {\sum_{m=1}^{N}S_{im}},

which acts to penalize long-range charge transfer seen with the QEq charge
equilibration scheme. In this equation, :math:`N` is the number of atoms in
the system and :math:`S_{ij}` is the overlap integral between atom :math:`i`
and atom :math:`j`.

The effect of an external electric field can be incorporated into the QTPIE
method by modifying the absolute or effective electronegativities of each
atom :ref:`(Chen) <qtpie-Chen>`. This fix models the effect of an external
electric field by using the effective electronegativity given in
:ref:`(Gergs) <Gergs>`:

.. math::
   \chi_{\mathrm{eff},i} = \frac{\sum_{j=1}^{N} (\chi_i - \chi_j + \phi_i - \phi_j) S_{ij}}
                                {\sum_{m=1}^{N}S_{im}},

where :math:`\phi_i` and :math:`\phi_j` are the electric
potentials at the positions of atom :math:`i` and :math:`j`
due to the external electric field.

This fix is typically used in conjunction with the ReaxFF force
field model as implemented in the :doc:`pair_style reaxff <pair_reaxff>`
command, but it can be used with any potential in LAMMPS, so long as it
defines and uses charges on each atom. For more technical details about the
charge equilibration performed by `fix qtpie/reaxff`, which is the same as in
:doc:`fix qeq/reaxff <fix_qeq_reaxff>` except for the use of
:math:`\chi_{\mathrm{eff},i}`, please refer to :ref:`(Aktulga) <qeq-Aktulga2>`.
To be explicit, this fix replaces :math:`\chi_k` of eq. 3 in
:ref:`(Aktulga) <qeq-Aktulga2>` with :math:`\chi_{\mathrm{eff},k}`.

This fix requires the absolute electronegativity, :math:`\chi`, in eV, the
self-Coulomb potential, :math:`\eta`, in eV, and the shielded Coulomb
constant, :math:`\gamma`, in :math:`\AA^{-1}`. If the *params* setting above
is the word "reaxff", then these are extracted from the
:doc:`pair_style reaxff <pair_reaxff>` command and the ReaxFF force field
file it reads in.  If a file name is specified for *params*, then the
parameters are taken from the specified file and the file must contain
one line for each atom type.  The latter form must be used when performing
QTPIE with a non-ReaxFF potential. Each line should be formatted as follows,
ensuring that the parameters are given in units of eV, eV, and :math:`\AA^{-1}`,
respectively:

.. parsed-literal::

   itype chi eta gamma

where *itype* is the atom type from 1 to Ntypes. Note that eta is
defined here as twice the eta value in the ReaxFF file.

The overlap integrals in the equation for :math:`\chi_{\mathrm{eff},i}`
are computed by using normalized 1s Gaussian type orbitals. The Gaussian
orbital exponents, :math:`\alpha`, that are needed to compute the overlap
integrals are taken from the file given by *gfile*.
This file must contain one line for each atom type and provide the Gaussian
orbital exponent for each atom type in units of inverse square Bohr radius.
Each line should be formatted as follows:

.. parsed-literal::

   itype alpha

Empty lines or any text following the pound sign (#) are ignored. An example
*gfile* for a system with two atom types is

.. parsed-literal::

    # An example gfile. Exponents are taken from Table 2.2 of Chen, J. (2009).
    # Theory and applications of fluctuating-charge models.
    # The units of the exponents are 1 / (Bohr radius)^2 .
    1  0.2240  # O
    2  0.5434  # H

The optional *maxiter* keyword allows changing the max number
of iterations in the linear solver. The default value is 200.

.. note::

   In order to solve the self-consistent equations for electronegativity
   equalization, LAMMPS imposes the additional constraint that all the
   charges in the fix group must add up to zero.  The initial charge
   assignments should also satisfy this constraint.  LAMMPS will print a
   warning if that is not the case.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  This fix computes a global scalar (the number of
iterations) and a per-atom vector (the effective electronegativity), which
can be accessed by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

This fix is invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the REAXFF package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

This fix does not correctly handle interactions involving multiple
periodic images of the same atom.  Hence, it should not be used for
periodic cell dimensions smaller than the non-bonded cutoff radius,
which is typically :math:`10~\AA` for ReaxFF simulations.

This fix may be used in combination with :doc:`fix efield <fix_efield>`
and will apply the external electric field during charge equilibration,
but there may be only one fix efield instance used and the electric field
must be applied to all atoms in the system. Consequently, `fix efield` must
be used with *group-ID* all and must not be used with the keyword *region*.
Equal-style variables can be used for electric field vector
components without any further settings. Atom-style variables can be used
for spatially-varying electric field vector components, but the resulting
electric potential must be specified as an atom-style variable using
the *potential* keyword for `fix efield`.

Related commands
""""""""""""""""

:doc:`pair_style reaxff <pair_reaxff>`, :doc:`fix qeq/reaxff <fix_qeq_reaxff>`,
:doc:`fix acks2/reaxff <fix_acks2_reaxff>`

Default
"""""""

maxiter 200

----------

.. _Rappe3:

**(Rappe)** Rappe and Goddard III, Journal of Physical Chemistry, 95,
3358-3363 (1991).

.. _qtpie-Chen:

**(Chen)** Chen, Jiahao. Theory and applications of fluctuating-charge models.
University of Illinois at Urbana-Champaign, 2009.

.. _Gergs:

**(Gergs)** Gergs, Dirkmann and Mussenbrock.
Journal of Applied Physics 123.24 (2018).

.. _qeq-Aktulga2:

**(Aktulga)** Aktulga, Fogarty, Pandit, Grama, Parallel Computing, 38,
245-259 (2012).
