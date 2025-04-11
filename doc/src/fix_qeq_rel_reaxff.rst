.. index:: fix qeq/rel/reaxff

fix qeq/rel/reaxff command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID qeq/rel/reaxff Nevery cutlo cuthi tolerance params gfile args

* ID, group-ID are documented in :doc:`fix <fix>` command
* qeq/rel/reaxff = style name of this fix command
* Nevery = perform QEqR every this many steps
* cutlo,cuthi = lo and hi cutoff for Taper radius
* tolerance = precision to which charges will be equilibrated
* params = reaxff or a filename
* gfile = the name of a file containing Gaussian orbital exponents
* one or more keywords or keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *scale* or *maxiter* or *nowarn*
       *scale* beta = set value of scaling factor *beta* (determines strength of electric polarization)
       *maxiter* N = limit the number of iterations to *N*
       *nowarn* = do not print a warning message if the maximum number of iterations is reached

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all qeq/rel/reaxff 1 0.0 10.0 1.0e-6 reaxff exp.qeqr
   fix 1 all qeq/rel/reaxff 1 0.0 10.0 1.0e-6 params.qeqr exp.qeqr scale 1.5 maxiter 500 nowarn

Description
"""""""""""

.. versionadded:: 19Nov2024

This fix implements the QEqR method for charge equilibration, which
differs from the QEq charge equilibration method :ref:`(Rappe and
Goddard) <Rappe4>` only in how external electric fields are accounted
for.  This fix therefore raises a warning when used without :doc:`fix
efield <fix_efield>` since :doc:`fix qeq/reaxff <fix_qeq_reaxff>` should
be used without an external electric field.  Charges are computed with
the QEqR method by minimizing the electrostatic energy of the system in
the same way as the QEq method but where the absolute electronegativity,
:math:`\chi_i`, of each atom in the QEq method is replaced with an
effective electronegativity given by

.. math::
   \chi_{\mathrm{r}i} = \chi_i + \frac{\sum_{j=1}^{N} \beta(\phi_i - \phi_j) S_{ij}}
                                      {\sum_{m=1}^{N}S_{im}},

where :math:`N` is the number of atoms in the system, :math:`\beta` is a
scaling factor, :math:`\phi_i` and :math:`\phi_j` are the electric
potentials at the positions of atoms :math:`i` and :math:`j` due to the
external electric field and :math:`S_{ij}` is the overlap integral
between atoms :math:`i` and :math:`j`.  This formulation is advantageous
over the method used by :doc:`fix qeq/reaxff <fix_qeq_reaxff>` to
account for an external electric field in that it permits periodic
boundaries in the direction of an external electric field and in that it
does not worsen long-range charge transfer seen with QEq.

This fix is typically used in conjunction with the ReaxFF force field
model as implemented in the :doc:`pair_style reaxff <pair_reaxff>`
command, but it can be used with any potential in LAMMPS, so long as it
defines and uses charges on each atom.  For more technical details about
the charge equilibration performed by *fix qeq/rel/reaxff*, which is the
same as in :doc:`fix qeq/reaxff <fix_qeq_reaxff>` except for the use of
:math:`\chi_{\mathrm{r}i}`, please refer to :ref:`(Aktulga)
<qeq-Aktulga3>`.  To be explicit, *fix qeq/rel/reaxff* replaces
:math:`\chi_k` of eq. 3 in :ref:`(Aktulga) <qeq-Aktulga3>` with
:math:`\chi_{\mathrm{r}k}` when an external electric field is applied.

This fix requires the absolute electronegativity, :math:`\chi`, in eV,
the self-Coulomb potential, :math:`\eta`, in eV, and the shielded
Coulomb constant, :math:`\gamma`, in :math:`\AA^{-1}`.  If the *params*
setting above is the word "reaxff", then these are extracted from the
:doc:`pair_style reaxff <pair_reaxff>` command and the ReaxFF force
field file it reads in.  If a file name is specified for *params*, then
the parameters are taken from the specified file and the file must
contain one line for each atom type.  The latter form must be used when
using this fix with a non-ReaxFF potential. Each line should be
formatted as follows, ensuring that the parameters are given in units of
eV, eV, and :math:`\AA^{-1}`, respectively:

.. parsed-literal::

   itype chi eta gamma

where *itype* is the atom type from 1 to Ntypes. Note that eta is
defined here as twice the eta value in the ReaxFF file.

The overlap integrals :math:`S_{ij}` are computed by using normalized 1s
Gaussian type orbitals. The Gaussian orbital exponents, :math:`\alpha`,
that are needed to compute the overlap integrals are taken from the file
given by *gfile*.  This file must contain one line for each atom type
and provide the Gaussian orbital exponent for each atom type in units of
inverse square Bohr radius.  Each line should be formatted as follows:

.. parsed-literal::

   itype alpha

Empty lines or any text following the pound sign (#) are ignored. An
example *gfile* for a system with two atom types is

.. parsed-literal::

    # An example gfile. Exponents are taken from Table 2.2 of Chen, J. (2009).
    # Theory and applications of fluctuating-charge models.
    # The units of the exponents are 1 / (Bohr radius)^2 .
    1  0.2240  # O
    2  0.5434  # H

The optional *scale* keyword sets the value of :math:`\beta` in the
equation for :math:`\chi_{\mathrm{r}i}`. The default value is 1.0.

The optional *maxiter* keyword allows changing the max number of
iterations in the linear solver. The default value is 200.

The optional *nowarn* keyword silences the warning message printed when
the maximum number of iterations is reached.  This can be useful for
comparing serial and parallel results where having the same fixed number
of iterations is desired, which can be achieved by using a very small
tolerance and setting *maxiter* to the desired number of iterations.

.. note::

   In order to solve the self-consistent equations for electronegativity
   equalization, LAMMPS imposes the additional constraint that all the
   charges in the fix group must add up to zero. The initial charge
   assignments should also satisfy this constraint. LAMMPS will print a
   warning if that is not the case.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  This fix computes a global scalar (the number of
iterations) and a per-atom vector (the effective electronegativity),
which can be accessed by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

This fix is invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the REAXFF package.  It is only enabled if LAMMPS
was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

This fix does not correctly handle interactions involving multiple
periodic images of the same atom.  Hence, it should not be used for
periodic cell dimensions smaller than the non-bonded cutoff radius,
which is typically :math:`10~\AA` for ReaxFF simulations.

This fix may be used in combination with :doc:`fix efield <fix_efield>`
and will apply the external electric field during charge equilibration,
but there may be only one fix efield instance used and the electric
field must be applied to all atoms in the system. Consequently, `fix
efield` must be used with *group-ID* all and must not be used with the
keyword *region*.  Equal-style variables can be used for electric field
vector components without any further settings. Atom-style variables can
be used for spatially-varying electric field vector components, but the
resulting electric potential must be specified as an atom-style variable
using the *potential* keyword for `fix efield`.

Related commands
""""""""""""""""

:doc:`pair_style reaxff <pair_reaxff>`, :doc:`fix qeq/reaxff <fix_qeq_reaxff>`,
:doc:`fix acks2/reaxff <fix_acks2_reaxff>`, :doc:`fix qtpie/reaxff <fix_qtpie_reaxff>`

Default
"""""""

scale = 1.0 and maxiter = 200

----------

.. _Rappe4:

**(Rappe)** Rappe and Goddard III, Journal of Physical Chemistry, 95,
3358-3363 (1991).

.. _qeq-Aktulga3:

**(Aktulga)** Aktulga, Fogarty, Pandit, Grama, Parallel Computing, 38,
245-259 (2012).
