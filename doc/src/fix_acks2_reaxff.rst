.. index:: fix acks2/reaxff
.. index:: fix acks2/reaxff/kk

fix acks2/reaxff command
========================

Accelerator Variants: *acks2/reaxff/kk*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID acks2/reaxff Nevery cutlo cuthi tolerance params args

* ID, group-ID are documented in :doc:`fix <fix>` command
* acks2/reaxff = style name of this fix command
* Nevery = perform ACKS2 every this many steps
* cutlo,cuthi = lo and hi cutoff for Taper radius
* tolerance = precision to which charges will be equilibrated
* params = reaxff or a filename
* one or more keywords or keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *maxiter* or *nowarn*
       *maxiter* N = limit the number of iterations to *N*
       *nowarn* = do not print a warning message if the maximum number of iterations was reached

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all acks2/reaxff 1 0.0 10.0 1.0e-6 reaxff
   fix 1 all acks2/reaxff 1 0.0 10.0 1.0e-6 param.acks2 maxiter 500

Description
"""""""""""

Perform the atom-condensed Kohn--Sham DFT to second order (ACKS2) charge
equilibration method as described in :ref:`(Verstraelen) <Verstraelen>`.
ACKS2 impedes unphysical long-range charge transfer sometimes seen with
QEq (e.g., for dissociation of molecules), at increased computational
cost.  It is typically used in conjunction with the ReaxFF force field
model as implemented in the :doc:`pair_style reaxff <pair_reaxff>`
command, but it can be used with any potential in LAMMPS, so long as it
defines and uses charges on each atom.  For more technical details about
the charge equilibration performed by fix acks2/reaxff, see the
:ref:`(O'Hearn) <O'Hearn>` paper.

The ACKS2 method minimizes the electrostatic energy of the system by
adjusting the partial charge on individual atoms based on interactions
with their neighbors.  It requires some parameters for each atom type.
If the *params* setting above is the word "reaxff", then these are
extracted from the :doc:`pair_style reaxff <pair_reaxff>` command and
the ReaxFF force field file it reads in.  If a file name is specified
for *params*, then the parameters are taken from the specified file
and the file must contain one line for each atom type.  The latter
form must be used when performing QeQ with a non-ReaxFF potential.
The lines should be formatted as follows:

.. parsed-literal::

   bond_softness
   itype chi eta gamma bcut

where the first line is the global parameter *bond_softness*. The
remaining 1 to Ntypes lines include *itype*, the atom type from 1 to
Ntypes, *chi*, the electronegativity in eV, *eta*, the self-Coulomb
potential in eV, *gamma*, the valence orbital exponent, and *bcut*, the
bond cutoff distance.  Note that these 4 quantities are also in the
ReaxFF potential file, except that eta is defined here as twice the eta
value in the ReaxFF file. Note that unlike the rest of LAMMPS, the units
of this fix are hard-coded to be :math:`\AA`, eV, and
electronic charge.

The optional *maxiter* keyword allows changing the max number
of iterations in the linear solver. The default value is 200.

.. versionchanged:: TBD

The optional *nowarn* keyword silences the warning message printed
when the maximum number of iterations was reached.  This can be
useful for simulations where the iteration limit is set
intentionally low.

.. versionchanged:: TBD

If the fix group does not contain all atoms, only the charges of the
atoms in the fix group are equilibrated.  The charges of all other atoms
are not changed, but they still contribute to the electrostatic
potential at the atoms in the fix group, i.e. they act as fixed external
charges.  All accelerated variants of this fix follow the same
convention.

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
iterations) for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

This fix is invoked during :doc:`energy minimization <minimize>`.

----------

.. include:: accel_styles.rst

----------

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
vector may only have components in non-periodic directions. Equal-style
variables can be used for electric field vector components without any further
settings. Atom-style variables can be used for spatially-varying electric field
vector components, but the resulting electric potential must be specified
as an atom-style variable using the *potential* keyword for `fix efield`.

The KOKKOS package version of this fix (*acks2/reaxff/kk*) requires
that the KOKKOS package is compiled for double precision
(``-D KOKKOS_PREC=double``, which is the default setting), since the
linear solver does not reliably converge with reduced precision.

Related commands
""""""""""""""""

:doc:`pair_style reaxff <pair_reaxff>`,
:doc:`fix qeq/reaxff <fix_qeq_reaxff>`,
:doc:`fix qtpie/reaxff <fix_qtpie_reaxff>`,
:doc:`fix qeq/rel/reaxff <fix_qeq_rel_reaxff>`

Default
"""""""

maxiter 200

----------

.. _O'Hearn:

**(O'Hearn)** O'Hearn, Alperen, Aktulga, SIAM J. Sci. Comput., 42(1), C1--C22 (2020).

.. _Verstraelen:

**(Verstraelen)** Verstraelen, Ayers, Speybroeck, Waroquier, J. Chem. Phys. 138, 074108 (2013).
