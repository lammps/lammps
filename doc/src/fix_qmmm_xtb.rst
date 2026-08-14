.. index:: fix qmmm/xtb

fix qmmm/xtb command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID qmmm/xtb elements E1 E2 ... Entypes cutoff distance keyword value ...

* ID and group-ID are documented in :doc:`fix <fix>`
* qmmm/xtb = style name of this fix command
* atoms in group-ID form the QM region; all other atoms form the MM region
* elements = one chemical symbol for every LAMMPS atom type
* cutoff = explicit QM/MM point-charge cutoff (distance units)
* zero or more keyword/value pairs may be appended

  .. parsed-literal::

       *charge* value = integer total charge of the QM region (default 0)
       *uhf* value = integer number of unpaired electrons (default 0)
       *method* value = gfn1 or gfn2 (default gfn2)
       *accuracy* value = xTB SCC accuracy multiplier (default 1.0e-3)
       *maxiter* value = maximum number of SCC iterations (default 250)
       *etemp* value = electronic temperature in K (default 300.0)
       *mmhardness* value = external-charge hardness in atomic units (default 0.0)
       *ewald_alpha* value = direct QM-image Ewald coefficient in 1/Angstrom
       *kmax* values = kmax_x kmax_y kmax_z for the direct QM-image Ewald sum
       *ksqmax* value = maximum squared integer k-vector for the direct sum

Examples
""""""""

.. code-block:: LAMMPS

   group qm id 1:12
   pair_style hybrid/overlay lj/cut 10.0 coul/long 8.0
   pair_coeff * * lj/cut 0.10 3.0
   pair_coeff * * coul/long
   kspace_style pppm/xtb 1.0e-7
   fix qmmm qm qmmm/xtb elements H O Na Cl cutoff 8.0 charge 0 uhf 0

Description
"""""""""""

.. versionadded:: TBD

This fix performs a non-covalent, electrostatically embedded GFN1-xTB or
GFN2-xTB QM/MM calculation in a three-dimensional periodic cell.  The
*method* keyword selects the Hamiltonian.  Atoms in the fix group are treated
by xTB.  MM charge sites within *cutoff* of any QM atom are passed to xTB as
explicit external point charges, including the analytic gradients on those
sites.

The periodic electrostatic partition combines explicit near-field MM charges,
a PPPM correction for the remaining MM potential, and a direct-Ewald response
for periodic QM images.  During every xTB SCC iteration the atom-dependent
Hamiltonian shift is

.. math::

   \phi_i^\mathrm{shift} =
   \left(\phi_i^\mathrm{MM,PPPM} - \phi_i^\mathrm{MM,near}\right)
   + \sum_j A_{ij}^\mathrm{QM,image} q_j.

The first term is fixed during an MD force evaluation.  LAMMPS obtains
the reciprocal MM potential from PPPM at the QM atom positions, then
subtracts the Gaussian-screened contribution of the explicit near MM
charges and adds the neutralizing-background term.  The second term is
a dense direct-Ewald response of the usually small QM region and is
updated from the current xTB Mulliken charges at every SCC iteration.

After SCC convergence, the fix stores the Mulliken charges in the
LAMMPS per-atom charge array for the QM atoms.  It retains the QM and
near-MM analytic gradients returned by xTB, adds the periodic QM-image
and QM/MM corrections, and removes the classical real-space and
reciprocal terms that would otherwise be counted twice.  The underlying
MM-MM force-field contribution is unchanged.

The scalar computed by this fix is an extensive energy correction.  It
is designed to be added to the ordinary LAMMPS potential energy; it is
not the isolated xTB energy.  The fix also contributes its current
first-stage virial estimate to global pressure.  Constant-pressure or
changing-box simulations are rejected because the direct QM-image
virial has not yet been validated against box finite differences.

The *elements* list maps every LAMMPS atom type to an atomic number.
``X`` or ``NULL`` may be used for an MM-only extra point-charge site.
Such a type cannot occur in the QM group.

The *mmhardness* behavior is:

* a value greater than :math:`10^{-6}` assigns that value to every MM
  point charge;
* a value between :math:`-10^{-6}` and :math:`10^{-6}` assigns the
  selected xTB method's hydrogen hardness to every MM point charge;
* a value less than :math:`-10^{-6}` uses the selected xTB method's chemical
  hardness of each MM element multiplied by ``abs(mmhardness)``.  Extra and
  virtual point sites use hydrogen's hardness.

When *ewald_alpha* is omitted it is set to
:math:`10/V^{1/3}`; *kmax* defaults to ``8 8 8`` and *ksqmax* defaults
to 64.  These parameters should be converged for the chosen cell and QM
region independently of the PPPM accuracy.

Classical force-field setup
"""""""""""""""""""""""""""

With ``kspace_style pppm/xtb``, the fix requires a distinct Coulomb-only
:doc:`pair_style coul/long <pair_coul>` sub-style whose cutoff equals
*cutoff*.  With ``kspace_style pppm/tip4p/xtb``, either the Coulomb-only
``pair_style tip4p/long`` or the standard combined
:doc:`pair_style lj/cut/tip4p/long <pair_lj_cut_tip4p>` may be used.  In the
combined style, the two reference evaluations have identical Lennard-Jones
terms, so those terms cancel from the QM/MM correction and the original MM
Lennard-Jones force and energy are preserved.

When the QM and MM regions use disjoint atom types, a Coulomb-only sub-style
can instead be assigned to MM-MM type pairs only.  For example, if types 1-5
are QM and types 6-7 are MM, use
``pair_coeff 6*7 6*7 tip4p/long`` and assign a separate ``lj/cut`` sub-style
to the required Lennard-Jones pairs.  The fix verifies the type partition and
the ``pair_style hybrid/overlay`` mappings before omitting the two pair
reference evaluations.  A type shared by the QM and MM regions or any other
partial Coulomb mapping is rejected.

For implicit TIP4P, the oxygen charge is embedded at the virtual M site and
the returned QM/MM force is redistributed to O/H/H with the same geometry
weights as the TIP4P pair and KSpace styles.  An implicit TIP4P molecule must
remain entirely in the MM region.  Water in the QM region must instead be
represented by ordinary real O/H atoms without the implicit TIP4P virtual-site
mapping.

The user is responsible for removing classical bonded, angle,
dihedral, improper, and short-range non-bonded terms that are replaced by
the QM Hamiltonian.  The initial implementation assumes no covalent
QM/MM boundary and has no link atoms.  QM/MM pairs must not have
``special_bonds`` Coulomb scaling, since xTB's explicit external-charge
embedding does not apply LAMMPS topology scaling.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information is written to binary restart files.  Reissue the fix
command after :doc:`read_restart <read_restart>`; xTB starts a new SCC
wave function and then reuses it between subsequent force evaluations.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by the
standard fix energy mechanism.  No fix-specific *fix_modify* keywords
are defined.

This fix computes a global scalar, accessible as ``f_ID``, containing
the extensive QM/MM energy correction.  It also updates ``q`` for QM
atoms with the converged Mulliken charges.  No per-atom energy, per-atom
virial, global vector, or per-atom array is provided.

Fixed-box atom-coordinate minimization is supported.  The xTB SCC solve,
MM potential projection, and periodic corrections are recomputed at every
trial geometry requested by the minimizer.  The QM/MM correction energy is
included in the minimizer objective by default.  Do not use
``fix_modify ID energy no`` during minimization; the fix rejects that
inconsistent configuration.

Energy minimization is more sensitive to numerical noise than dynamics.
Converge the xTB *accuracy*, PPPM accuracy, and direct-Ewald parameters for
the requested force tolerance.  Also avoid placing MM sites close to the
explicit *cutoff*, since changing the explicit point-charge set during a line
search can make the numerical objective less smooth.  Cell optimization with
``fix box/relax`` remains unsupported.

Restrictions
""""""""""""

This fix is part of the QMMM-XTB package.  It is available only when
LAMMPS is built with that package and compatible libxtb 6.7 or newer Fortran
modules; see the :ref:`QMMM-XTB build instructions <qmmm-xtb>`.

The initial implementation has these restrictions:

* GFN1-xTB and GFN2-xTB are supported;
* only ``real`` and ``metal`` units are supported;
* the global :doc:`dielectric <dielectric>` setting must remain 1.0;
* orthogonal and triclinic boxes are supported, but the box must remain
  fixed and periodic in all three dimensions;
* :doc:`kspace_style pppm/xtb <kspace_style>` and ``pppm/tip4p/xtb`` are
  supported; other PPPM variants, including dispersion, slab, and wire
  styles, are not supported;
* the QM group and its atom IDs must remain fixed and the QM region must
  be compact relative to the periodic cell;
* covalent QM/MM boundaries and link atoms are not supported;
* fixed-box atom-coordinate minimization is supported, but r-RESPA and
  ``fix box/relax`` are not;
* constant-pressure fixes and changing-box fixes are not supported.

Related commands
""""""""""""""""

:doc:`fix qmmm <fix_qmmm>`, :doc:`fix mdi/qmmm <fix_mdi_qmmm>`,
:doc:`pair_style coul/long <pair_coul>`, ``pair_style tip4p/long``,
:doc:`pair_style lj/cut/tip4p/long <pair_lj_cut_tip4p>`,
:doc:`kspace_style pppm/xtb <kspace_style>`, ``kspace_style pppm/tip4p/xtb``

Default
"""""""

The optional keyword defaults are ``method gfn2``, ``charge 0``, ``uhf 0``,
``accuracy 1.0e-3``, ``maxiter 250``, ``etemp 300.0``,
``mmhardness 0.0``, ``kmax 8 8 8``, and ``ksqmax 64``.  The direct
QM-image Ewald coefficient defaults to :math:`10/V^{1/3}`.
