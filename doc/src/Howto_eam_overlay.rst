Two-band (and multi-band) EAM potentials with hybrid/overlay
============================================================

Some embedded-atom-method (EAM) potentials add a *second* embedding term to
the standard EAM energy expression in order to capture effects that a single
density and embedding function cannot.  The best known examples are the
*two-band* (or *double-band*, "2BM") models for Fe-Cr, which add an
"s-band" embedding term on top of the usual "d-band" EAM so that the sign of
the alloy mixing enthalpy changes with chromium concentration
:ref:`(Olsson) <Olsson2band>`, :ref:`(Bonny) <Bonny2band>`.

This Howto shows that such a model can be run with the existing
:doc:`pair_style eam/fs <pair_eam>` and
:doc:`pair_style hybrid/overlay <pair_hybrid>` commands, **without modifying
or recompiling LAMMPS**, and how to build the required potential files.

The two-band energy expression
-------------------------------

In a two-band model the energy of atom *i* is

.. math::

   E_i = \frac{1}{2} \sum_{j \neq i} V(r_{ij})
         + F^{d}\!\left(\rho^{d}_i\right)
         + F^{s}\!\left(\rho^{s}_i\right), \qquad
   \rho^{d}_i = \sum_{j \neq i} \phi^{d}(r_{ij}), \quad
   \rho^{s}_i = \sum_{j \neq i} \phi^{s}(r_{ij}) .

This is an ordinary EAM (pair term :math:`V`, d-band density
:math:`\phi^{d}`, d-band embedding :math:`F^{d}`) plus one *extra* embedding
term :math:`F^{s}` acting on a second, independent density
:math:`\rho^{s}`.  Grouping the terms,

.. math::

   E_i = \underbrace{\left[\,\frac{1}{2} \sum_{j} V + F^{d}(\rho^{d})\,\right]}_{\text{a standard eam/fs}}
       + \underbrace{\left[\,F^{s}(\rho^{s})\,\right]}_{\text{eam/fs with }V\equiv 0}

shows that the total energy is the sum of two EAM potentials.  Adding
energies, forces, and virials of two potentials acting on the same atoms is
exactly what :doc:`pair_style hybrid/overlay <pair_hybrid>` does.

Why the s-band file must be *eam/fs*
------------------------------------

In the two-band models the second (s-band) density is non-zero only between
*unlike* atoms; it vanishes between atoms of the same element
:ref:`(Bonny) <Bonny2band>`:

.. math::

   \phi^{s}_{AA} = \phi^{s}_{BB} = 0, \qquad \phi^{s}_{AB} = \phi^{s}_{BA} \neq 0 .

A density that depends on the element of *both* the source and the target
atom can only be expressed in the Finnis-Sinclair form, i.e. with
:doc:`pair_style eam/fs <pair_eam>`, whose density functions
:math:`\rho_{\alpha\beta}(r)` are indexed per element pair.  It *cannot* be
written as an :doc:`eam/alloy <pair_eam>` potential, which has a single
density per element summed over all neighbors.  Therefore the s-band file is
always an *eam/fs* file in which the like-element density functions are set
to zero and only the unlike-element density is non-zero.  Its pair-potential
(:math:`r\,\phi`) section is set entirely to zero, because the pair term
belongs to the d-band file.

The d-band file is an ordinary binary EAM (the pair term plus the d-band
density and embedding).  Because its density depends only on the neighbor's
element it could equally be written as *eam/alloy*; using *eam/fs* for both
files keeps the recipe uniform.

The recipe
----------

.. code-block:: LAMMPS

   pair_style hybrid/overlay eam/fs eam/fs
   pair_coeff * * eam/fs 1 FeCr.dband.eam.fs Fe Cr   # pair V + d-band F_d(rho_d)
   pair_coeff * * eam/fs 2 FeCr.sband.eam.fs Fe Cr   # s-band F_s(rho_s), no pair

The numeric ``1`` and ``2`` after ``eam/fs`` are required because the same
sub-style is used twice; see :doc:`pair_hybrid <pair_hybrid>`.  Each
*eam/fs* instance keeps its own density and embedding-derivative arrays and
performs its own communication, so the two densities
:math:`\rho^{d}` and :math:`\rho^{s}` never mix.  ``hybrid/overlay`` sums the
two contributions, and because the s-band file has no pair term there is no
double counting: the combined energy is exactly
:math:`\tfrac{1}{2}\sum V + F^{d}(\rho^{d}) + F^{s}(\rho^{s})`.

The two files are independent *eam/fs* files and may use different radial
grids and cutoffs; the neighbor list is built with the larger of the two
cutoffs.  Because both files store the same element masses, you may also set
masses explicitly with the :doc:`mass <mass>` command to avoid any ambiguity.

This generalizes to more than two bands: overlay *N* *eam/fs* files to add
*N* embedding terms.  The same construction works for any model of the form
"standard EAM plus one or more additive embedding terms."

Building the s-band file
------------------------

The helper script ``tools/eam_generate/two_band_fecr.py`` writes the s-band
*eam/fs* file (cross-only density, :math:`F^{s}` embedding, zero pair term)
and can also emit an illustrative demo d-band file plus a ready-to-run input
deck:

.. code-block:: bash

   cd tools/eam_generate
   # self-contained runnable demo (writes both files + a sample input):
   ./two_band_fecr.py --demo --prefix FeCr_demo

   # real use: add an s-band to your own trusted Fe-Cr d-band eam/fs file:
   ./two_band_fecr.py --dband MyFeCr.eam.fs --preset olsson-vasp --prefix FeCr_2bm

The script implements the s-band density as the square of a 4s Slater
function with a smooth cutoff and the s-band embedding as
:math:`F^{s}(\rho)=c_1\sqrt{\rho}+c_2\rho^2+c_3\rho^4`, matching the forms
used in the published models.  It ships parameter presets taken from the
literature, but note that the absolute normalization of :math:`\rho^{s}` and
the matching :math:`F^{s}` coefficients are convention-sensitive: validate a
preset against a known result (for instance the concentration dependence of
the mixing enthalpy) before using it for production work.  The default
``--demo`` parameters are self-consistent and illustrative only; they
exercise the overlay mechanics but are not a published fit.

A worked check
--------------

Running the demo input produced by ``--demo`` and decomposing the energy with
:doc:`compute pair <compute_pair>` confirms the construction:

.. code-block:: LAMMPS

   compute ed all pair eam/fs 1 epair   # d-band sub-style energy
   compute es all pair eam/fs 2 epair   # s-band sub-style energy
   thermo_style custom step pe c_ed c_es

The total potential energy is exactly ``c_ed + c_es``.  For a *pure* element
the s-band energy ``c_es`` is identically zero, because the like-element
s-density is zero, so :math:`\rho^{s}=0` and :math:`F^{s}(0)=0`.  Introducing
unlike neighbors (an alloy) makes ``c_es`` non-zero: the s-band term samples
the local concentration, which is the entire purpose of the model.

Restrictions and notes
----------------------

* The s-band embedding of these models typically has a :math:`\sqrt{\rho}`
  term, whose slope diverges as :math:`\rho \to 0`.  This is harmless here:
  LAMMPS tabulates :math:`F^{s}` as a spline (finite at every grid point) and
  the force contribution is multiplied by the s-density derivative, which is
  identically zero for like pairs.  A pure-element atom therefore contributes
  no spurious s-band force.

* This Howto reproduces the two-band *functional form* exactly.  Reproducing a
  specific published *parameterization* additionally requires the correct
  d-band/pair files (for Fe-Cr, the Mendelev iron potential plus the refit
  chromium and cross interactions) and validated s-band coefficients.

----------

.. _Olsson2band:

**(Olsson)** Olsson, Wallenius, Domain, Nordlund, Malerba, Physical Review B,
72, 214119 (2005).

.. _Bonny2band:

**(Bonny)** Bonny, Pasianot, Terentyev, Malerba, Philosophical Magazine, 91,
1724 (2011).
