.. index:: fix electrode/conp
.. index:: fix electrode/conq
.. index:: fix electrode/thermo
.. index:: fix electrode/conp/intel
.. index:: fix electrode/conq/intel
.. index:: fix electrode/thermo/intel

fix electrode/conp command
==========================

Accelerator Variant: *electrode/conp/intel*

fix electrode/conq command
==========================

Accelerator Variant: *electrode/conq/intel*

fix electrode/thermo command
============================

Accelerator Variant: *electrode/thermo/intel*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID style args keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *electrode/conp* or *electrode/conq* or *electrode/thermo*
* args = arguments used by a particular style

  .. parsed-literal::

       *electrode/conp* args = potential eta
       *electrode/conq* args = charge eta
       *electrode/thermo* args = potential eta *temp* values
            potential = electrode potential
            charge = electrode charge
            eta = reciprocal width of electrode charge smearing
            *temp* values = T_v tau_v rng_v
                T_v = temperature of thermo-potentiostat
                tau_v = time constant of thermo-potentiostat
                rng_v = integer used to initialize random number generator

* zero or more keyword/value pairs may be appended
* keyword = *algo* or *symm* or *couple* or *etypes* or *ffield* or *write_mat* or *write_inv* or *read_mat* or *read_inv*

.. parsed-literal::

    *algo* values = *mat_inv* or *mat_cg* tol or *cg* tol
        specify the algorithm used to compute the electrode charges
    *symm* value = *on* or *off*
        turn on/off charge neutrality constraint for the electrodes
    *couple* values = group-ID val
        group-ID = group of atoms treated as additional electrode
        val = electric potential or charge on this electrode
    *etypes* value = *on* or *off*
        turn on/off type-based optimized neighbor lists (electrode and electrolyte types may not overlap)
    *ffield* value = *on* or *off*
        turn on/off finite-field implementation
    *write_mat* value = filename
        filename = file to which to write elastance matrix
    *write_inv* value = filename
        filename = file to which to write inverted matrix
    *read_mat* value = filename
        filename = file from which to read elastance matrix
    *read_inv* value = filename
        filename = file from which to read inverted matrix

Examples
""""""""

.. code-block:: LAMMPS

   fix fxconp bot electrode/conp -1.0 1.805 couple top 1.0 couple ref 0.0 write_inv inv.csv symm on
   fix fxconp electrodes electrode/conq 0.0 1.805 algo cg 1e-5
   fix fxconp bot electrode/thermo -1.0 1.805 temp 298 100 couple top 1.0

Description
"""""""""""

The *electrode* fixes implement the constant potential method (CPM)
(:ref:`Siepmann <Siepmann>`, :ref:`Reed <Reed3>`), and modern variants,
to accurately model electrified, conductive electrodes. This is
primarily useful for studying electrode-electrolyte interfaces,
especially at high potential differences or ionicities, with non-planar
electrodes such as nanostructures or nanopores, and to study dynamic
phenomena such as charging or discharging time scales or conductivity or
ionic diffusivities.

Each *electrode* fix allows users to set additional electrostatic
relationships between the specified groups which model useful
electrostatic configurations:

* *electrode/conp* sets potentials or potential differences between electrodes

  *  (resulting in changing electrode total charges)

* *electrode/conq* sets the total charge on each electrode

  *  (resulting in changing electrode potentials)

* *electrode/thermo* sets a thermopotentiostat
  :ref:`(Deissenbeck)<Deissenbeck>` between two electrodes

  * (resulting in changing charges and potentials with appropriate
     average potential difference and thermal variance)

The first group-ID provided to each fix specifies the first electrode
group, and more group(s) are added using the *couple* keyword for each
additional group.  While *electrode/thermo* only accepts two groups,
*electrode/conp* and *electrode/conq* accept any number of groups, up to
LAMMPS's internal restrictions (see Restrictions below). Electrode
groups must not overlap, i.e.  the fix will issue an error if any
particle is detected to belong to at least two electrode groups.

CPM involves updating charges on groups of electrode particles, per time
step, so that the system's total energy is minimized with respect to
those charges.  From basic electrostatics, this is equivalent to making
each group conductive, or imposing an equal electrostatic potential on
every particle in the same group (hence the name CPM).  The charges are
usually modelled as a Gaussian distribution to make the charge-charge
interaction matrix invertible (:ref:`Gingrich <Gingrich>`).  The keyword
*eta* specifies the distribution's width in units of inverse length.

.. versionadded:: 22Dec2022

Three algorithms are available to minimize the energy, varying in how
matrices are pre-calculated before a run to provide computational
speedup. These algorithms can be selected using the keyword *algo*:

* *algo mat_inv* pre-calculates the capacitance matrix and obtains the
  charge configuration in one matrix-vector calculation per time step

* *algo mat_cg* pre-calculates the elastance matrix (inverse of
  capacitance matrix) and obtains the charge configuration using a
  conjugate gradient solver in multiple matrix-vector calculations per
  time step

* *algo cg* does not perform any pre-calculation and obtains the charge
  configuration using a conjugate gradient solver and multiple
  calculations of the electric potential per time step.

For both *cg* methods, the command must specify the conjugate gradient
tolerance. *fix electrode/thermo* currently only supports the *mat_inv*
algorithm.

The keyword *symm* can be set *on* (or *off*) to turn on (or turn off)
the capacitance matrix constraint that sets total electrode charge to be
zero.  This has slightly different effects for each *fix electrode*
variant.  For *fix electrode/conp*, with *symm off*, the potentials
specified are absolute potentials, but the charge configurations
satisfying them may add up to an overall non-zero, varying charge for
the electrodes (and thus the simulation box). With *symm on*, the total
charge over all electrode groups is constrained to zero, and potential
differences rather than absolute potentials are the physically relevant
quantities.

For *fix electrode/conq*, with *symm off*, overall neutrality is
explicitly obeyed or violated by the user input (which is not
checked!). With *symm on*, overall neutrality is ensured by ignoring the
user-input charge for the last listed electrode (instead, its charge
will always be minus the total sum of all other electrode charges). For
*fix electrode/thermo*, overall neutrality is always automatically
imposed for any setting of *symm*, but *symm on* allows finite-field
mode (*ffield on*, described below) for faster simulations.

For all three fixes, any potential (or charge for *conq*) can be
specified as an equal-style variable prefixed with "v\_". For example,
the following code will ramp the potential difference between electrodes
from 0.0V to 2.0V over the course of the simulation:

.. code-block:: LAMMPS

   fix fxconp bot electrode/conp 0.0 1.805 couple top v_v symm on
   variable v equal ramp(0.0, 2.0)

Note that these fixes only parse their supplied variable name when
starting a run, and so these fixes will accept equal-style variables
defined *after* the fix definition, including variables dependent on the
fix's own output. This is useful, for example, in the fix's internal
finite-field commands (see below).  For an advanced example of this see
the in.conq2 input file in the directory
``examples/PACKAGES/electrode/graph-il``.

This fix necessitates the use of a long range solver that calculates and
provides the matrix of electrode-electrode interactions and a vector of
electrode-electrolyte interactions.  The Kspace styles
*ewald/electrode*, *pppm/electrode* and *pppm/electrode/intel* are
created specifically for this task :ref:`(Ahrens-Iwers) <Ahrens-Iwers>`.

For systems with non-periodic boundaries in one or two directions dipole
corrections are available with the :doc:`kspace_modify <kspace_modify>`.
For ewald/electrode a two-dimensional Ewald summation :ref:`(Hu) <Hu>`
can be used by setting "slab ew2d":

.. code-block:: LAMMPS

   kspace_modify slab <slab_factor>
   kspace_modify wire <wire_factor>
   kspace_modify slab ew2d

Two implementations for the calculation of the elastance matrix are
available with pppm and can be selected using the *amat onestep/twostep*
keyword.  *onestep* is the default; *twostep* can be faster for large
electrodes and a moderate mesh size but requires more memory.

.. code-block:: LAMMPS

   kspace_modify amat onestep/twostep

For all versions of the fix, the keyword-value *ffield on* enables the
finite-field mode (:ref:`Dufils <Dufils>`, :ref:`Tee <Tee>`), which uses
an electric field across a periodic cell instead of non-periodic
boundary conditions to impose a potential difference between the two
electrodes bounding the cell. The fix (with name *fix-ID*) detects which
of the two electrodes is "on top" (has the larger maximum *z*-coordinate
among all particles).  Assuming the first electrode group is on top, it
then issues the following commands internally:

.. code-block:: LAMMPS

   variable fix-ID_ffield_zfield equal (f_fix-ID[2]-f_fix-ID[1])/lz
   efield fix-ID_efield all efield 0.0 0.0 v_fix-ID_ffield_zfield

which implements the required electric field as the potential difference
divided by cell length.  The internal commands use variable so that the
electric field will correctly vary with changing potentials in the
correct way (for example with equal-style potential difference or with
*fix electrode/conq*).  This keyword requires two electrodes and will
issue an error with any other number of electrodes. This keyword
requires electroneutrality to be imposed (*symm on*) and will issue an
error otherwise.

.. versionchanged:: 22Dec2022

For all versions of the fix, the keyword-value *etypes on* enables
type-based optimized neighbor lists. With this feature enabled, LAMMPS
provides the fix with an occasional neighbor list restricted to
electrode-electrode interactions for calculating the electrode matrix,
and a perpetual neighbor list restricted to electrode-electrolyte
interactions for calculating the electrode potentials, using particle
types to list only desired interactions, and typically resulting in
5--10\% less computational time.  Without this feature the fix will
simply use the active pair style's neighbor list.  This feature cannot
be enabled if any electrode particle has the same type as any
electrolyte particle (which would be unusual in a typical simulation)
and the fix will issue an error in that case.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix currently does not write any information to restart files.

The *fix_modify tf* option enables the Thomas-Fermi metallicity model
(:ref:`Scalfi <Scalfi>`) and allows parameters to be set for each atom type.

.. code-block:: LAMMPS

   fix_modify ID tf type length voronoi


If this option is used parameters must be set for all atom types of the
electrode.

The *fix_modify timer* option turns on (off) additional timer outputs in the log
file, for code developers to track optimization.

.. code-block:: LAMMPS

   fix_modify ID timer on/off

----------

These fixes compute a global (extensive) scalar, a global (intensive)
vector, and a global array, which can be accessed by various
:doc:`output commands <Howto_output>`.

The global scalar outputs the energy added to the system by this fix,
which is the negative of the total charge on each electrode multiplied
by that electrode's potential.

The global vector outputs the potential on each electrode (and thus has
*N* entries if the fix manages *N* electrode groups), in :doc:`units
<units>` of electric field multiplied by distance (thus volts for *real*
and *metal* units).  The electrode groups' ordering follows the order in
which they were input in the fix command using *couple*. The global
vector output is useful for *fix electrode/conq* and *fix
electrode/thermo*, where potential is dynamically updated based on
electrolyte configuration instead of being directly set.

The global array has *N* rows and *2N+1* columns, where the fix manages
*N* electrode groups managed by the fix. For the *I*-th row of the
array, the elements are:

* array[I][1] = total charge that group *I* would have had *if it were
  at 0 V applied potential* * array[I][2 to *N* + 1] = the *N* entries
  of the *I*-th row of the electrode capacitance matrix (definition
  follows) * array[I][*N* + 2 to *2N* + 1] = the *N* entries of the
  *I*-th row of the electrode elastance matrix (the inverse of the
  electrode capacitance matrix)

The :math:`N \times N` electrode capacitance matrix, denoted :math:`\mathbf{C}`
in the following equation, summarizes how the total charge induced on each
electrode (:math:`\mathbf{Q}` as an *N*-vector) is related to the potential on
each electrode, :math:`\mathbf{V}`, and the charge-at-0V :math:`\mathbf{Q}_{0V}`
(which is influenced by the local electrolyte structure):

.. math::

   \mathbf{Q} = \mathbf{Q}_{0V} + \mathbf{C} \cdot \mathbf{V}

The charge-at-0V, electrode capacitance and elastance matrices are internally
used to calculate the potentials required to induce the specified total
electrode charges in *fix electrode/conq* and *fix electrode/thermo*. With the
*symm on* option, the electrode capacitance matrix would be singular, and thus
its last row is replaced with *N* copies of its top-left entry
(:math:`\mathbf{C}_{11}`) for invertibility.

The global array output is mainly useful for quickly determining the 'vacuum
capacitance' of the system (capacitance with only electrodes, no electrolyte),
and can also be used for advanced simulations setting the potential as some
function of the charge-at-0V (such as the ``in.conq2`` example mentioned above).

Please cite :ref:`(Ahrens-Iwers2022) <Ahrens-Iwers2>` in any publication that
uses this implementation.  Please cite also the publication on the combination
of the CPM with PPPM if you use *pppm/electrode* :ref:`(Ahrens-Iwers)
<Ahrens-Iwers>`.

----------

Restrictions
""""""""""""

For algorithms that use a matrix for the electrode-electrode
interactions, positions of electrode particles have to be immobilized at
all times.

With *ffield off* (i.e. the default), the box geometry is expected to be
*z*-non-periodic (i.e. *boundary p p f*), and this fix will issue an
error if the box is *z*-periodic. With *ffield on*, the box geometry is
expected to be *z*-periodic, and this fix will issue an error if the box
is *z*-non-periodic.

The parallelization for the fix works best if electrode atoms are evenly
distributed across processors. For a system with two electrodes at the bottom
and top of the cell this can be achieved with *processors * * 2*, or with the
line

.. code-block:: LAMMPS

   if "$(extract_setting(world_size) % 2) == 0" then "processors * * 2"

which avoids an error if the script is run on an odd number of
processors (such as on just one processor for testing).

The fix creates an additional group named *[fix-ID]_group* which is the
union of all electrode groups supplied to LAMMPS. This additional group
counts towards LAMMPS's limitation on the total number of groups
(currently 32), which may not allow scripts that use that many groups to
run with this fix.

The matrix-based algorithms (*algo mat_inv* and *algo mat_cg*) currently
store an interaction matrix (either elastance or capacitance) of *N* by
*N* doubles for each MPI process. This memory requirement may be
prohibitive for large electrode groups.  The fix will issue a warning if
it expects to use more than 0.5 GiB of memory.

Default
"""""""

The default keyword-option settings are *algo mat_inv*, *symm off*,
*etypes off* and *ffield off*.

----------

.. include:: accel_styles.rst

----------

.. _Siepmann:

**(Siepmann)** Siepmann and Sprik, J. Chem. Phys. 102, 511 (1995).

.. _Reed3:

**(Reed)** Reed *et al.*, J. Chem. Phys. 126, 084704 (2007).

.. _Deissenbeck:

**(Deissenbeck)** Deissenbeck *et al.*, Phys. Rev. Letters 126, 136803 (2021).

.. _Gingrich:

**(Gingrich)** Gingrich, `MSc thesis` <https://gingrich.chem.northwestern.edu/papers/ThesiswCorrections.pdf>` (2010).

.. _Ahrens-Iwers:

**(Ahrens-Iwers)** Ahrens-Iwers and Meissner, J. Chem. Phys. 155, 104104 (2021).

.. _Hu:

**(Hu)** Hu, J. Chem. Theory Comput. 10, 5254 (2014).

.. _Dufils:

**(Dufils)** Dufils *et al.*, Phys. Rev. Letters 123, 195501 (2019).

.. _Tee:

**(Tee)** Tee and Searles, J. Chem. Phys. 156, 184101 (2022).

.. _Scalfi:

**(Scalfi)** Scalfi *et al.*, J. Chem. Phys., 153, 174704 (2020).

.. _Ahrens-Iwers2:

**(Ahrens-Iwers2022)** Ahrens-Iwers *et al.*, J. Chem. Phys. 157, 084801 (2022).
