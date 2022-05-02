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

.. parsed-literal::

   fix ID group-ID electrode/conp potential eta keyword values ...
   fix ID group-ID electrode/conq charge eta keyword values ...
   fix ID group-ID electrode/thermo potential eta temp T_v tau_v rng_v keyword values ...

* ID, group-ID are documented in fix command
* mode = electrode/conp or electrode/conq or electrode/thermo
* potential = electrode potential
* charge = electrode charge
* eta = reciprocal width of electrode charge smearing
* T_v = temperature of thermo-potentiostat
* tau_v = time constant of thermo-potentiostat
* rng_v = integer used to initialize random number generator

.. parsed-literal::

    *symm(etry) on/off*
        turn on/off charge neutrality constraint for the electrodes
    *couple group-ID value*
        group-ID = group of atoms treated as additional electrode
        value = electric potential or charge on this electrode
    *etypes values = atom types*
        specify atom types exclusive to the electrode for optimized neighbor lists
    *ffield on/off*
        turn on/off finite-field implementation
    *write_mat filename*
        write elastance matrix to file
    *write_inv filename*
        write inverted matrix to file
    *read_mat filename*
        read elastance matrix from file
    *read_inv filename*
        read inverted matrix from file

Examples
""""""""

.. code-block:: LAMMPS

   fix fxconp bot electrode/conp -1.0 1.805 couple top 1.0 couple ref 0.0 write_inv inv.csv symm on
   fix fxconp electrodes electrode/conq 0.0 1.805
   fix fxconp bot electrode/thermo -1.0 1.805 temp 298 100 couple top 1.0

Description
"""""""""""

fix electrode/conp mode implements a constant potential method (CPM)
(:ref:`Siepmann <Siepmann>`, :ref:`Reed <Reed3>`). Charges of groups specified
via group-ID and optionally with the `couple` keyword are adapted to meet their respective
potential at every time step.  An arbitrary number of electrodes can be set but
the respective groups may not overlap. Electrode charges have a Gaussian charge
distribution with reciprocal width eta. The energy minimization is achieved via
matrix inversion :ref:`(Wang) <Wang5>`.

fix electrode/conq enforces a total charge specified in the input on each electrode. The energy is
minimized w.r.t. the charge distribution within the electrode.

fix electrode/thermo implements a thermo-potentiostat :ref:`(Deissenbeck)
<Deissenbeck>`. Temperature and time constant of the thermo-potentiostat need
to be specified using the temp keyword. Currently, only two electrodes are possible with
this style.

This fix necessitates the use of a long range solver that calculates and provides the matrix
of electrode-electrode interactions and a vector of electrode-electrolyte
interactions.  The Kspace styles *ewald/electrode*, *pppm/electrode* and
*pppm/electrode/intel* are created specifically for this task
:ref:`(Ahrens-Iwers) <Ahrens-Iwers>`.

For systems with non-periodic boundaries in one or two directions dipole
corrections are available with the :doc:`kspace_modify <kspace_modify>`.  For
ewald/electrode a two-dimensional Ewald summation :ref:`(Hu) <Hu>` can be used
by setting "slab ew2d":

.. code-block:: LAMMPS

   kspace_modify slab <slab_factor>
   kspace_modify wire <wire_factor>
   kspace_modify slab ew2d

Two implementations for the calculation of the elastance matrix are available
with pppm and can be selected using the *amat onestep/twostep* keyword.
*onestep* is the default; *twostep* can be faster for large electrodes and a
moderate mesh size but requires more memory.

.. code-block:: LAMMPS

   kspace_modify amat onestep/twostep


The *fix_modify tf* option enables the Thomas-Fermi metallicity model
(:ref:`Scalfi <Scalfi>`) and allows parameters to be set for each atom type.

.. code-block:: LAMMPS

   fix_modify ID tf type length voronoi


If this option is used parameters must be set for all atom types of the electrode.

The *fix_modify timer* option turns on (off) additional timer outputs in the log
file, for code developers to track optimization.

.. code-block:: LAMMPS

   fix_modify ID timer on/off

----------

These fixes compute a global array which can be accessed by various
:doc:`output commands <Howto_output>`. The array has *N* rows and *2N+2*
columns, where *N* is the number of electrode groups managed by the fix. For the
*I*-th row of the array, the elements are:

* array[I][1] = potential instantaneously applied to group *I*
* array[I][1] = total charge that group *I* would have had *if it were at 0 V
  applied potential*
* array[I][2 to *N* + 1] = the *N* entries of the *I*-th row of the electrode
  capacitance matrix (definition follows)
* array[I][*N* + 2 to *2N* + 1] = the *N* entries of the *I*-th row of the electrode
  elastance matrix (the inverse of the electrode capacitance matrix)

The first column is handy for *fix electrode/conq* and *fix
electrode/thermo* simulations, to output the potentials dynamically calculated and
imposed by the fix (see following formula). The potential is in units of electric field times
distance, which gives volts for most (but not all) :doc:`units <units>` settings.

The :math:`N \times N` electrode capacitance matrix, denoted :math:`\mathbf{C}`
in the following equation, summarizes how the total charge induced on each electrode
(:math:`\mathbf{Q}` as an *N*-vector) is related to the potential on each
electrode, :math:`\mathbf{V}`, and the charge-at-0V :math:`\mathbf{Q}_{0V}`
(which is influenced by the local electrolyte structure):

.. math::

   \mathbf{Q} = \mathbf{Q}_{0V} + \mathbf{C} \cdot \mathbf{V}

The charge-at-0V, electrode capacitance and elastance matrices are internally used to
calculate the potentials required to induce the specified total electrode charges
in *fix electrode/conq* and *fix electrode/thermo*. With the *symm on* option,
the electrode capacitance matrix would be singular, and thus its last row is
replaced with *N* copies of its top-left entry (:math:`\mathbf{C}_{11}`) for
invertibility.

----------

Restrictions
""""""""""""

Positions of electrode particles have to be immobilized at all times.

The parallelization for the fix works best if electrode atoms are evenly
distributed across processors. For a system with two electrodes at the bottom
and top of the cell this can be achieved with *processors * * 2*, or with the
line

.. code-block:: LAMMPS

   if "$(extract_setting(world_size) % 2) == 0" then "processors * * 2"

which avoids an error if the script is run on an odd number of processors (such
as on just one processor for testing).

----------

.. include:: accel_styles.rst

----------

.. _Siepmann:

**(Siepmann)** Siepmann and Sprik, J. Chem. Phys. 102, 511 (1995).

.. _Reed3:

**(Reed)** Reed *et al.*, J. Chem. Phys. 126, 084704 (2007).

.. _Wang5:

**(Wang)** Wang *et al.*, J. Chem. Phys. 141, 184102 (2014).

.. _Deissenbeck:

**(Deissenbeck)** Deissenbeck *et al.*, Phys. Rev. Letters 126, 136803 (2021).

.. _Ahrens-Iwers:

**(Ahrens-Iwers)** Ahrens-Iwers and Meissner, J. Chem. Phys. 155, 104104 (2021).

.. _Hu:

**(Hu)** Hu, J. Chem. Theory Comput. 10, 5254 (2014).

.. _Scalfi:

**(Scalfi)** Scalfi *et al.*, J. Chem. Phys., 153, 174704 (2020).

