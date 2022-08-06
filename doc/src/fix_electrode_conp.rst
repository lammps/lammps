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
* keyword = *symm* or *couple* or *etypes* or *ffield* or *write_mat* or *write_inv* or *read_mat* or *read_inv*

.. parsed-literal::

    *symm* value = *on* or *off*
        turn on/off charge neutrality constraint for the electrodes
    *couple* values = group-ID val
        group-ID = group of atoms treated as additional electrode
        val = electric potential or charge on this electrode
    *etypes* values = type
        type = atom type (can be a range) exclusive to the electrode for optimized neighbor lists
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

The *fix_modify set* options allow calculated quantities to be accessed via
internal variables. Currently four types of quantities can be accessed:

.. code-block:: LAMMPS

   fix-modify ID set v group-ID variablename
   fix-modify ID set qsb group-ID variablename
   fix-modify ID set mc group-ID1 group-ID2 variablename
   fix-modify ID set me group-ID1 group-ID2 variablename

One use case is to output the potential that is internally calculated and
applied to each electrode group by *fix electrode/conq* or *fix electrode/thermo*.
For that case the *v* option makes *fix electrode* update the variable
*variablename* with the potential applied to group *group-ID*, where *group-ID*
must be a group whose charges are updated by *fix electrode* and *variablename*
must be an internal-style variable:

.. code-block:: LAMMPS

   fix conq bot electrode/conq -1.0 1.979 couple top 1.0
   variable vbot internal 0.0
   fix_modify conq set v bot vbot

The *qsb* option similarly outputs the total updated charge of the group if its
potential were 0.0V. The *mc* option requires two *group-IDs*, and outputs the
entry \{*group-ID1*, *group-ID2*\} of the (symmetric) *macro-capacitance* matrix
(MC) which relates the electrodes' applied potentials (V), total charges (Q), and
total charges at 0.0 V (Qsb):

.. math::

   \mathbf{Q} = \mathbf{Q}_{SB} + \mathbf{MC} \cdot \mathbf{V}

Lastly, the *me* option also requires two *group-IDs* and outputs the entry
\{*group-ID1*, *group-ID2*\} of the *macro-elastance* matrix, which is the
inverse of the macro-capacitance matrix. (As the names denote, the
macro-capacitance matrix gives electrode charges from potentials, and the
macro-elastance matrix gives electrode potentials from charges).

.. warning::

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

