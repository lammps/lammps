.. index:: compute mliap

compute mliap command
====================

Syntax
""""""

.. code-block:: LAMMPS

   compute mliap

Examples
""""""""

.. code-block:: LAMMPS

   compute mliap model linear nelems nparams descriptor sna InP.mliap.descriptor

Description
"""""""""""

Compute style *mliap* provides a general interface to the gradient
of machine-learning interatomic potentials w.r.t. model parameters. 
It is used primarily for calculating the gradient of energy, force, and
stress components w.r.t. model parameters, which is useful when training
:doc:`MLIAP pair_style <pair_mliap>` to match target data.
It provides separate 
definitions of the interatomic potential functional form (*model*)
and the geometric quantities that characterize the atomic positions
(*descriptor*). By defining *model* and *descriptor* separately, 
it is possible to use many different models with a given descriptor,
or many different descriptors with a given model. Currently, the 
compute supports just two models, *linear* and *quadratic*,
and one descriptor, *sna*, the SNAP descriptor used by 
:doc:`pair_style snap <pair_snap>`, including the linear, quadratic,
and chem variants. Work is currently underway to extend
the interface to handle neural network energy models,
and it is also straightforward to add new descriptor styles.

The compute *mliap* command must be followed by two keywords
*model* and *descriptor* in either order. 

The *model* keyword is followed by a model style, currently limited to
either *linear* or *quadratic*. In both cases,
this is followed by two arguments. *nelems* is the number of elements.
It must be equal to the number of LAMMPS atom types. *nparams*
is the number of parameters per element for this model i.e.
the number of paramter gradients for each element. Note these definitions
are identical to those of *nelems* and *nparams* in the 
:doc:`pair_style mliap <pair_mliap>` model file.
 
The *descriptor* keyword is followed by a descriptor style, and additional arguments.
Currently the only descriptor style is *sna*, indicating the bispectrum component 
descriptors used by the Spectral Neighbor Analysis Potential (SNAP) potentials of 
:doc:`pair_style snap <pair_snap>`.
The \'p\' in SNAP is dropped, because keywords that match pair_styles are silently stripped 
out by the LAMMPS command parser. A single additional argument specifies the descriptor filename 
containing the parameters and setting used by the SNAP descriptor. 
The descriptor filename usually ends in the *.mliap.descriptor* extension.
The format of this file is identical to descriptor file in the 
:doc:`pair_style mliap <pair_mliap>`, and is described in detail
there. 

.. note::

The number of LAMMPS atom types (and the value of *nelems* in the model)
must match the value of *nelems* in the descriptor file. 

Compute *mliap* calculates a global array containing gradient information.
The number of columns in the array is the :math:`nelems \times nparams + 1`.
The first row of the array contain the derivative of potential energy w.r.t to
each parameter and each element. The last six rows
of the array contain the corresponding derivatives of the
virial stress tensor, listed in Voigt notation: *pxx*, *pyy*, *pzz*,
*pyz*, *pxz*, *pxy*. In between 3\*\ *N* rows containing the derivatives
of the force components. 

The element in the last column of each row contains
the potential energy, force, or stress, according to the row.
These quantities correspond to the user-specified reference potential
that must be subtracted from the target data when fitting SNAP.
The potential energy calculation uses the built in compute *thermo_pe*.
The stress calculation uses a compute called *snap_press* that is
automatically created behind the scenes, according to the following
command:

.. code-block:: LAMMPS

   compute snap_press all pressure NULL virial

See section below on output for a detailed explanation of the data
layout in the global array.

Atoms not in the group do not contribute to this compute. 
Neighbor atoms not in the group do not contribute to this compute.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e. each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently.

.. note::

   If you have a bonded system, then the settings of
   :doc:`special_bonds <special_bonds>` command can remove pairwise
   interactions between atoms in the same bond, angle, or dihedral.  This
   is the default setting for the :doc:`special_bonds <special_bonds>`
   command, and means those pairwise interactions do not appear in the
   neighbor list.  Because this fix uses the neighbor list, it also means
   those pairs will not be included in the calculation.  One way to get
   around this, is to write a dump file, and use the :doc:`rerun <rerun>`
   command to compute the bispectrum components for snapshots in the dump
   file.  The rerun script can use a :doc:`special_bonds <special_bonds>`
   command that includes all pairs in the neighbor list.

----------

**Output info:**

Compute *mliap* evaluates a global array.
The columns are arranged into
*nelems* blocks, listed in order of element *I*\ . Each block
contains one column for each of the *nparams* model parameters. 
A final column contains the corresponding energy, force component
on an atom, or virial stress component. The rows of the array appear
in the following order:

* 1 row: Derivatives of potential energy w.r.t. each parameter of each element.
* 3\*\ *N* rows: Derivatives of force components. x, y, and z components of 
force on atom *i* appearing in consecutive rows. The atoms are sorted based on atom ID.
* 6 rows: Derivatives of virial stress tensor  w.r.t. each parameter of each element.
The ordering of the rows follows Voigt notation: *pxx*, *pyy*, *pzz*,
*pyz*, *pxz*, *pxy*.

These values can be accessed by any command that uses per-atom values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

Restrictions
""""""""""""

This compute is part of the MLIAP package.  It is only enabled if
LAMMPS was built with that package.  In addition, building LAMMPS with the MLIAP package
requires building LAMMPS with the SNAP package.
See the :doc:`Build package <Build_package>` doc page for more info.
doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_style snap <pair_mliap>`

**Default:** none
