.. index:: compute mliap

compute mliap command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID mliap ... keyword values ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* mliap = style name of this compute command
* two or more keyword/value pairs must be appended
* keyword = *model* or *descriptor* or *gradgradflag*

  .. parsed-literal::

       *model* values = style
         style = *linear* or *quadratic* or *mliappy*
       *descriptor* values = style filename
         style = *sna*
         filename = name of file containing descriptor definitions
       *gradgradflag* value = 0/1
         toggle gradgrad method for force gradient

Examples
""""""""

.. code-block:: LAMMPS

   compute mliap model linear descriptor sna Ta06A.mliap.descriptor

Description
"""""""""""

Compute style *mliap* provides a general interface to the gradient
of machine-learning interatomic potentials with respect to model parameters.
It is used primarily for calculating the gradient of energy, force, and
stress components with respect to model parameters, which is useful when
training :doc:`mliap pair_style <pair_mliap>` models to match target data.
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

The *model* keyword is followed by the model style (*linear*,
*quadratic* or *mliappy*).  The *mliappy* model is only available if
LAMMPS is built with the *mliappy* Python module. There are
:ref:`specific installation instructions <mliap>` for that module.

The *descriptor* keyword is followed by a descriptor style, and
additional arguments.  The compute currently supports two descriptor
styles *sna* and *so3*, but it is is straightforward to add additional
descriptor styles.  The SNAP descriptor style *sna* is the same as that
used by :doc:`pair_style snap <pair_snap>`, including the linear,
quadratic, and chem variants.  A single additional argument specifies
the descriptor filename containing the parameters and setting used by
the SNAP descriptor.  The descriptor filename usually ends in the
*.mliap.descriptor* extension.  The format of this file is identical to
the descriptor file in the :doc:`pair_style mliap <pair_mliap>`, and is
described in detail there.

.. note::

   The number of LAMMPS atom types (and the value of *nelems* in the model)
   must match the value of *nelems* in the descriptor file.

Compute *mliap* calculates a global array containing gradient information.
The number of columns in the array is *nelems* :math:`\times` *nparams* + 1.
The first row of the array contain the derivative of potential energy with
respect to. to each parameter and each element. The last six rows
of the array contain the corresponding derivatives of the
virial stress tensor, listed in Voigt notation: *pxx*, *pyy*, *pzz*,
*pyz*, *pxz*, and *pxy*. In between the energy and stress rows are
the :math:`3N` rows containing the derivatives of the force components.
See section below on output for a detailed description of how
rows and columns are ordered.

The element in the last column of each row contains
the potential energy, force, or stress, according to the row.
These quantities correspond to the user-specified reference potential
that must be subtracted from the target data when training a model.
The potential energy calculation uses the built in compute *thermo_pe*.
The stress calculation uses a compute called *mliap_press* that is
automatically created behind the scenes, according to the following
command:

.. code-block:: LAMMPS

   compute mliap_press all pressure NULL virial

See section below on output for a detailed explanation of the data
layout in the global array.

The optional keyword *gradgradflag* controls how the force
gradient is calculated. A value of 1 requires that the model provide
the matrix of double gradients of energy with respect to both parameters
and descriptors. For the linear and quadratic models this matrix is
sparse and so is easily calculated and stored. For other models, this
matrix may be prohibitively expensive to calculate and store.
A value of 0 requires that the descriptor provide the derivative
of the descriptors with respect to the position of every neighbor atom.
This is not optimal for linear and quadratic models, but may be
a better choice for more complex models.

Atoms not in the group do not contribute to this compute.
Neighbor atoms not in the group do not contribute to this compute.
The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e., each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently.

.. note::

   If the user-specified reference potentials includes bonded and
   non-bonded pairwise interactions, then the settings of
   :doc:`special_bonds <special_bonds>` command can remove pairwise
   interactions between atoms in the same bond, angle, or dihedral.  This
   is the default setting for the :doc:`special_bonds <special_bonds>`
   command, and means those pairwise interactions do not appear in the
   neighbor list.  Because this fix uses the neighbor list, it also means
   those pairs will not be included in the calculation. The :doc:`rerun <rerun>`
   command is not an option here, since the reference potential is required
   for the last column of the global array. A work-around is to prevent
   pairwise interactions from being removed by explicitly adding a
   *tiny* positive value for every pairwise interaction that would otherwise be
   set to zero in the :doc:`special_bonds <special_bonds>` command.

----------

Output info
"""""""""""

Compute *mliap* evaluates a global array.  The columns are arranged into
*nelems* blocks, listed in order of element *I*\ . Each block
contains one column for each of the *nparams* model parameters.
A final column contains the corresponding energy, force component
on an atom, or virial stress component. The rows of the array appear
in the following order:

* 1 row: Derivatives of potential energy with respect to each parameter of each element.
* :math:`3N` rows: Derivatives of force components; the *x*, *y*, and *z*
  components of the force on atom *i* appear in consecutive rows. The atoms are
  sorted based on atom ID.
* 6 rows: Derivatives of the virial stress tensor with respect to each
  parameter of each element. The ordering of the rows follows Voigt notation:
  *pxx*, *pyy*, *pzz*, *pyz*, *pxz*, *pxy*.

These values can be accessed by any command that uses a global array
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options. To see how this command
can be used within a Python workflow to train machine-learning interatomic
potentials, see the examples in `FitSNAP <https://github.com/FitSNAP/FitSNAP>`_.

Restrictions
""""""""""""

This compute is part of the ML-IAP package.  It is only enabled if
LAMMPS was built with that package. In addition, building LAMMPS with
the ML-IAP package requires building LAMMPS with the ML-SNAP package.
The *mliappy* model also requires building LAMMPS with the PYTHON
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`pair_style mliap <pair_mliap>`

Default
"""""""

The keyword defaults are gradgradflag = 1
