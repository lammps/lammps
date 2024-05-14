.. index:: pair_style metatensor

pair_style metatensor command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style metatensor model_path ... keyword values ...

* model_path = path to the file containing the exported metatensor model
* keyword = *device* or *extensions* or *check_consistency*

  .. parsed-literal::

       *device* values = device_name
         device_name = name of the Torch device to use for the calculations
       *extensions* values = directory
         directory = path to a directory containing TorchScript extensions as
         shared libraries. If the model uses extensions, we will try to load
         them from this directory first
       *check_consistency* values = on or off
         set this to on/off to enable/disable internal consistency checks,
         verifying both the data passed by LAMMPS to the model, and the data
         returned by the model to LAMMPS.

Examples
""""""""

.. code-block:: LAMMPS

   pair_style metatensor exported-model.pt device cuda extensions /home/user/torch-extensions/
   pair_style metatensor soap-gap.pt check_consistency on
   pair_coeff * * 6 8 1

Description
"""""""""""

Pair style *metatensor* provides access to models following `metatensor's
atomistic models <https://docs.metatensor.org/latest/atomistic/index.html>`
interface; and enable using such models as interatomic potentials to drive a
LAMMPS simulation. The models can be fully defined and trained by the user using
Python code, or be existing pre-trained models. The interface can be used with
any type of machine learning model, as long as the implementation of the model
is compatible with TorchScript.

The only required argument for *pair_style metatensor* is the path to the model
file, which should be an exported metatensor model.

Optionally, users can define which torch *device* (e.g. cpu, cuda, cuda:0,
*etc.*) should be used to run the model. If this is not given, the code will run
on the best available device. If the model uses custom TorchScript operators
defined in a TorchScript extension, the shared library defining these extensions
will be searched in the *extensions* path, and loaded before trying to load the
model itself. Finally, *check_consistency* can be set to *on* or *off* to enable
(respectively disable) additional internal consistency checks in the data being
passed from LAMMPS to the model and back.

A single *pair_coeff* command should be used with the *metatensor* style,
specifying the mapping from LAMMPS types to the atomic types the model can
handle. The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
This is followed by a list of N arguments that specify the mapping of metatensor
atomic types to LAMMPS types, where N is the number of LAMMPS atom types.

.. See the :doc:`pair_coeff <pair_coeff>` page for alternate ways
.. to specify the path for the *model* and *extensions*.


Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support the :doc:`pair_modify <pair_modify>` shift,
table, and tail options.

This pair style does not write its information to :doc:`binary restart files
<restart>`, since it is stored in model files.  Thus, you need to re-specify the
pair_style and pair_coeff commands in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the :doc:`run_style
respa <run_style>` command.  It does not support the *inner*, *middle*, *outer*
keywords.

----------

Restrictions
""""""""""""

This pair style is part of the ML-METATENSOR package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>`
page for more info.


Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none
