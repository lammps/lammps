.. index:: pair_style mliap

pair_style mliap command
========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style mliap ... keyword values ...

* two keyword/value pairs must be appended
* keyword = *model* or *descriptor*

  .. parsed-literal::

       *model* values = style filename
         style = *linear* or *quadratic*
         filename = name of file containing model definitions
       *descriptor* values = style filename
         style = *sna*
         filename = name of file containing descriptor definitions

Examples
""""""""

.. code-block:: LAMMPS

   pair_style mliap model linear InP.mliap.model descriptor sna InP.mliap.descriptor
   pair_style mliap model quadratic W.mliap.model descriptor sna W.mliap.descriptor
   pair_coeff * * In P

Description
"""""""""""

Pair style *mliap* provides a general interface to families of 
machine-learning interatomic potentials. It allows separate 
definitions of the interatomic potential functional form (*model*)
and the geometric quantities that characterize the atomic positions
(*descriptor*). By defining *model* and *descriptor* separately,
it is possible to use many different models with a given descriptor,
or many different descriptors with a given model. Currently, the pair_style
supports just two models, *linear* and *quadratic*,
and one descriptor, *sna*, the SNAP descriptor used by :doc:`pair_style snap <pair_snap>`, including the linear, quadratic,
and chem variants. Work is currently underway to extend
the interface to handle neural network energy models,
and it is also straightforward to add new descriptor styles.
In order to train a model, it is useful to know the gradient or derivative
of energy, force, and stress w.r.t. model parameters. This information
can be accessed using the related :doc:`compute mliap <compute_mliap>` command.

The pair_style *mliap* command must be followed by two keywords
*model* and *descriptor* in either order. A single
*pair_coeff* command is also required. The first 2 arguments
must be \* \* so as to span all LAMMPS atom types.
This is followed by a list of N arguments
that specify the mapping of MLIAP
element names to LAMMPS atom types,
where N is the number of LAMMPS atom types.

The *model* keyword is followed by a model style, currently limited to
either *linear* or *quadratic*. In both cases,
this is followed by a single argument specifying the model filename containing the 
parameters for a set of elements. 
The model filename usually ends in the *.mliap.model* extension.
It may contain parameters for many elements. The only requirement is that it
contain at least those element names appearing in the
*pair_coeff* command.

The top of the model file can contain any number of blank and comment lines (start with #),
but follows a strict format after that. The first non-blank non-comment
line must contain two integers:

* nelems  = Number of elements
* nparams = Number of parameters

This is followed by one block for each of the *nelem* elements.
Each block consists of *nparams* parameters, one per line.
Note that this format is similar, but not identical to that used
for the :doc:`pair_style snap <pair_snap>` coefficient file.
Specifically, the line containing the element weight and radius is omitted,
since these are handled by the *descriptor*.

The *descriptor* keyword is followed by a descriptor style, and additional arguments.
Currently the only descriptor style is *sna*, indicating the bispectrum component
descriptors used by the Spectral Neighbor Analysis Potential (SNAP) potentials of
:doc:`pair_style snap <pair_snap>`.
The \'p\' in SNAP is dropped, because keywords that match pair_styles are silently stripped
out by the LAMMPS command parser. A single additional argument specifies the descriptor filename
containing the parameters and setting used by the SNAP descriptor.
The descriptor filename usually ends in the *.mliap.descriptor* extension.

The SNAP descriptor file closely follows the format of the
:doc:`pair_style snap <pair_snap>` parameter file.
The file can contain blank and comment lines (start
with #) anywhere. Each non-blank non-comment line must contain one
keyword/value pair. The required keywords are *rcutfac* and
*twojmax*\ . There are many optional keywords that are described
on the :doc:`pair_style snap <pair_snap>` doc page.
In addition, the SNAP descriptor file must contain
the *nelems*, *elems*, *radelems*, and *welems* keywords.
The *nelems* keyword specifies the number of elements
provided in the other three keywords.
The *elems* keyword is followed by a list of *nelems*
element names that must include the element
names appearing in the *pair_coeff* command,
but can contain other names too.
Similarly, the *radelems* and *welems* keywords are
followed by lists of *nelems* numbers giving the element radius
and element weight of each element. Obviously, the order
in which the elements are listed must be consistent for all
three keywords.

See the :doc:`pair_coeff <pair_coeff>` doc page for alternate ways
to specify the path for these *model* and *descriptor* files.

**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS with
user-specifiable parameters as described above.  You never need to
specify a pair_coeff command with I != J arguments for this style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.

----------

Restrictions
""""""""""""

This style is part of the MLIAP package.  It is only enabled if LAMMPS
was built with that package. In addition, building LAMMPS with the MLIAP package
requires building LAMMPS with the SNAP package.
See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_style snap  <pair_snap>`, :doc:`compute mliap <compute_mliap>`

**Default:** none
