.. index:: pair_style python

pair_style python command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style python cutoff

cutoff = global cutoff for interactions in python potential classes

Examples
""""""""

.. code-block:: LAMMPS

   pair_style python 2.5
   pair_coeff * * py_pot.LJCutMelt lj

   pair_style hybrid/overlay coul/long 12.0 python 12.0
   pair_coeff * * coul/long
   pair_coeff * * python py_pot.LJCutSPCE OW NULL

Description
"""""""""""

The *python* pair style provides a way to define pairwise additive
potential functions as python script code that is loaded into LAMMPS
from a python file which must contain specific python class definitions.
This allows to rapidly evaluate different potential functions without
having to modify and re-compile LAMMPS. Due to python being an
interpreted language, however, the performance of this pair style is
going to be significantly slower (often between 20x and 100x) than
corresponding compiled code. This penalty can be significantly reduced
through generating tabulations from the python code through the
:doc:`pair_write <pair_write>` command, which is supported by this style.

Only a single :doc:`pair_coeff <pair_coeff>` command is used with the
*python* pair style which specifies a python class inside a python module
or a file that LAMMPS will look up in the current directory, a folder
pointed to by the ``LAMMPS_POTENTIALS`` environment variable or somewhere
in your python path.  A single python module can hold multiple python pair
class definitions.  The class definitions itself have to follow specific
rules that are explained below.

Atom types in the python class are specified through symbolic
constants, typically strings. These are mapped to LAMMPS atom types by
specifying N additional arguments after the class name in the
pair_coeff command, where N must be the number of currently defined
atom types:

As an example, imagine a file *py_pot.py* has a python potential class
names *LJCutMelt* with parameters and potential functions for a two
Lennard-Jones atom types labeled as 'LJ1' and 'LJ2'. In your LAMMPS
input and you would have defined 3 atom types, out of which the first
two are supposed to be using the 'LJ1' parameters and the third the
'LJ2' parameters, then you would use the following pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * py_pot.LJCutMelt LJ1 LJ1 LJ2

The first two arguments **must** be \* \* so as to span all LAMMPS atom
types.  The first two LJ1 arguments map LAMMPS atom types 1 and 2 to
the LJ1 atom type in the LJCutMelt class of the py_pot.py file.  The
final LJ2 argument maps LAMMPS atom type 3 to the LJ2 atom type the
python file.  If a mapping value is specified as NULL, the mapping is
not performed, any pair interaction with this atom type will be
skipped. This can be used when a *python* potential is used as part of
the *hybrid* or *hybrid/overlay* pair style. The NULL values are then
placeholders for atom types that will be used with other potentials.

----------

The python potential file has to start with the following code:

.. code-block:: python

   from __future__ import print_function

   class LAMMPSPairPotential(object):
       def __init__(self):
           self.pmap=dict()
           self.units='lj'
       def map_coeff(self,name,ltype):
           self.pmap[ltype]=name
       def check_units(self,units):
           if (units != self.units):
              raise Exception("Conflicting units: %s vs. %s" % (self.units,units))

Any classes with definitions of specific potentials have to be derived
from this class and should be initialize in a similar fashion to the
example given below.

.. note::

   The class constructor has to set up a data structure containing
   the potential parameters supported by this class.  It should also
   define a variable *self.units* containing a string matching one of the
   options of LAMMPS' :doc:`units <units>` command, which is used to
   verify, that the potential definition in the python class and in the
   LAMMPS input match.

Here is an example for a single type Lennard-Jones potential class
*LJCutMelt* in reduced units, which defines an atom type *lj* for
which the parameters epsilon and sigma are both 1.0:

.. code-block:: python

   class LJCutMelt(LAMMPSPairPotential):
       def __init__(self):
           super(LJCutMelt,self).__init__()
           # set coeffs: 48*eps*sig**12, 24*eps*sig**6,
           #              4*eps*sig**12,  4*eps*sig**6
           self.units = 'lj'
           self.coeff = {'lj'  : {'lj'  : (48.0,24.0,4.0,4.0)}}

The class also has to provide two methods for the computation of the
potential energy and forces, which have be named *compute_force*,
and *compute_energy*, which both take 3 numerical arguments:

* rsq   = the square of the distance between a pair of atoms (float)
* itype = the (numerical) type of the first atom
* jtype = the (numerical) type of the second atom

This functions need to compute the force and the energy, respectively,
and use the result as return value. The functions need to use the
*pmap* dictionary to convert the LAMMPS atom type number to the symbolic
value of the internal potential parameter data structure. Following
the *LJCutMelt* example, here are the two functions:

.. code-block:: python

      def compute_force(self,rsq,itype,jtype):
           coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
           r2inv  = 1.0/rsq
           r6inv  = r2inv*r2inv*r2inv
           lj1 = coeff[0]
           lj2 = coeff[1]
           return (r6inv * (lj1*r6inv - lj2))*r2inv

       def compute_energy(self,rsq,itype,jtype):
           coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
           r2inv  = 1.0/rsq
           r6inv  = r2inv*r2inv*r2inv
           lj3 = coeff[2]
           lj4 = coeff[3]
           return (r6inv * (lj3*r6inv - lj4))

.. note::

   for consistency with the C++ pair styles in LAMMPS, the
   *compute_force* function follows the conventions of the Pair::single()
   methods and does not return the full force, but the force scaled by
   the distance between the two atoms, so this value only needs to be
   multiplied by delta x, delta y, and delta z to conveniently obtain the
   three components of the force vector between these two atoms.

----------

.. note::

   The evaluation of scripted python code will slow down the
   computation pair-wise interactions quite significantly. However, this
   can be largely worked around through using the python pair style not
   for the actual simulation, but to generate tabulated potentials on the
   fly using the :doc:`pair_write <pair_write>` command. Please see below
   for an example LAMMPS input of how to build a table file:

.. code-block:: LAMMPS

   pair_style python 2.5
   pair_coeff * * py_pot.LJCutMelt lj
   shell rm -f melt.table
   pair_write  1 1 2000 rsq 0.01 2.5 lj1_lj2.table lj

Note that it is strongly recommended to try to **delete** the potential
table file before generating it. Since the *pair_write* command will
always **append** to a table file, while pair style table will use the
**first match**\ . Thus when changing the potential function in the python
class, the table pair style will still read the old variant unless the
table file is first deleted.

After switching the pair style to *table*, the potential tables need
to be assigned to the LAMMPS atom types like this:

.. code-block:: LAMMPS

   pair_style      table linear 2000
   pair_coeff      1  1  melt.table lj

This can also be done for more complex systems.  Please see the
*examples/python* folders for a few more examples.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Mixing of potential parameters has to be handled inside the provided
python module. The python pair style simply assumes that force and
energy computation can be correctly performed for all pairs of atom
types as they are mapped to the atom type labels inside the python
potential class.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the PYTHON package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_write <pair_write>`,
:doc:`pair style table <pair_table>`

Default
"""""""

none
