Call Python from a LAMMPS input script
======================================

LAMMPS has several commands which can be used to invoke Python
code directly from an input script:

* :doc:`python <python>`
* :doc:`variable python <variable>`
* :doc:`fix python/invoke <fix_python_invoke>`
* :doc:`pair_style python <pair_python>`

The :doc:`python <python>` command which can be used to define and
execute a Python function that you write the code for.  The Python
function can also be assigned to a LAMMPS python-style variable via
the :doc:`variable <variable>` command.  Each time the variable is
evaluated, either in the LAMMPS input script itself, or by another
LAMMPS command that uses the variable, this will trigger the Python
function to be invoked.

The Python code for the function can be included directly in the input
script or in an auxiliary file.  The function can have arguments which
are mapped to LAMMPS variables (also defined in the input script) and
it can return a value to a LAMMPS variable.  This is thus a mechanism
for your input script to pass information to a piece of Python code,
ask Python to execute the code, and return information to your input
script.

Note that a Python function can be arbitrarily complex.  It can import
other Python modules, instantiate Python classes, call other Python
functions, etc.  The Python code that you provide can contain more
code than the single function.  It can contain other functions or
Python classes, as well as global variables or other mechanisms for
storing state between calls from LAMMPS to the function.

The Python function you provide can consist of "pure" Python code that
only performs operations provided by standard Python.  However, the
Python function can also "call back" to LAMMPS through its
Python-wrapped library interface, in the manner described in the
:doc:`Python run <Python_run>` doc page.  This means it can issue LAMMPS
input script commands or query and set internal LAMMPS state.  As an
example, this can be useful in an input script to create a more
complex loop with branching logic, than can be created using the
simple looping and branching logic enabled by the :doc:`next <next>` and
:doc:`if <if>` commands.

See the :doc:`python <python>` doc page and the :doc:`variable <variable>`
doc page for its python-style variables for more info, including
examples of Python code you can write for both pure Python operations
and callbacks to LAMMPS.

The :doc:`fix python/invoke <fix_python_invoke>` command can execute
Python code at selected timesteps during a simulation run.

The :doc:`pair_style python <pair_python>` command allows you to define
pairwise potentials as python code which encodes a single pairwise
interaction.  This is useful for rapid development and debugging of a
new potential.

To use any of these commands, you only need to build LAMMPS with the
PYTHON package installed:

.. parsed-literal::

   make yes-python
   make machine

Note that this will link LAMMPS with the Python library on your
system, which typically requires several auxiliary system libraries to
also be linked.  The list of these libraries and the paths to find
them are specified in the lib/python/Makefile.lammps file.  You need
to insure that file contains the correct information for your version
of Python and your machine to successfully build LAMMPS.  See the
lib/python/README file for more info.

If you want to write Python code with callbacks to LAMMPS, then you
must also follow the steps summarized in the :doc:`Python run <Python_run>` doc page.  I.e. you must build LAMMPS as a shared
library and insure that Python can find the python/lammps.py file and
the shared library.
