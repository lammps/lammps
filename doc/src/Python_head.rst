Use Python with LAMMPS
**********************

These doc pages describe various ways that LAMMPS and Python can be
used together.


.. toctree::
   :maxdepth: 1

   Python_overview
   Python_run
   Python_shlib
   Python_install
   Python_mpi
   Python_test
   Python_library
   Python_pylammps
   Python_examples
   Python_call

If you're not familiar with `Python <http://www.python.org>`_, it's a
powerful scripting and programming language which can do most
everything that lower-level languages like C or C++ can do in fewer
lines of code.  The only drawback is slower execution speed.  Python
is also easy to use as a "glue" language to drive a program through
its library interface, or to hook multiple pieces of software
together, such as a simulation code plus a visualization tool, or to
run a coupled multiscale or multiphysics model.

See the :doc:`Howto\_couple <Howto_couple>` doc page for more ideas about
coupling LAMMPS to other codes.  See the :doc:`Howto library <Howto_library>` doc page for a description of the LAMMPS
library interface provided in src/library.h and src/library.h.  That
interface is exposed to Python either when calling LAMMPS from Python
or when calling Python from a LAMMPS input script and then calling
back to LAMMPS from Python code.  The library interface is designed to
be easy to add functionality to.  Thus the Python interface to LAMMPS
is also easy to extend as well.

If you create interesting Python scripts that run LAMMPS or
interesting Python functions that can be called from a LAMMPS input
script, that you think would be generally useful, please post them as
a pull request to our `GitHub site <https://github.com/lammps/lammps>`_,
and they can be added to the LAMMPS distribution or webpage.
