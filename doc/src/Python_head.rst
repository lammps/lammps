Use Python with LAMMPS
**********************

These pages describe various ways that LAMMPS and Python can be used
together.

.. toctree::
   :maxdepth: 1

   Python_overview
   Python_install
   Python_run
   Python_module
   Python_ext
   Python_call
   Python_formats
   Python_examples
   Python_error
   Python_trouble

If you are not familiar with `Python <https://www.python.org>`_, it is a
powerful scripting and programming language which can do almost
everything that compiled languages like C, C++, or Fortran can do in
fewer lines of code. It also comes with a large collection of add-on
modules for many purposes (either bundled or easily installed from
Python code repositories).  The major drawback is slower execution speed
of the script code compared to compiled programming languages.  But when
the script code is interfaced to optimized compiled code, performance
can be on par with a standalone executable, for as long as the scripting
is restricted to high-level operations.  Thus Python is also convenient
to use as a "glue" language to "drive" a program through its library
interface, or to hook multiple pieces of software together, such as a
simulation code and a visualization tool, or to run a coupled
multi-scale or multi-physics model.

See the :doc:`Howto_couple` page for more ideas about coupling LAMMPS
to other codes.  See the :doc:`Library` page for a description of the
LAMMPS library interfaces.  That interface is exposed to Python either
when calling LAMMPS from Python or when calling Python from a LAMMPS
input script and then calling back to LAMMPS from Python code.  The
C-library interface is designed to be easy to add functionality to,
thus the Python interface to LAMMPS is easy to extend as well.

If you create interesting Python scripts that run LAMMPS or
interesting Python functions that can be called from a LAMMPS input
script, that you think would be generally useful, please post them as
a pull request to our `GitHub site <https://github.com/lammps/lammps>`_,
and they can be added to the LAMMPS distribution or web page.
