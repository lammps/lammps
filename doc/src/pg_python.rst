The ``lammps`` Python module
****************************

.. py:module:: lammps

The LAMMPS Python interface is implemented as a module called
:py:mod:`lammps` in the ``lammps.py`` file in the ``python`` folder of
the LAMMPS source code distribution.  After compilation of LAMMPS, the
module can be installed into a Python system folder or a user folder
with ``make install-python``.  Components of the module can then loaded
into a Python session with the ``import`` command.

There are multiple Python interface classes in the :py:mod:`lammps` module:

- the :py:class:`lammps <lammps.lammps>` class. This is a wrapper around
  the C-library interface and its member functions try to replicate the
  :doc:`C-library API <pg_library>` closely.  This is the most
  feature-complete Python API.
- the :py:class:`PyLammps <lammps.PyLammps>` class. This is a more high-level
  and more Python style class implemented on top of the
  :py:class:`lammps <lammps.lammps>` class.
- the :py:class:`IPyLammps <lammps.IPyLammps>` class is derived from
  :py:class:`PyLammps <lammps.PyLammps>` and adds embedded graphics
  features to conveniently include LAMMPS into `Jupyter
  <https://jupyter.org/>`_ notebooks.

.. _mpi4py_url: https://mpi4py.readthedocs.io


