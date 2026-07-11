Per-atom properties
===================

Similar to what is described in :doc:`Library_atoms`, the instances of
:py:class:`lammps <lammps.lammps>` can be used to extract atom quantities
and modify some of them.

In some cases the data returned is a direct reference to the original data
inside LAMMPS cast to ``ctypes`` pointers. Where possible, the wrappers will
determine the ``ctypes`` data type and cast pointers accordingly. If
``numpy`` is installed arrays can also be extracted as numpy arrays, which
will access the C arrays directly and have the correct dimensions to protect
against invalid accesses.

.. warning::

   When accessing per-atom data,
   please note that this data is the per-processor local data and indexed
   accordingly. These arrays can change sizes and order at every neighbor list
   rebuild and atom sort event as atoms are migrating between subdomains.

.. code-block:: python

   from lammps import lammps

   lmp = lammps()
   lmp.file("in.sysinit")


   # Read/Write access via ctypes
   nlocal = lmp.extract_global("nlocal")
   x = lmp.extract_atom("x")

   for i in range(nlocal):
      print("(x,y,z) = (", x[i][0], x[i][1], x[i][2], ")")

   # Read/Write access via NumPy arrays
   atom_id = L.numpy.extract_atom("id")
   atom_type = L.numpy.extract_atom("type")
   x = L.numpy.extract_atom("x")
   v = L.numpy.extract_atom("v")
   f = L.numpy.extract_atom("f")

   # set position in 2D simulation
   x[0] = (1.0, 0.0)

   # set position in 3D simulation
   x[0] = (1.0, 0.0, 1.)

   lmp.close()


**Methods**:

* :py:meth:`extract_atom() <lammps.lammps.extract_atom()>`: extract a per-atom quantity

**Numpy Methods**:

* :py:meth:`numpy.extract_atom() <lammps.numpy_wrapper.numpy_wrapper.extract_atom()>`: extract a per-atom quantity as numpy array
