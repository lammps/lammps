Per-atom properties
===================

Similar to what is described in :doc:`Library_atoms`, the instances of
:py:class:`lammps <lammps.lammps>`, :py:class:`PyLammps <lammps.PyLammps>`, or
:py:class:`IPyLammps <lammps.IPyLammps>` can be used to extract atom quantities
and modify some of them. The main difference between the interfaces is how the information
is exposed.

While the :py:class:`lammps <lammps.lammps>` is just a thin layer that wraps C API calls,
:py:class:`PyLammps <lammps.PyLammps>` and :py:class:`IPyLammps <lammps.IPyLammps>` expose
information as objects and properties.

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
   rebuild and atom sort event as atoms are migrating between sub-domains.

.. tabs::

   .. tab:: lammps API

      .. code-block:: python

         from lammps import lammps

         lmp = lammps()
         lmp.file("in.sysinit")

         nlocal = lmp.extract_global("nlocal")
         x = lmp.extract_atom("x")

         for i in range(nlocal):
            print("(x,y,z) = (", x[i][0], x[i][1], x[i][2], ")")

         lmp.close()

      **Methods**:

      * :py:meth:`extract_atom() <lammps.lammps.extract_atom()>`: extract a per-atom quantity

      **Numpy Methods**:

      * :py:meth:`numpy.extract_atom() <lammps.numpy_wrapper.numpy_wrapper.extract_atom()>`: extract a per-atom quantity as numpy array

   .. tab:: PyLammps/IPyLammps API

      All atoms in the current simulation can be accessed by using the :py:attr:`PyLammps.atoms <lammps.PyLammps.atoms>` property.
      Each element of this list is a :py:class:`Atom <lammps.Atom>` or :py:class:`Atom2D <lammps.Atom2D>` object. The attributes of
      these objects provide access to their data (id, type, position, velocity, force, etc.):

      .. code-block:: Python

         # access first atom
         L.atoms[0].id
         L.atoms[0].type

         # access second atom
         L.atoms[1].position
         L.atoms[1].velocity
         L.atoms[1].force

      Some attributes can be changed:

      .. code-block:: Python

         # set position in 2D simulation
         L.atoms[0].position = (1.0, 0.0)

         # set position in 3D simulation
         L.atoms[0].position = (1.0, 0.0, 1.0)

