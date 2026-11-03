Neighbor list access
====================

Access to neighbor lists is handled through a couple of wrapper classes
that allows one to treat it like either a python list or a NumPy array.  The
access procedure is similar to that of the C-library interface: use one
of the "find" functions to look up the index of the neighbor list in the
global table of neighbor lists and then get access to the neighbor list
data. The code sample below demonstrates reading the neighbor list data
using the NumPy access method.

.. note::

   The neighbor arrays are part of the neighbor list itself, so the
   neighbor indices are stored the same way LAMMPS uses them internally.
   For systems with :doc:`special neighbors <special_bonds>` (pairs of
   atoms connected by bonds, angles, or dihedrals), the two highest bits
   of a neighbor index may be set to record whether that pair is a 1-2,
   1-3, or 1-4 neighbor.  Such indices are very large or negative numbers
   and must have those bits removed before they can be used as local atom
   indices, e.g. with ``j = nlist[jj] & 0x1FFFFFFF`` (or ``nlist &
   0x1FFFFFFF`` for the whole NumPy array).  The value of ``(nlist[jj]
   >> 30) & 3`` is 1, 2, or 3 for 1-2, 1-3, or 1-4 neighbors,
   respectively, and 0 for all other neighbors.

   No bits are set for special neighbors with weighting factors of 1.0
   for both Lennard-Jones and Coulomb interactions.  Special neighbors
   with weighting factors of 0.0 for both are removed from the neighbor
   list instead, except for the cases listed in the :doc:`special_bonds
   <special_bonds>` command documentation (e.g. when a KSpace style is
   used), where the bits are set for all special neighbors.

.. code-block:: python

   from lammps import lammps
   import numpy as np

   lmp = lammps()
   lmp.commands_string("""
   region box block -2 2 -2 2 -2 2
   lattice fcc 1.0
   create_box 1 box
   create_atoms 1 box
   mass 1 1.0
   pair_style lj/cut 2.5
   pair_coeff 1 1 1.0 1.0
   run 0 post no""")

   # look up the neighbor list
   nlidx = lmp.find_pair_neighlist('lj/cut')
   nl = lmp.numpy.get_neighlist(nlidx)
   tags = lmp.extract_atom('id')
   print("half neighbor list with {} entries".format(nl.size))
   # print neighbor list contents
   for i in range(0,nl.size):
       idx, nlist  = nl.get(i)
       print("\natom {} with ID {} has {} neighbors:".format(idx,tags[idx],nlist.size))
       if nlist.size > 0:
           for n in np.nditer(nlist):
               print("  atom {} with ID {}".format(n,tags[n]))

Another example for extracting a full neighbor list without evaluating a
potential is shown below.

.. code-block:: python

    from lammps import lammps
    import numpy as np

    lmp = lammps()
    lmp.commands_string("""
    newton off
    region box block -2 2 -2 2 -2 2
    lattice fcc 1.0
    create_box 1 box
    create_atoms 1 box
    mass 1 1.0
    pair_style zero 1.0 full
    pair_coeff * *
    run 0 post no""")

    # look up the neighbor list
    nlidx = lmp.find_pair_neighlist('zero')
    nl = lmp.numpy.get_neighlist(nlidx)
    tags = lmp.extract_atom('id')
    print("full neighbor list with {} entries".format(nl.size))
    # print neighbor list contents
    for i in range(0,nl.size):
        idx, nlist  = nl.get(i)
        print("\natom {} with ID {} has {} neighbors:".format(idx,tags[idx],nlist.size))
        if nlist.size > 0:
            for n in np.nditer(nlist):
                pass
                print("  atom {} with ID {}".format(n,tags[n]))

**Methods:**

* :py:meth:`lammps.get_neighlist() <lammps.lammps.get_neighlist()>`: Get neighbor list for given index
* :py:meth:`lammps.get_neighlist_size() <lammps.lammps.get_neighlist_size()>`: Get number of elements in neighbor list
* :py:meth:`lammps.get_neighlist_element_neighbors() <lammps.lammps.get_neighlist_element_neighbors()>`: Get element in neighbor list and its neighbors

* :py:meth:`lammps.find_pair_neighlist() <lammps.lammps.find_pair_neighlist()>`: Find neighbor list of pair style
* :py:meth:`lammps.find_fix_neighlist() <lammps.lammps.find_pair_neighlist()>`: Find neighbor list of fix style
* :py:meth:`lammps.find_compute_neighlist() <lammps.lammps.find_pair_neighlist()>`: Find neighbor list of compute style


**NumPy Methods:**

* :py:meth:`lammps.numpy.get_neighlist() <lammps.numpy_wrapper.numpy_wrapper.get_neighlist()>`: Get neighbor list for given index, which uses NumPy arrays for its element neighbor arrays
* :py:meth:`lammps.numpy.get_neighlist_element_neighbors() <lammps.numpy_wrapper.numpy_wrapper.get_neighlist_element_neighbors()>`: Get element in neighbor list and its neighbors (as a numpy array)
