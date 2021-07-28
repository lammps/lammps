Neighbor list access
====================

Access to neighbor lists is handled through a couple of wrapper classes
that allows to treat it like either a python list or a NumPy array.  The
access procedure is similar to that of the C-library interface: use one
of the "find" functions to look up the index of the neighbor list in the
global table of neighbor lists and then get access to the neighbor list
data. The code sample below demonstrates reading the neighbor list data
using the NumPy access method.

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

**Methods:**

* :py:meth:`lammps.get_neighlist() <lammps.lammps.get_neighlist()>`: Get neighbor list for given index
* :py:meth:`lammps.get_neighlist_size()`: Get number of elements in neighbor list
* :py:meth:`lammps.get_neighlist_element_neighbors()`: Get element in neighbor list and its neighbors

* :py:meth:`lammps.find_pair_neighlist() <lammps.lammps.find_pair_neighlist()>`: Find neighbor list of pair style
* :py:meth:`lammps.find_fix_neighlist() <lammps.lammps.find_pair_neighlist()>`: Find neighbor list of pair style
* :py:meth:`lammps.find_compute_neighlist() <lammps.lammps.find_pair_neighlist()>`: Find neighbor list of pair style


**NumPy Methods:**

* :py:meth:`lammps.numpy.get_neighlist() <lammps.numpy_wrapper.numpy_wrapper.get_neighlist()>`: Get neighbor list for given index, which uses NumPy arrays for its element neighbor arrays
* :py:meth:`lammps.numpy.get_neighlist_element_neighbors() <lammps.numpy_wrapper.numpy_wrapper.get_neighlist_element_neighbors()>`: Get element in neighbor list and its neighbors (as numpy array)
