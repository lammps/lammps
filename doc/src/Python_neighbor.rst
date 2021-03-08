Neighbor list access
====================

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
