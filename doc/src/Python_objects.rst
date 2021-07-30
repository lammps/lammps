Compute, fixes, variables
*************************

This section documents accessing or modifying data from objects like
computes, fixes, or variables in LAMMPS using the :py:mod:`lammps` module.

.. tabs::

   .. tab:: lammps API

      For :py:meth:`lammps.extract_compute() <lammps.lammps.extract_compute()>` and
      :py:meth:`lammps.extract_fix() <lammps.lammps.extract_fix()>`, the global, per-atom,
      or local data calculated by the compute or fix can be accessed. What is returned
      depends on whether the compute or fix calculates a scalar or vector or array.
      For a scalar, a single double value is returned.  If the compute or fix calculates
      a vector or array, a pointer to the internal LAMMPS data is returned, which you can
      use via normal Python subscripting.

      The one exception is that for a fix that calculates a
      global vector or array, a single double value from the vector or array
      is returned, indexed by I (vector) or I and J (array).  I,J are
      zero-based indices.
      See the :doc:`Howto output <Howto_output>` page for a discussion of
      global, per-atom, and local data, and of scalar, vector, and array
      data types.  See the doc pages for individual :doc:`computes <compute>`
      and :doc:`fixes <fix>` for a description of what they calculate and
      store.

      For :py:meth:`lammps.extract_variable() <lammps.lammps.extract_variable()>`,
      an :doc:`equal-style or atom-style variable <variable>` is evaluated and
      its result returned.

      For equal-style variables a single ``c_double`` value is returned and the
      group argument is ignored.  For atom-style variables, a vector of
      ``c_double`` is returned, one value per atom, which you can use via normal
      Python subscripting. The values will be zero for atoms not in the
      specified group.

      :py:meth:`lammps.numpy.extract_compute() <lammps.numpy_wrapper.numpy_wrapper.extract_compute()>`,
      :py:meth:`lammps.numpy.extract_fix() <lammps.numpy_wrapper.numpy_wrapper.extract_fix()>`, and
      :py:meth:`lammps.numpy.extract_variable() <lammps.numpy_wrapper.numpy_wrapper.extract_variable()>` are
      equivalent NumPy implementations that return NumPy arrays instead of ``ctypes`` pointers.

      The :py:meth:`lammps.set_variable() <lammps.lammps.set_variable()>` method sets an
      existing string-style variable to a new string value, so that subsequent LAMMPS
      commands can access the variable.

      **Methods**:

      * :py:meth:`lammps.extract_compute() <lammps.lammps.extract_compute()>`: extract value(s) from a compute
      * :py:meth:`lammps.extract_fix() <lammps.lammps.extract_fix()>`: extract value(s) from a fix
      * :py:meth:`lammps.extract_variable() <lammps.lammps.extract_variable()>`: extract value(s) from a variable
      * :py:meth:`lammps.set_variable() <lammps.lammps.set_variable()>`: set existing named string-style variable to value

      **NumPy Methods**:

      * :py:meth:`lammps.numpy.extract_compute() <lammps.numpy_wrapper.numpy_wrapper.extract_compute()>`: extract value(s) from a compute, return arrays as numpy arrays
      * :py:meth:`lammps.numpy.extract_fix() <lammps.numpy_wrapper.numpy_wrapper.extract_fix()>`: extract value(s) from a fix, return arrays as numpy arrays
      * :py:meth:`lammps.numpy.extract_variable() <lammps.numpy_wrapper.numpy_wrapper.extract_variable()>`: extract value(s) from a variable, return arrays as numpy arrays


   .. tab:: PyLammps/IPyLammps API

      PyLammps and IPyLammps classes currently do not add any additional ways of
      retrieving information out of computes and fixes. This information can still be accessed by using the lammps API:

      .. code-block:: python

         L.lmp.extract_compute(...)
         L.lmp.extract_fix(...)
         # OR
         L.lmp.numpy.extract_compute(...)
         L.lmp.numpy.extract_fix(...)

      LAMMPS variables can be both defined and accessed via the :py:class:`PyLammps <lammps.PyLammps>` interface.

      To define a variable you can use the :doc:`variable <variable>` command:

      .. code-block:: Python

         L.variable("a index 2")

      A dictionary of all variables is returned by the :py:attr:`PyLammps.variables <lammps.PyLammps.variables>` property:

      you can access an individual variable by retrieving a variable object from the
      ``L.variables`` dictionary by name

      .. code-block:: Python

         a = L.variables['a']

      The variable value can then be easily read and written by accessing the value
      property of this object.

      .. code-block:: Python

         print(a.value)
         a.value = 4
