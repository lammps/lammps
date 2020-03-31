LAMMPS Command Classes
**********************

In LAMMPS most commands are implemented as classes
derived from the :ref:`Pointers <lammps_ns_pointers>`
class which have to contain a ``command()`` method.
Creating the class instance corresponds to an
input file command and its arguments are passed
to the ``command()`` which will then execute the
command.

----------

.. _lammps_ns_run:
.. doxygenclass:: LAMMPS_NS::Run
   :project: progguide
   :members:
   :protected-members:
   :no-link:

----------

.. _lammps_ns_write_data:
.. doxygenclass:: LAMMPS_NS::WriteData
   :project: progguide
   :members:
   :protected-members:
   :no-link:

