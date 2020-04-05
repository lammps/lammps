LAMMPS Command Classes
**********************

In LAMMPS most commands are implemented as classes
derived from the :cpp:class:`LAMMPS_NS::Pointers`
class which have to contain a ``command()`` method.
Creating the class instance corresponds to an
input file command and its arguments are passed
to the ``command()`` which will then execute the
command.

----------

.. doxygenclass:: LAMMPS_NS::Run
   :project: progguide
   :members:
   :protected-members:

----------

.. doxygenclass:: LAMMPS_NS::WriteData
   :project: progguide
   :members:
   :protected-members:

