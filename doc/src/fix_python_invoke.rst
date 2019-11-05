.. index:: fix python/invoke

fix python/invoke command
=========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID python/invoke N callback function_name

* ID, group-ID are ignored by this fix
* python/invoke = style name of this fix command
* N = execute every N steps
* callback = *post\_force* or *end\_of\_step*
  
  .. parsed-literal::
  
       *post_force* = callback after force computations on atoms every N time steps
       *end_of_step* = callback after every N time steps



Examples
""""""""


.. parsed-literal::

   python post_force_callback here """
   from lammps import lammps

   def post_force_callback(lammps_ptr, vflag):
       lmp = lammps(ptr=lammps_ptr)
       # access LAMMPS state using Python interface
   """

   python end_of_step_callback here """
   def end_of_step_callback(lammps_ptr):
       lmp = lammps(ptr=lammps_ptr)
       # access LAMMPS state using Python interface
   """

   fix pf  all python/invoke 50 post_force post_force_callback
   fix eos all python/invoke 50 end_of_step end_of_step_callback

Description
"""""""""""

This fix allows you to call a Python function during a simulation run.
The callback is either executed after forces have been applied to atoms
or at the end of every N time steps.

Callback functions must be declared in the global scope of the
active Python interpreter. This can either be done by defining it
inline using the python command or by importing functions from other
Python modules. If LAMMPS is driven using the library interface from
Python, functions defined in the driving Python interpreter can also
be executed.

Each callback is given a pointer object as first argument. This can be
used to initialize an instance of the lammps Python interface, which
gives access to the LAMMPS state from Python.

.. warning::

   While you can access the state of LAMMPS via library functions
   from these callbacks, trying to execute input script commands will in the best
   case not work or in the worst case result in undefined behavior.

Restrictions
""""""""""""


This fix is part of the PYTHON package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Building LAMMPS with the PYTHON package will link LAMMPS with the
Python library on your system.  Settings to enable this are in the
lib/python/Makefile.lammps file.  See the lib/python/README file for
information on those settings.

Related commands
""""""""""""""""

:doc:`python command <python>`


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
