.. index:: fix python/move

fix python/move command
=======================

Syntax
""""""

.. parsed-literal::

   fix python/move pymodule.CLASS

pymodule.CLASS = use class **CLASS** in module/file **pymodule** to compute how to move atoms

Examples
""""""""

.. code-block:: LAMMPS

   fix  1 all python/move py_nve.NVE
   fix  1 all python/move py_nve.NVE_OPT

Description
"""""""""""

The *python/move* fix style provides a way to define ways how particles
are moved during an MD run from python script code, that is loaded from
a file into LAMMPS and executed at the various steps where other fixes
can be executed. This python script must contain specific python class
definitions.

This allows to implement complex position updates and also modified
time integration methods. Due to python being an interpreted language,
however, the performance of this fix can be moderately to significantly
slower than the corresponding C++ code. For specific cases, this
performance penalty can be limited through effective use of NumPy.

----------

The python module file has to start with the following code:

.. code-block:: python

   from __future__ import print_function
   import lammps
   import ctypes
   import traceback
   import numpy as np
   #
   class LAMMPSFix(object):
       def __init__(self, ptr, group_name="all"):
           self.lmp = lammps.lammps(ptr=ptr)
           self.group_name = group_name
   #
   class LAMMPSFixMove(LAMMPSFix):
       def __init__(self, ptr, group_name="all"):
           super(LAMMPSFixMove, self).__init__(ptr, group_name)
   #
       def init(self):
           pass
   #
       def initial_integrate(self, vflag):
           pass
   #
       def final_integrate(self):
           pass
   #
       def initial_integrate_respa(self, vflag, ilevel, iloop):
           pass
   #
       def final_integrate_respa(self, ilevel, iloop):
           pass
   #
       def reset_dt(self):
           pass

Any classes implementing new atom motion functionality have to be
derived from the **LAMMPSFixMove** class, overriding the available
methods as needed.

Examples for how to do this are in the *examples/python* folder.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This pair style is part of the PYTHON package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`, :doc:`fix python/invoke <fix_python_invoke>`

Default
"""""""

none
