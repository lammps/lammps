.. index:: compute_modify

compute_modify command
======================

Syntax
""""""

.. code-block:: LAMMPS

   compute_modify compute-ID keyword value ...

* compute-ID = ID of the compute to modify
* one or more keyword/value pairs may be listed
* keyword = *extra/dof* or *dynamic/dof*

  .. parsed-literal::

       *extra/dof* value = N
         N = # of extra degrees of freedom to subtract
       *dynamic/dof* value = *yes* or *no*
         yes/no = do or do not re-compute the number of degrees of freedom (DOF) contributing to the temperature

Examples
""""""""

.. code-block:: LAMMPS

   compute_modify myTemp extra/dof 0
   compute_modify newtemp dynamic/dof yes extra/dof 600

Description
"""""""""""

Modify one or more parameters of a previously defined compute.  Not
all compute styles support all parameters.

The *extra/dof* keyword refers to how many degrees of freedom are
subtracted (typically from :math:`3N`) as a normalizing factor in a
temperature computation.  Only computes that compute a temperature use
this option.  The default is 2 or 3 for :doc:`2d or 3d systems
<dimension>` which is a correction factor for an ensemble of velocities
with zero total linear momentum. For compute temp/partial, if one or
more velocity components are excluded, the value used for *extra/dof* is
scaled accordingly. You can use a negative number for the *extra/dof*
parameter if you need to add degrees-of-freedom.  See the :doc:`compute
temp/asphere <compute_temp_asphere>` command for an example.

The *dynamic/dof* keyword determines whether the number
of atoms :math:`N` in the compute group and their associated degrees
of freedom (DOF) are re-computed each time a temperature is computed.
Only compute styles that calculate a temperature use this option.  By
default, :math:`N` and their DOF are assumed to be constant.  If you
are adding atoms or molecules to the system (see the :doc:`fix pour
<fix_pour>`, :doc:`fix deposit <fix_deposit>`, and :doc:`fix gcmc
<fix_gcmc>` commands) or expect atoms or molecules to be lost
(e.g. due to exiting the simulation box or via :doc:`fix evaporate
<fix_evaporate>`), then this option should be used to ensure the
temperature is correctly normalized.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute <compute>`

Default
"""""""

The option defaults are extra/dof = 2 or 3 for 2d or 3d systems,
respectively, and dynamic/dof = *no*.
