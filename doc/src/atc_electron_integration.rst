.. index:: fix_modify AtC extrinsic electron_integration

fix_modify AtC extrinsic electron_integration command
=====================================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> extrinsic electron_integration <integration_type> [<num_subcycle_steps>]

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* extrinsic electron_integration = name of the AtC sub-command
* integration\_type = *explicit* or *implicit* or *steadydescriptor*
* num\_subcycle\_steps = number of subcycle steps for the electron time integration (optional)


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC extrinsic electron_integration implicit
   fix_modify AtC extrinsic electron_integration explicit 100

Description
"""""""""""

Switches between integration schemes for the electron temperature. The
number of subcycling steps used to integrate the electron temperature for
one LAMMPS timestep can be manually adjusted to capture fast electron
dynamics.

Restrictions
""""""""""""

For use only with the two\_temperature type of the AtC fix (see
:doc:`fix atc <fix_atc>` command)

Default
"""""""

*implicit* and *subcycle\_steps* = 1
