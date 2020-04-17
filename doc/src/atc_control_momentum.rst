.. index:: fix_modify AtC control momentum

fix_modify AtC control momentum command
=======================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> control <physics_type> <solution_parameter> <value>
   fix_modify AtC control momentum none
   fix_modify AtC control momentum rescale <frequency>
   fix_modify AtC control momentum glc_displacement
   fix_modify AtC control momentum glc_velocity
   fix_modify AtC control momentum hoover
   fix_modify AtC control momentum flux [faceset face_set_id, interpolate]
   
* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* control = name of the AtC sub-command
* physics_type = *thermal* or *momentum*
* solution_parameter = *max_iterations* or *tolerance*
* value = solution_parameter value
* momentum option = *none* or *rescale* or *glc_displacement* or *glc_velocity* *hoover* or *flux*
* frequency = time step frequency for applying displacement and velocity rescaling
* faceset_id = id of boundary face set (optional, only for *faceset*)


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC control momentum none
   fix_modify AtC control momentum flux faceset bndy_faces
   fix_modify AtC control momentum glc_velocity

Description
"""""""""""

The general version of *control* sets the numerical parameters for the
matrix solvers used in the specified control algorithm.  Many solution
approaches require iterative solvers, and these methods enable users to
provide the maximum number of iterations and the relative tolerance.

The *control momentum* version sets the momentum exchange mechanism from
the finite elements to the atoms, managed through a control algorithm.
*rescale* computes a scale factor for each atom to match the finite
element temperature.  *hoover* is a Gaussian least-constraint isokinetic
thermostat enforces that the nodal restricted atomic temperature matches
the finite element temperature.  *flux* is a similar mode, but rather
adds energy to the atoms based on conservation of energy.

*correction_max_iterations* sets the maximum number of iterations to
compute the second order in time correction term for lambda with the
fractional step method. The method uses the same tolerance as the
controller's matrix solver.

Restrictions
""""""""""""

Only for be used with the specific controllers *thermal* or *momentum*.
They are ignored if a lumped solution is requested.

*control momentum* is only for be used with specific transfers: elastic
*rescale* not valid with time filtering activated

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC control thermal <atc_control_thermal>`

Default
"""""""

- *max_iterations* is the number of rows in the matrix.
- *tolerance* is 1.0e-10.
