.. index:: fix_modify AtC control thermal

fix_modify AtC control thermal command
======================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> control <physics_type> <solution_parameter> <value>
   fix_modify <AtC fixID> control thermal <control_type> <optional_args>
   fix_modify <AtC fixID> control thermal rescale <frequency>
   fix_modify <AtC fixID> control thermal flux <boundary_integration_type> <faceset_id>
   fix_modify <AtC fixID> control thermal correction_max_iterations <max_iterations>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* control = name of the AtC sub-command
* physics_type = *thermal* or *momentum*
* solution_parameter = *max_iterations* or *tolerance*
* value = solution_parameter value
* thermal control_type = *none* or *rescale* or *hoover* or *flux*
* frequency = time step frequency for applying velocity rescaling
* boundary_integration_type = *faceset* or *interpolate* (optional)
* faceset_id = id of boundary face set (optional, only for *faceset*)
* correction_max_iterations = maximum number of iterations that will be used by iterative matrix solvers for *thermal* physics type

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC control thermal none
   fix_modify AtC control thermal rescale 10
   fix_modify AtC control thermal hoover
   fix_modify AtC control thermal flux
   fix_modify AtC control thermal flux faceset bndy_faces
   fix_modify AtC control thermal correction_max_iterations 10

Description
"""""""""""

The general version of *control* sets the numerical parameters for the
matrix solvers used in the specified control algorithm.  Many solution
approaches require iterative solvers, and these methods enable users to
provide the maximum number of iterations and the relative tolerance.

The *control thermal* version sets the energy exchange mechanism from
the finite elements to the atoms, managed through a control algorithm.
*rescale* computes a scale factor for each atom to match the finite
element temperature.  *hoover* is a Gaussian least-constraint isokinetic
thermostat enforces that the nodal restricted atomic temperature matches
the finite element temperature.  *flux* is a similar mode, but rather
adds energy to the atoms based on conservation of energy. *hoover* and
*flux* allow the prescription of sources or fixed temperatures on the
atoms.

*correction_max_iterations* sets the maximum number of iterations to
compute the second order in time correction term for lambda with the
fractional step method. The method uses the same tolerance as the
controller's matrix solver.

Restrictions
""""""""""""

Only for be used with the specific controllers *thermal* or *momentum*.
They are ignored if a lumped solution is requested.

*control thermal* is only for be used with specific transfers: thermal (*rescale*\ , *hoover*\ , *flux*\ ), *two_temperature* (*flux*\ ).
*rescale* not valid with time filtering activated

*correction_max_iterations* is only for use with *thermal* physics using
the fractional step method.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC control momentum <atc_control_momentum>`

Default
"""""""

- *max_iterations* is the number of rows in the matrix.
- *tolerance* is 1.0e-10.
- *rescale* frequency is 1
- *flux* boundary_integration_type is *interpolate*
- *correction_max_iterations* is 20
