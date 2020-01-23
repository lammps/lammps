.. index:: fix atc

fix atc command
===============

Syntax
""""""


.. parsed-literal::

   fix <fixID> <group> atc <type> <parameter_file>

* fixID = name of fix
* group = name of group fix is to be applied
* type = *thermal* or *two\_temperature* or *hardy* or *field*

.. parsed-literal::

    *thermal* = thermal coupling with fields: temperature
    *two_temperature* = electron-phonon coupling with field: temperature and electron_temperature
    *hardy* = on-the-fly post-processing using kernel localization functions (see "related" section for possible fields)
    *field* = on-the-fly post-processing using mesh-based localization functions (see "related" section for possible fields)

* parameter\_file = name of the file with material parameters. Note: Neither hardy nor field requires a parameter file


Examples
""""""""


.. parsed-literal::

   fix AtC internal atc thermal Ar_thermal.dat
   fix AtC internal atc two_temperature Ar_ttm.mat
   fix AtC internal atc hardy
   fix AtC internal atc field

Description
"""""""""""

This fix is the beginning to creating a coupled FE/MD simulation and/or an on-the-fly estimation of continuum fields. The coupled versions of this fix do Verlet integration and the post-processing does not. After instantiating this fix, several other fix\_modify commands will be needed to set up the problem, e.g. define the finite element mesh and prescribe initial and boundary conditions.

.. image:: JPG/atc_nanotube.jpg
   :align: center


.. parsed-literal::

   The following coupling example is typical, but non-exhaustive:
    # ... commands to create and initialize the MD system

    # initial fix to designate coupling type and group to apply it to
    # tag group physics material_file
    fix AtC internal atc thermal Ar_thermal.mat

    # create a uniform 12 x 2 x 2 mesh that covers region contain the group
    # nx ny nz region periodicity
    fix_modify AtC mesh create 12 2 2 mdRegion f p p

    # specify the control method for the type of coupling
    # physics control_type
    fix_modify AtC thermal control flux

    # specify the initial values for the empirical field "temperature"
    # field node_group value
    fix_modify AtC initial temperature all 30

    # create an output stream for nodal fields
    # filename output_frequency
    fix_modify AtC output atc_fe_output 100

    run 1000

likewise for this post-processing example:


.. parsed-literal::

    # ... commands to create and initialize the MD system

    # initial fix to designate post-processing and the group to apply it to
    # no material file is allowed nor required
    fix AtC internal atc hardy

    # for hardy fix, specific kernel function (function type and range) to # be used as a localization function
    fix AtC kernel quartic_sphere 10.0

    # create a uniform 1 x 1 x 1 mesh that covers region contain the group
    # with periodicity this effectively creats a system average
    fix_modify AtC mesh create 1 1 1 box p p p

    # change from default lagrangian map to eulerian
    # refreshed every 100 steps
    fix_modify AtC atom_element_map eulerian 100

    # start with no field defined
    # add mass density, potential energy density, stress and temperature
    fix_modify AtC fields add density energy stress temperature

    # create an output stream for nodal fields
    # filename output_frequency
    fix_modify AtC output nvtFE 100 text

    run 1000

the mesh's linear interpolation functions can be used as the localization function
by using the field option:


.. parsed-literal::

    fix AtC internal atc field
    fix_modify AtC mesh create 1 1 1 box p p p
    ...

Note coupling and post-processing can be combined in the same simulations using separate fixes.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  The :doc:`fix_modify <fix_modify>` options
relevant to this fix are listed below.  No global scalar or vector or
per-atom quantities are stored by this fix for access by various
:doc:`output commands <Howto_output>`.  No parameter of this fix can be
used with the *start/stop* keywords of the :doc:`run <run>` command.
This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


Thermal and two\_temperature (coupling) types use a Verlet time-integration algorithm. The hardy type does not contain its own time-integrator and must be used with a separate fix that does contain one, e.g. nve, nvt, etc. In addition, currently:

* the coupling is restricted to thermal physics
* the FE computations are done in serial on each processor.

Related commands
""""""""""""""""

After specifying this fix in your input script, several other :doc:`fix_modify <fix_modify>` commands are used to setup the problem, e.g. define the finite element mesh and prescribe initial and boundary conditions.

fix\_modify commands for setup:

* `fix\_modify AtC mesh create <USER/atc/man_mesh_create.html>`_
* `fix\_modify AtC mesh quadrature <USER/atc/man_mesh_quadrature.html>`_
* `fix\_modify AtC mesh read <USER/atc/man_mesh_read.html>`_
* `fix\_modify AtC mesh write <USER/atc/man_mesh_write.html>`_
* `fix\_modify AtC mesh create\_nodeset <USER/atc/man_mesh_create_nodeset.html>`_
* `fix\_modify AtC mesh add\_to\_nodeset <USER/atc/man_mesh_add_to_nodeset.html>`_
* `fix\_modify AtC mesh create\_faceset box <USER/atc/man_mesh_create_faceset_box.html>`_
* `fix\_modify AtC mesh create\_faceset plane <USER/atc/man_mesh_create_faceset_plane.html>`_
* `fix\_modify AtC mesh create\_elementset <USER/atc/man_mesh_create_elementset.html>`_
* `fix\_modify AtC mesh delete\_elements <USER/atc/man_mesh_delete_elements.html>`_
* `fix\_modify AtC mesh nodeset\_to\_elementset <USER/atc/man_mesh_nodeset_to_elementset.html>`_
* `fix\_modify AtC boundary <USER/atc/man_boundary.html>`_
* `fix\_modify AtC internal\_quadrature <USER/atc/man_internal_quadrature.html>`_
* `fix\_modify AtC time\_integration (thermal) <USER/atc/man_thermal_time_integration.html>`_
* `fix\_modify AtC time\_integration (momentum) <USER/atc/man_momentum_time_integration.html>`_
* `fix\_modify AtC extrinsic electron\_integration <USER/atc/man_electron_integration.html>`_
* `fix\_modify AtC internal\_element\_set <USER/atc/man_internal_element_set.html>`_
* `fix\_modify AtC decomposition <USER/atc/man_decomposition.html>`_

fix\_modify commands for boundary and initial conditions:

* `fix\_modify AtC initial <USER/atc/man_initial.html>`_
* `fix\_modify AtC fix <USER/atc/man_fix_nodes.html>`_
* `fix\_modify AtC unfix <USER/atc/man_unfix_nodes.html>`_
* `fix\_modify AtC fix\_flux <USER/atc/man_fix_flux.html>`_
* `fix\_modify AtC unfix\_flux <USER/atc/man_unfix_flux.html>`_
* `fix\_modify AtC source <USER/atc/man_source.html>`_
* `fix\_modify AtC remove\_source <USER/atc/man_remove_source.html>`_

fix\_modify commands for control and filtering:

* `fix\_modify AtC control <USER/atc/man_control.html>`_
* `fix\_modify AtC control thermal <USER/atc/man_control_thermal.html>`_
* `fix\_modify AtC control thermal correction\_max\_iterations <USER/atc/man_control_thermal_correction_max_iterations.html>`_
* `fix\_modify AtC control momentum <USER/atc/man_control_momentum.html>`_
* `fix\_modify AtC control localized\_lambda <USER/atc/man_localized_lambda.html>`_
* `fix\_modify AtC control lumped\_lambda\_solve <USER/atc/man_lumped_lambda_solve.html>`_
* `fix\_modify AtC control mask\_direction <USER/atc/man_mask_direction.html>`_ control
* `fix\_modify AtC filter <USER/atc/man_time_filter.html>`_
* `fix\_modify AtC filter scale <USER/atc/man_filter_scale.html>`_
* `fix\_modify AtC filter type <USER/atc/man_filter_type.html>`_
* `fix\_modify AtC equilibrium\_start <USER/atc/man_equilibrium_start.html>`_
* `fix\_modify AtC extrinsic exchange <USER/atc/man_extrinsic_exchange.html>`_
* `fix\_modify AtC poisson\_solver <USER/atc/man_poisson_solver.html>`_

fix\_modify commands for output:

* `fix\_modify AtC output <USER/atc/man_output.html>`_
* `fix\_modify AtC output nodeset <USER/atc/man_output_nodeset.html>`_
* `fix\_modify AtC output elementset <USER/atc/man_output_elementset.html>`_
* `fix\_modify AtC output boundary\_integral <USER/atc/man_boundary_integral.html>`_
* `fix\_modify AtC output contour\_integral <USER/atc/man_contour_integral.html>`_
* `fix\_modify AtC mesh output <USER/atc/man_mesh_output.html>`_
* `fix\_modify AtC write\_restart <USER/atc/man_write_restart.html>`_
* `fix\_modify AtC read\_restart <USER/atc/man_read_restart.html>`_

fix\_modify commands for post-processing:

* `fix\_modify AtC kernel <USER/atc/man_hardy_kernel.html>`_
* `fix\_modify AtC fields <USER/atc/man_hardy_fields.html>`_
* `fix\_modify AtC grdients <USER/atc/man_hardy_gradients.html>`_
* `fix\_modify AtC rates <USER/atc/man_hardy_rates.html>`_
* `fix\_modify AtC computes <USER/atc/man_hardy_computes.html>`_
* `fix\_modify AtC on\_the\_fly <USER/atc/man_hardy_on_the_fly.html>`_
* `fix\_modify AtC pair\_interactions/bond\_interactions <USER/atc/man_pair_interactions.html>`_
* `fix\_modify AtC sample\_frequency <USER/atc/man_sample_frequency.html>`_
* `fix\_modify AtC set <USER/atc/man_set.html>`_

miscellaneous fix\_modify commands:

* `fix\_modify AtC atom\_element\_map <USER/atc/man_atom_element_map.html>`_
* `fix\_modify AtC atom\_weight <USER/atc/man_atom_weight.html>`_
* `fix\_modify AtC write\_atom\_weights <USER/atc/man_write_atom_weights.html>`_
* `fix\_modify AtC reset\_time <USER/atc/man_reset_time.html>`_
* `fix\_modify AtC reset\_atomic\_reference\_positions <USER/atc/man_reset_atomic_reference_positions.html>`_
* `fix\_modify AtC fe\_md\_boundary <USER/atc/man_fe_md_boundary.html>`_
* `fix\_modify AtC boundary\_faceset <USER/atc/man_boundary_faceset.html>`_
* `fix\_modify AtC consistent\_fe\_initialization <USER/atc/man_consistent_fe_initialization.html>`_
* `fix\_modify AtC mass\_matrix <USER/atc/man_mass_matrix.html>`_
* `fix\_modify AtC material <USER/atc/man_material.html>`_
* `fix\_modify AtC atomic\_charge <USER/atc/man_atomic_charge.html>`_
* `fix\_modify AtC source\_integration <USER/atc/man_source_integration.html>`_
* `fix\_modify AtC temperature\_definition <USER/atc/man_temperature_definition.html>`_
* `fix\_modify AtC track\_displacement <USER/atc/man_track_displacement.html>`_
* `fix\_modify AtC boundary\_dynamics <USER/atc/man_boundary_dynamics.html>`_
* `fix\_modify AtC add\_species <USER/atc/man_add_species.html>`_
* `fix\_modify AtC add\_molecule <USER/atc/man_add_molecule.html>`_
* `fix\_modify AtC remove\_species <USER/atc/man_remove_species.html>`_
* `fix\_modify AtC remove\_molecule <USER/atc/man_remove_molecule.html>`_

Note: a set of example input files with the attendant material files are included with this package

Default
"""""""
None


----------


For detailed exposition of the theory and algorithms please see:

.. _Wagner:



**(Wagner)** Wagner, GJ; Jones, RE; Templeton, JA; Parks, MA, "An atomistic-to-continuum coupling method for heat transfer in solids." Special Issue of Computer Methods and Applied Mechanics (2008) 197:3351.

.. _Zimmeman2004:



**(Zimmerman2004)** Zimmerman, JA; Webb, EB; Hoyt, JJ;. Jones, RE; Klein, PA; Bammann, DJ, "Calculation of stress in atomistic simulation." Special Issue of Modelling and Simulation in Materials Science and Engineering (2004), 12:S319.

.. _Zimmerman2010:



**(Zimmerman2010)** Zimmerman, JA; Jones, RE; Templeton, JA, "A material frame approach for evaluating continuum variables in atomistic simulations." Journal of Computational Physics (2010), 229:2364.

.. _Templeton2010:



**(Templeton2010)** Templeton, JA; Jones, RE; Wagner, GJ, "Application of a field-based method to spatially varying thermal transport problems in molecular dynamics." Modelling and Simulation in Materials Science and Engineering (2010), 18:085007.

.. _Jones:



**(Jones)** Jones, RE; Templeton, JA; Wagner, GJ; Olmsted, D; Modine, JA, "Electron transport enhanced molecular dynamics for metals and semi-metals." International Journal for Numerical Methods in Engineering (2010), 83:940.

.. _Templeton2011:



**(Templeton2011)** Templeton, JA; Jones, RE; Lee, JW; Zimmerman, JA; Wong, BM, "A long-range electric field solver for molecular dynamics based on atomistic-to-continuum modeling." Journal of Chemical Theory and Computation (2011), 7:1736.

.. _Mandadapu:



**(Mandadapu)** Mandadapu, KK; Templeton, JA; Lee, JW, "Polarization as a field variable from molecular dynamics simulations." Journal of Chemical Physics (2013), 139:054115.

Please refer to the standard finite element (FE) texts, e.g. T.J.R Hughes " The finite element method ", Dover 2003, for the basics of FE simulation.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
