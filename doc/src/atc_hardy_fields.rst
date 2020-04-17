.. index:: fix_modify AtC fields

fix_modify AtC fields command
=============================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> fields <all|none>
   fix_modify <AtC fixID> fields <add|delete> <list_of_fields>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* fields = name of the AtC sub-command
* *all* or *none* = output all or no fields
* *add* or *delete* = add or delete the listed output fields
* list_of_fields = one or more of the fields listed below:

  - density : mass per unit volume
  - displacement : displacement vector
  - momentum : momentum per unit volume
  - velocity : defined by momentum divided by density
  - projected_velocity : simple kernel estimation of atomic velocities
  - temperature : temperature derived from the relative atomic kinetic energy
  - kinetic_temperature : temperature derived from the full kinetic energy
  - number_density : simple kernel estimation of number of atoms per unit volume
  - stress : Cauchy stress tensor for eulerian analysis (atom_element_map), or first Piola-Kirchhoff stress tensor for lagrangian analysis
  - transformed_stress : first Piola-Kirchhoff stress tensor for eulerian analysis (atom_element_map), or Cauchy stress tensor for lagrangian analysis
  - heat_flux : spatial heat flux vector for eulerian, or referential heat flux vector for lagrangian
  - potential_energy : potential energy per unit volume
  - kinetic_energy : kinetic energy per unit volume
  - thermal_energy : thermal energy (kinetic energy - continuum kinetic energy) per unit volume
  - internal_energy : total internal energy (potential + thermal) per unit volume
  - energy : total energy (potential + kinetic) per unit volume
  - number_density : number of atoms per unit volume
  - eshelby_stress : configurational stress (energy-momentum) tensor defined by [Eshelby]_
  - vacancy_concentration : volume fraction of vacancy content
  - type_concentration : volume fraction of a specific atom type
  
Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC fields add velocity temperature


Description
"""""""""""

Allows modification of the fields calculated and output by the AtC
transfer class.  The commands are cumulative, e.g.:

.. code-block:: LAMMPS

   fix_modify AtC fields none
   fix_modify AtC fields add velocity temperature

will only output the velocity and temperature fields.

Restrictions
""""""""""""

Must be used with :doc:`fix atc hardy <fix_atc>`.  Currently, the stress
and heat flux formulas are only correct for central force potentials,
e.g. Lennard-Jones and EAM but not Stillinger-Weber.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC gradients <atc_hardy_gradients>`
- :doc:`fix_modify AtC rates <atc_hardy_rates>`
- :doc:`fix_modify AtC computes <atc_hardy_computes>`


Default
"""""""

By default, no fields are output.

References
""""""""""

.. [Eshelby] J.D. Eshelby, Philos. Trans. Royal Soc. London A, Math. Phys. Sci., Vol. 244, No. 877 (1951) pp. 87-112; J. Elasticity, Vol. 5, Nos. 3-4 (1975) pp. 321-335]
