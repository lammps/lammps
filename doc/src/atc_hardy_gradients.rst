.. index:: fix_modify AtC gradients

fix_modify AtC gradients command
================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> gradients <add|delete> <list_of_fields>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* gradients = name of the AtC sub-command
* *add* or *delete* = select whether to add or delete calculation of gradients for the listed output fields
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

   fix_modify AtC gradients add temperature velocity stress
   fix_modify AtC gradients delete velocity


Description
"""""""""""

Requests calculation and output of gradients of the fields from the AtC
transfer class.  These gradients will be with regard to spatial or
material coordinate for Eulerian or Lagrangian analysis, respectively,
as specified by :doc:`fix_modify AtC atom_element_map <atc_atom_element_map>`


Restrictions
""""""""""""

Must be used with :doc:`fix atc hardy <fix_atc>`.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC atom_element_map <atc_atom_element_map>`
- :doc:`fix_modify AtC fields <atc_hardy_fields>`
- :doc:`fix_modify AtC rates <atc_hardy_rates>`

Default
"""""""

None.
