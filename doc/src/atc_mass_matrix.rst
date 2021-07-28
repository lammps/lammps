.. index:: fix_modify AtC mass_matrix

fix_modify AtC mass_matrix command
==================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> mass_matrix <fe|md_fe>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* mass_matrix = name of the AtC sub-command
* *fe* or *md_fe* = activate/deactivate using the FE mass matrix in the MD region

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC mass_matrix fe

Description
"""""""""""

Determines whether AtC uses the FE mass matrix based on Gaussian
quadrature or based on atomic quadrature in the MD region. This is
useful for fully overlapping simulations to improve efficiency.

Restrictions
""""""""""""

Should not be used unless the FE region is contained within the MD
region, otherwise the method will be unstable and inaccurate.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

*md_fe*
