.. index:: fix_modify AtC control lumped_lambda_solve

fix_modify AtC control lumped_lambda_solve command
==================================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> control lumped_lambda_solve <on|off>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* control lumped_lambda_solve = name of the AtC sub-command
* *on* or *off* = Toggles state of lumped matrix

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC control lumped_lambda_solve on

Description
"""""""""""

Command select whether to use or not use lumped matrix for lambda solve.

Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""
- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

off
