.. index:: fix_modify AtC control localized_lambda

fix_modify AtC control localized_lambda command
===============================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> control localized_lambda <on|off>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* control localized_lambda = name of the AtC sub-command
* *on* or *off* = Toggles state of localization algorithm

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC control localized_lambda on

Description
"""""""""""

Turns the localization algorithms *on* or *off* for control algorithms
to restrict the influence of FE coupling or boundary conditions to a
region near the boundary of the MD region.  Control algorithms will not
affect atoms in elements not possessing faces on the boundary of the
region.  Flux-based control is localized via row-sum lumping while
quantity control is done by solving a truncated matrix equation.

Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""
- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

off
