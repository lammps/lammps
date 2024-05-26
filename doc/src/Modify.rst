Modifying & extending LAMMPS
****************************

LAMMPS has a modular design, so that it is easy to modify or extend with
new functionality.  In fact, about 95% of its source code is optional.
The following pages give basic instructions on adding new features to
LAMMPS.  More in-depth explanations and documentation of individual
functions and classes are given in :doc:`Developer`.

If you add a new feature to LAMMPS and think it will be of general
interest to other users, we encourage you to submit it for inclusion in
LAMMPS. This process is explained in the following three pages:

* :doc:`how to prepare and submit your code <Modify_contribute>`
* :doc:`requirements for submissions <Modify_requirements>`
* :doc:`style guidelines <Modify_style>`

A summary description of various types of styles in LAMMPS follows.
A discussion of implementing specific styles from scratch is given
in :doc:`writing new styles <Developer_write>`.

.. toctree::
   :maxdepth: 1

   Modify_overview
   Modify_contribute
   Modify_requirements
   Modify_style

.. toctree::
   :maxdepth: 1

   Modify_atom
   Modify_pair
   Modify_bond
   Modify_compute
   Modify_fix
   Modify_command
   Modify_dump
   Modify_kspace
   Modify_min
   Modify_region
   Modify_body
   Modify_gran_sub_mod
   Modify_thermo
   Modify_variable
