Modifying & extending LAMMPS
****************************

LAMMPS is designed in a modular fashion and to be easy to modify or
extend with new functionality.  In fact, about 95% of its source code
are optional.  The following pages give basic instructions on what
is required when adding new styles of different kinds to LAMMPS.

If you add a new feature to LAMMPS and think it will be of general
interest to other users, we encourage you to submit it for inclusion in
LAMMPS as a pull request on our `GitHub site
<https://github.com/lammps/lammps>`_, after reading about :doc:`how to
prepare your code for submission <Modify_contribute>` and :doc:`the
style requirements and recommendations <Modify_style>`.

.. toctree::
   :maxdepth: 1

   Modify_overview
   Modify_contribute
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
