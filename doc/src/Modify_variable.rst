Variable options
================

There is one class that computes and stores :doc:`variable <variable>`
information in LAMMPS; see the file variable.cpp.  The value
associated with a variable can be periodically printed to the screen
via the :doc:`print <print>`, :doc:`fix print <fix_print>`, or
:doc:`thermo_style custom <thermo_style>` commands.  Variables of style
"equal" can compute complex equations that involve the following types
of arguments:

.. parsed-literal::

   thermo keywords = ke, vol, atoms, ...
   other variables = v_a, v_myvar, ...
   math functions = div(x,y), mult(x,y), add(x,y), ...
   group functions = mass(group), xcm(group,x), ...
   atom values = x[123], y[3], vx[34], ...
   compute values = c_mytemp[0], c_thermo_press[3], ...

Adding keywords for the :doc:`thermo_style custom <thermo_style>`
command (which can then be accessed by variables) is discussed on the
:doc:`Modify thermo <Modify_thermo>` doc page.

Adding a new math function of one or two arguments can be done by
editing one section of the Variable::evaluate() method.  Search for
the word "customize" to find the appropriate location.

Adding a new group function can be done by editing one section of the
Variable::evaluate() method.  Search for the word "customize" to find
the appropriate location.  You may need to add a new method to the
Group class as well (see the group.cpp file).

Accessing a new atom-based vector can be done by editing one section
of the Variable::evaluate() method.  Search for the word "customize"
to find the appropriate location.

Adding new :doc:`compute styles <compute>` (whose calculated values can
then be accessed by variables) is discussed on the :doc:`Modify compute <Modify_compute>` doc page.
