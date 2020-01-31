Thermodynamic output options
============================

There is one class that computes and prints thermodynamic information
to the screen and log file; see the file thermo.cpp.

There are two styles defined in thermo.cpp: "one" and "multi".  There
is also a flexible "custom" style which allows the user to explicitly
list keywords for quantities to print when thermodynamic info is
output.  See the :doc:`thermo_style <thermo_style>` command for a list
of defined quantities.

The thermo styles (one, multi, etc) are simply lists of keywords.
Adding a new style thus only requires defining a new list of keywords.
Search for the word "customize" with references to "thermo style" in
thermo.cpp to see the two locations where code will need to be added.

New keywords can also be added to thermo.cpp to compute new quantities
for output.  Search for the word "customize" with references to
"keyword" in thermo.cpp to see the several locations where code will
need to be added.

Note that the :doc:`thermo_style custom <thermo>` command already allows
for thermo output of quantities calculated by :doc:`fixes <fix>`,
:doc:`computes <compute>`, and :doc:`variables <variable>`.  Thus, it may
be simpler to compute what you wish via one of those constructs, than
by adding a new keyword to the thermo command.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
