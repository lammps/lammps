.. index:: fix plumed

fix plumed command
==================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID plumed keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* plumed = style name of this fix command
* keyword = *plumedfile* or *outfile*
  
  .. parsed-literal::
  
       *plumedfile* arg = name of PLUMED input file to use (default: NULL)
       *outfile* arg = name of file on which to write the PLUMED log (default: NULL)



Examples
""""""""

fix pl all plumed all plumed plumedfile plumed.dat outfile p.log

Description
"""""""""""

This fix instructs LAMMPS to call the `PLUMED <plumedhome_>`_ library, which
allows one to perform various forms of trajectory analysis on the fly
and to also use methods such as umbrella sampling and metadynamics to
enhance the sampling of phase space.

The documentation included here only describes the fix plumed command
itself.  This command is LAMMPS specific, whereas most of the
functionality implemented in PLUMED will work with a range of MD codes,
and when PLUMED is used as a stand alone code for analysis.  The full
`documentation for PLUMED <plumeddocs_>`_ is available online and included
in the PLUMED source code.  The PLUMED library development is hosted at
`https://github.com/plumed/plumed2 <https://github.com/plumed/plumed2>`_
A detailed discussion of the code can be found in :ref:`(PLUMED) <PLUMED>`.

There is an example input for using this package with LAMMPS in the
examples/USER/plumed directory.


----------


The command to make LAMMPS call PLUMED during a run requires two keyword
value pairs pointing to the PLUMED input file and an output file for the
PLUMED log. The user must specify these arguments every time PLUMED is
to be used.  Furthermore, the fix plumed command should appear in the
LAMMPS input file **after** relevant input parameters (e.g. the timestep)
have been set.

The *group-ID* entry is ignored. LAMMPS will always pass all the atoms
to PLUMED and there can only be one instance of the plumed fix at a
time. The way the plumed fix is implemented ensures that the minimum
amount of information required is communicated.  Furthermore, PLUMED
supports multiple, completely independent collective variables, multiple
independent biases and multiple independent forms of analysis.  There is
thus really no restriction in functionality by only allowing only one
plumed fix in the LAMMPS input.

The *plumedfile* keyword allows the user to specify the name of the
PLUMED input file.  Instructions as to what should be included in a
plumed input file can be found in the `documentation for PLUMED <plumeddocs_>`_

The *outfile* keyword allows the user to specify the name of a file in
which to output the PLUMED log.  This log file normally just repeats the
information that is contained in the input file to confirm it was
correctly read and parsed.  The names of the files in which the results
are stored from the various analysis options performed by PLUMED will
be specified by the user in the PLUMED input file.

**Restart, fix\_modify, output, run start/stop, minimize info:**

When performing a restart of a calculation that involves PLUMED you must
include a RESTART command in the PLUMED input file as detailed in the
`PLUMED documentation <plumeddocs_>`_.  When the restart command is found in
the PLUMED input PLUMED will append to the files that were generated in
the run that was performed previously.  No part of the PLUMED restart
data is included in the LAMMPS restart files.  Furthermore, any history
dependent bias potentials that were accumulated in previous calculations
will be read in when the RESTART command is included in the PLUMED
input.

The :doc:`fix_modify <fix_modify>` *energy* option is not supported by
this fix.

Nothing is computed by this fix that can be accessed by any of the
:doc:`output commands <Howto_output>` within LAMMPS.  All the quantities
of interest can be output by commands that are native to PLUMED,
however.

Restrictions
""""""""""""


This fix is part of the USER-PLUMED package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

There can only be one plumed fix active at a time.

Related commands
""""""""""""""""

:doc:`fix smd <fix_smd>`
:doc:`fix colvars <fix_colvars>`

Default
"""""""

The default options are plumedfile = NULL and outfile = NULL


----------


.. _PLUMED:



**(PLUMED)** G.A. Tribello, M. Bonomi, D. Branduardi, C. Camilloni and G. Bussi, Comp. Phys. Comm 185, 604 (2014)

.. _plumeddocs: http://www.plumed.org/doc.html



.. _plumedhome: http://www.plumed.org/




.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
