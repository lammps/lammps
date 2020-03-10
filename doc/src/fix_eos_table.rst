.. index:: fix eos/table

fix eos/table command
=====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID eos/table style file N keyword

* ID, group-ID are documented in :doc:`fix <fix>` command
* eos/table = style name of this fix command
* style = *linear* = method of interpolation
* file = filename containing the tabulated equation of state
* N = use N values in *linear* tables
* keyword = name of table keyword corresponding to table file

Examples
""""""""

.. parsed-literal::

   fix 1 all eos/table linear eos.table 100000 KEYWORD

Description
"""""""""""

Fix *eos/table* applies a tabulated mesoparticle equation of state to
relate the particle internal energy (u\_i) to the particle internal
temperature (dpdTheta\_i).

Fix *eos/table* creates interpolation tables of length *N* from
internal energy values listed in a file as a function of internal
temperature.

The interpolation tables are created by fitting cubic splines to the
file values and interpolating energy values at each of *N* internal
temperatures, and vice versa.  During a simulation, these tables are
used to interpolate internal energy or temperature values as needed.
The interpolation is done with the *linear* style.

For the *linear* style, the internal temperature is used to find 2
surrounding table values from which an internal energy is computed by
linear interpolation, and vice versa.

The filename specifies a file containing tabulated internal
temperature and internal energy values.  The keyword specifies a
section of the file.  The format of this file is described below.

----------

The format of a tabulated file is as follows (without the
parenthesized comments):

.. parsed-literal::

   # EOS TABLE                (one or more comment or blank lines)

   KEYWORD                    (keyword is first text on line)
   N 500                      (N  parameter)
                              (blank)
   1   1.00 0.000             (index, internal temperature, internal energy)
   2   1.02 0.001
   ...
   500 10.0 0.500

A section begins with a non-blank line whose 1st character is not a
"#"; blank lines or lines starting with "#" can be used as comments
between sections.  The first line begins with a keyword which
identifies the section.  The line can contain additional text, but the
initial text must match the argument specified in the fix command.

The next line lists the number of table entries.  The parameter "N" is
required and its value is the number of table entries that follow.
Note that this may be different than the *N* specified in the :doc:`fix eos/table <fix_eos_table>` command.  Let Ntable = *N* in the fix
command, and Nfile = "N" in the tabulated file.  What LAMMPS does is a
preliminary interpolation by creating splines using the Nfile
tabulated values as nodal points.  It uses these to interpolate as
needed to generate energy and temperature values at Ntable different
points.  The resulting tables of length Ntable are then used as
described above, when computing energy and temperature relationships.
This means that if you want the interpolation tables of length Ntable
to match exactly what is in the tabulated file (with effectively no
preliminary interpolation), you should set Ntable = Nfile.

Following a blank line, the next N lines list the tabulated values.
On each line, the 1st value is the index from 1 to N, the 2nd value is
the internal temperature (in temperature units), the 3rd value is the
internal energy (in energy units).

Note that the internal temperature and internal energy values must
increase from one line to the next.

Note that one file can contain many sections, each with a tabulated
potential.  LAMMPS reads the file section by section until it finds
one that matches the specified keyword.

----------

Restrictions
""""""""""""

This command is part of the USER-DPD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This command also requires use of the :doc:`atom_style dpd <atom_style>`
command.

The equation of state must be a monotonically increasing function.

An error will occur if the internal temperature or internal energies
are not within the table cutoffs.

Related commands
""""""""""""""""

:doc:`fix shardlow <fix_shardlow>`, :doc:`pair dpd/fdt <pair_dpd_fdt>`

**Default:** none
