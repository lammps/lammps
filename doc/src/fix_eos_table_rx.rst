.. index:: fix eos/table/rx

fix eos/table/rx command
========================

fix eos/table/rx/kk command
===========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID eos/table/rx style file1 N keyword ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* eos/table/rx = style name of this fix command
* style = *linear* = method of interpolation
* file1 = filename containing the tabulated equation of state
* N = use N values in *linear* tables
* keyword = name of table keyword corresponding to table file
* file2 = filename containing the heats of formation of each species (optional)
* deltaHf = heat of formation for a single species in energy units (optional)
* energyCorr = energy correction in energy units (optional)
* tempCorrCoeff = temperature correction coefficient (optional)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all eos/table/rx linear eos.table 10000 KEYWORD thermo.table
   fix 1 all eos/table/rx linear eos.table 10000 KEYWORD 1.5
   fix 1 all eos/table/rx linear eos.table 10000 KEYWORD 1.5 0.025 0.0

Description
"""""""""""

Fix *eos/table/rx* applies a tabulated mesoparticle equation
of state to relate the concentration-dependent particle internal
energy (:math:`u_i`) to the particle internal temperature (:math:`\theta_i`).

The concentration-dependent particle internal energy (:math:`u_i`) is
computed according to the following relation:

.. math::

   U_{i} = \displaystyle\sum_{j=1}^{m} c_{i,j}(u_{j} + \Delta H_{f,j}) + \frac{3k_{b}T}{2} + Nk_{b}T \\

where *m* is the number of species, :math:`c_{i,j}` is the
concentration of species *j* in particle *i*\ , :math:`u_j` is the
internal energy of species j, :math:`\Delta H_{f,j} is the heat of
formation of species *j*\ , N is the number of molecules represented
by the coarse-grained particle, :math:`k_b` is the Boltzmann constant,
and *T* is the temperature of the system.  Additionally, it is
possible to modify the concentration-dependent particle internal
energy relation by adding an energy correction, temperature-dependent
correction, and/or a molecule-dependent correction.  An energy
correction can be specified as a constant (in energy units).  A
temperature correction can be specified by multiplying a temperature
correction coefficient by the internal temperature.  A molecular
correction can be specified by by multiplying a molecule correction
coefficient by the average number of product gas particles in the
coarse-grain particle.

Fix *eos/table/rx* creates interpolation tables of length *N* from *m*
internal energy values of each species :math:`u_j` listed in a file as a
function of internal temperature.  During a simulation, these tables
are used to interpolate internal energy or temperature values as needed.
The interpolation is done with the *linear* style.  For the *linear* style,
the internal temperature is used to find 2 surrounding table values from
which an internal energy is computed by linear interpolation.  A secant
solver is used to determine the internal temperature from the internal energy.

The first filename specifies a file containing tabulated internal
temperature and *m* internal energy values for each species :math:`u_j`.
The keyword specifies a section of the file.  The format of this
file is described below.

The second filename specifies a file containing heat of formation
:math:`\Delta H_{f,j}` for each species.

In cases where the coarse-grain particle represents a single molecular
species (i.e., no reactions occur and fix *rx* is not present in the
input file), fix *eos/table/rx* can be applied in a similar manner to
fix *eos/table* within a non-reactive DPD simulation.  In this case,
the heat of formation filename is replaced with the heat of formation
value for the single species.  Additionally, the energy correction and
temperature correction coefficients may also be specified as fix
arguments.

----------

The format of a tabulated file is as follows (without the
parenthesized comments):

.. parsed-literal::

   # EOS TABLE                (one or more comment or blank lines)

   KEYWORD                    (keyword is first text on line)
   N 500 h2 no2 n2 ... no     (N  parameter species1 species2 ... speciesN)
                              (blank)
   1   1.00 0.000 ... 0.0000  (index, internal temperature, internal energy of species 1, ..., internal energy of species m)
   2   1.02 0.001 ... 0.0002
   ...
   500 10.0 0.500 ... 1.0000

A section begins with a non-blank line whose 1st character is not a
"#"; blank lines or lines starting with "#" can be used as comments
between sections.  The first line begins with a keyword which
identifies the section.  The line can contain additional text, but the
initial text must match the argument specified in the fix command.

The next line lists the number of table entries and the species names
that correspond with all the species listed in the reaction equations
through the *fix rx* command.
The parameter "N" is required and its value is the number of table
entries that follow.  Let Nfile = "N" in the tabulated file.
What LAMMPS does is a preliminary interpolation by creating splines
using the Nfile tabulated values as nodal points.

Following a blank line, the next N lines list the tabulated values.
On each line, the 1st value is the index from 1 to N, the 2nd value is
the internal temperature (in temperature units), the 3rd value until
the *m+3* value are the internal energies of the m species (in energy units).

Note that all internal temperature and internal energy values must
increase from one line to the next.

Note that one file can contain many sections, each with a tabulated
potential.  LAMMPS reads the file section by section until it finds
one that matches the specified keyword.

----------

The format of a heat of formation file is as follows (without the
parenthesized comments):

.. parsed-literal::

   # HEAT OF FORMATION TABLE  (one or more comment or blank lines)

                              (blank)
   h2      0.00               (species name, heat of formation)
   no2     0.34
   n2      0.00
   ...
   no      0.93

Note that the species can be listed in any order.  The tag that is
used as the species name must correspond with the tags used to define
the reactions with the :doc:`fix rx <fix_rx>` command.

Alternatively, corrections to the EOS can be included by specifying
three additional columns that correspond to the energy correction,
the temperature correction coefficient and molecule correction
coefficient.  In this case, the format of the file is as follows:

.. parsed-literal::

   # HEAT OF FORMATION TABLE     (one or more comment or blank lines)

                                 (blank)
   h2      0.00 1.23 0.025  0.0  (species name, heat of formation, energy correction, temperature correction coefficient, molecule correction coefficient)
   no2     0.34 0.00 0.000 -1.76
   n2      0.00 0.00 0.000 -1.76
   ...
   no      0.93 0.00 0.000 -1.76

----------

Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.

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

:doc:`fix rx <fix_rx>`,
:doc:`pair dpd/fdt <pair_dpd_fdt>`

**Default:** none
