Auxiliary tools
***************

LAMMPS is designed to be a computational kernel for performing
molecular dynamics computations.  Additional pre- and post-processing
steps are often necessary to setup and analyze a simulation.  A list
of such tools can be found on the `LAMMPS webpage <lws_>`_ at these links:

* `Pre/Post processing <http://lammps.sandia.gov/prepost.html>`_
* `Offsite LAMMPS packages & tools <http://lammps.sandia.gov/offsite.html>`_
* `Pizza.py toolkit <pizza_>`_

The last link for `Pizza.py <pizza_>`_ is a Python-based tool developed at
Sandia which provides tools for doing setup, analysis, plotting, and
visualization for LAMMPS simulations.

.. _lws: http://lammps.sandia.gov
.. _pizza: http://pizza.sandia.gov
.. _python: http://www.python.org



Additional tools included in the LAMMPS distribution are described on
this page.

Note that many users write their own setup or analysis tools or use
other existing codes and convert their output to a LAMMPS input format
or vice versa.  The tools listed here are included in the LAMMPS
distribution as examples of auxiliary tools.  Some of them are not
actively supported by the LAMMPS developers, as they were contributed
by LAMMPS users.  If you have problems using them, we can direct you
to the authors.

The source code for each of these codes is in the tools sub-directory
of the LAMMPS distribution.  There is a Makefile (which you may need
to edit for your platform) which will build several of the tools which
reside in that directory.  Most of them are larger packages in their
own sub-directories with their own Makefiles and/or README files.


----------


Pre-processing tools
====================

+-----------------------------+------------------------+----------------------+----------------------------------+----------------------------------+-----------------------------+
| :ref:`amber2lmp <amber>`    | :ref:`ch2lmp <charmm>` | :ref:`chain <chain>` | :ref:`createatoms <createatoms>` | :ref:`drude <drude>`             | :ref:`eam database <eamdb>` |
+-----------------------------+------------------------+----------------------+----------------------------------+----------------------------------+-----------------------------+
| :ref:`eam generate <eamgn>` | :ref:`eff <eff>`       | :ref:`ipp <ipp>`     | :ref:`micelle2d <micelle>`       | :ref:`moltemplate <moltemplate>` | :ref:`msi2lmp <msi>`        |
+-----------------------------+------------------------+----------------------+----------------------------------+----------------------------------+-----------------------------+
| :ref:`polybond <polybond>`  |                        |                      |                                  |                                  |                             |
+-----------------------------+------------------------+----------------------+----------------------------------+----------------------------------+-----------------------------+

Post-processing tools
=====================

+--------------------------+----------------------------+------------------------+--------------------------+-------------------------------+-----------------------------+
| :ref:`amber2lmp <amber>` | :ref:`binary2txt <binary>` | :ref:`ch2lmp <charmm>` | :ref:`colvars <colvars>` | :ref:`eff <eff>`              | :ref:`fep <fep>`            |
+--------------------------+----------------------------+------------------------+--------------------------+-------------------------------+-----------------------------+
| :ref:`lmp2arc <arc>`     | :ref:`lmp2cfg <cfg>`       | :ref:`matlab <matlab>` | :ref:`phonon <phonon>`   | :ref:`pymol_asphere <pymol>`  | :ref:`python <pythontools>` |
+--------------------------+----------------------------+------------------------+--------------------------+-------------------------------+-----------------------------+
| :ref:`reax <reax_tool>`  | :ref:`replica <replica>`   | :ref:`smd <smd>`       | :ref:`spin <spin>`       | :ref:`xmgrace <xmgrace>`      |                             |
+--------------------------+----------------------------+------------------------+--------------------------+-------------------------------+-----------------------------+

Miscellaneous tools
===================

+--------------------------+----------------------+-------------------+--------------------+---------------------------------------+
| :ref:`doxygen <doxygen>` | :ref:`emacs <emacs>` | :ref:`i-pi <ipi>` | :ref:`kate <kate>` | :ref:`singularity <singularity_tool>` |
+--------------------------+----------------------+-------------------+--------------------+---------------------------------------+
| :ref:`vim <vim>`         |                      |                   |                    |                                       |
+--------------------------+----------------------+-------------------+--------------------+---------------------------------------+


----------


Tool descriptions
=================

.. _amber:

amber2lmp tool
--------------------------

The amber2lmp sub-directory contains two Python scripts for converting
files back-and-forth between the AMBER MD code and LAMMPS.  See the
README file in amber2lmp for more information.

These tools were written by Keir Novik while he was at Queen Mary
University of London.  Keir is no longer there and cannot support
these tools which are out-of-date with respect to the current LAMMPS
version (and maybe with respect to AMBER as well).  Since we don't use
these tools at Sandia, you will need to experiment with them and make
necessary modifications yourself.


----------


.. _binary:

binary2txt tool
----------------------------

The file binary2txt.cpp converts one or more binary LAMMPS dump file
into ASCII text files.  The syntax for running the tool is


.. parsed-literal::

   binary2txt file1 file2 ...

which creates file1.txt, file2.txt, etc.  This tool must be compiled
on a platform that can read the binary file created by a LAMMPS run,
since binary files are not compatible across all platforms.


----------


.. _charmm:

ch2lmp tool
------------------------

The ch2lmp sub-directory contains tools for converting files
back-and-forth between the CHARMM MD code and LAMMPS.

They are intended to make it easy to use CHARMM as a builder and as a
post-processor for LAMMPS. Using charmm2lammps.pl, you can convert a
PDB file with associated CHARMM info, including CHARMM force field
data, into its LAMMPS equivalent. Support for the CMAP correction of
CHARMM22 and later is available as an option. This tool can also add
solvent water molecules and Na+ or Cl- ions to the system.
Using lammps2pdb.pl you can convert LAMMPS atom dumps into PDB files.

See the README file in the ch2lmp sub-directory for more information.

These tools were created by Pieter in't Veld (pjintve at sandia.gov)
and Paul Crozier (pscrozi at sandia.gov) at Sandia.

CMAP support added and tested by Xiaohu Hu (hux2 at ornl.gov) and
Robert A. Latour (latourr at clemson.edu), David Hyde-Volpe, and
Tigran Abramyan, (Clemson University) and
Chris Lorenz (chris.lorenz at kcl.ac.uk), King's College London.


----------


.. _chain:

chain tool
----------------------

The file chain.f creates a LAMMPS data file containing bead-spring
polymer chains and/or monomer solvent atoms.  It uses a text file
containing chain definition parameters as an input.  The created
chains and solvent atoms can strongly overlap, so LAMMPS needs to run
the system initially with a "soft" pair potential to un-overlap it.
The syntax for running the tool is


.. parsed-literal::

   chain < def.chain > data.file

See the def.chain or def.chain.ab files in the tools directory for
examples of definition files.  This tool was used to create the system
for the :doc:`chain benchmark <Speed_bench>`.


----------


.. _colvars:

colvars tools
---------------------------

The colvars directory contains a collection of tools for post-processing
data produced by the colvars collective variable library.
To compile the tools, edit the makefile for your system and run "make".

Please report problems and issues the colvars library and its tools
at: https://github.com/colvars/colvars/issues

abf\_integrate:

MC-based integration of multidimensional free energy gradient
Version 20110511


.. parsed-literal::

   Syntax: ./abf_integrate < filename > [-n < nsteps >] [-t < temp >] [-m [0\|1] (metadynamics)] [-h < hill_height >] [-f < variable_hill_factor >]

The LAMMPS interface to the colvars collective variable library, as
well as these tools, were created by Axel Kohlmeyer (akohlmey at
gmail.com) at ICTP, Italy.


----------


.. _createatoms:

createatoms tool
----------------------------------

The tools/createatoms directory contains a Fortran program called
createAtoms.f which can generate a variety of interesting crystal
structures and geometries and output the resulting list of atom
coordinates in LAMMPS or other formats.

See the included Manual.pdf for details.

The tool is authored by Xiaowang Zhou (Sandia), xzhou at sandia.gov.


----------


.. _doxygen:

doxygen tool
--------------------------

The tools/doxygen directory contains a shell script called
doxygen.sh which can generate a call graph and API lists using
the `Doxygen software <http://doxygen.org>`_.

See the included README file for details.

The tool is authored by Nandor Tamaskovics, numericalfreedom at googlemail.com.


----------


.. _drude:

drude tool
----------------------

The tools/drude directory contains a Python script called
polarizer.py which can add Drude oscillators to a LAMMPS
data file in the required format.

See the header of the polarizer.py file for details.

The tool is authored by Agilio Padua and Alain Dequidt: agilio.padua
at univ-bpclermont.fr, alain.dequidt at univ-bpclermont.fr


----------


.. _eamdb:

eam database tool
-----------------------------

The tools/eam\_database directory contains a Fortran program that will
generate EAM alloy setfl potential files for any combination of 16
elements: Cu, Ag, Au, Ni, Pd, Pt, Al, Pb, Fe, Mo, Ta, W, Mg, Co, Ti,
Zr.  The files can then be used with the :doc:`pair_style eam/alloy <pair_eam>` command.

The tool is authored by Xiaowang Zhou (Sandia), xzhou at sandia.gov,
and is based on his paper:

X. W. Zhou, R. A. Johnson, and H. N. G. Wadley, Phys. Rev. B, 69,
144113 (2004).


----------


.. _eamgn:

eam generate tool
-----------------------------

The tools/eam\_generate directory contains several one-file C programs
that convert an analytic formula into a tabulated :doc:`embedded atom method (EAM) <pair_eam>` setfl potential file.  The potentials they
produce are in the potentials directory, and can be used with the
:doc:`pair_style eam/alloy <pair_eam>` command.

The source files and potentials were provided by Gerolf Ziegenhain
(gerolf at ziegenhain.com).


----------


.. _eff:

eff tool
------------------

The tools/eff directory contains various scripts for generating
structures and post-processing output for simulations using the
electron force field (eFF).

These tools were provided by Andres Jaramillo-Botero at CalTech
(ajaramil at wag.caltech.edu).


----------


.. _emacs:

emacs tool
----------------------

The tools/emacs directory contains an Emacs Lisp add-on file for GNU Emacs
that enables a lammps-mode for editing input scripts when using GNU Emacs,
with various highlighting options set up.

These tools were provided by Aidan Thompson at Sandia
(athomps at sandia.gov).


----------


.. _fep:

fep tool
------------------

The tools/fep directory contains Python scripts useful for
post-processing results from performing free-energy perturbation
simulations using the USER-FEP package.

The scripts were contributed by Agilio Padua (Universite Blaise
Pascal Clermont-Ferrand), agilio.padua at univ-bpclermont.fr.

See README file in the tools/fep directory.


----------


.. _ipi:

i-pi tool
-------------------

The tools/i-pi directory contains a version of the i-PI package, with
all the LAMMPS-unrelated files removed.  It is provided so that it can
be used with the :doc:`fix ipi <fix_ipi>` command to perform
path-integral molecular dynamics (PIMD).

The i-PI package was created and is maintained by Michele Ceriotti,
michele.ceriotti at gmail.com, to interface to a variety of molecular
dynamics codes.

See the tools/i-pi/manual.pdf file for an overview of i-PI, and the
:doc:`fix ipi <fix_ipi>` doc page for further details on running PIMD
calculations with LAMMPS.


----------


.. _ipp:

ipp tool
------------------

The tools/ipp directory contains a Perl script ipp which can be used
to facilitate the creation of a complicated file (say, a lammps input
script or tools/createatoms input file) using a template file.

ipp was created and is maintained by Reese Jones (Sandia), rjones at
sandia.gov.

See two examples in the tools/ipp directory.  One of them is for the
tools/createatoms tool's input file.


----------


.. _kate:

kate tool
--------------------

The file in the tools/kate directory is an add-on to the Kate editor
in the KDE suite that allow syntax highlighting of LAMMPS input
scripts.  See the README.txt file for details.

The file was provided by Alessandro Luigi Sellerio
(alessandro.sellerio at ieni.cnr.it).


----------


.. _arc:

lmp2arc tool
----------------------

The lmp2arc sub-directory contains a tool for converting LAMMPS output
files to the format for Accelrys' Insight MD code (formerly
MSI/Biosym and its Discover MD code).  See the README file for more
information.

This tool was written by John Carpenter (Cray), Michael Peachey
(Cray), and Steve Lustig (Dupont).  John is now at the Mayo Clinic
(jec at mayo.edu), but still fields questions about the tool.

This tool was updated for the current LAMMPS C++ version by Jeff
Greathouse at Sandia (jagreat at sandia.gov).


----------


.. _cfg:

lmp2cfg tool
----------------------

The lmp2cfg sub-directory contains a tool for converting LAMMPS output
files into a series of \*.cfg files which can be read into the
`AtomEye <http://mt.seas.upenn.edu/Archive/Graphics/A>`_ visualizer.  See
the README file for more information.

This tool was written by Ara Kooser at Sandia (askoose at sandia.gov).


----------


.. _matlab:

matlab tool
------------------------

The matlab sub-directory contains several `MATLAB <matlabhome_>`_ scripts for
post-processing LAMMPS output.  The scripts include readers for log
and dump files, a reader for EAM potential files, and a converter that
reads LAMMPS dump files and produces CFG files that can be visualized
with the `AtomEye <http://mt.seas.upenn.edu/Archive/Graphics/A>`_
visualizer.

See the README.pdf file for more information.

These scripts were written by Arun Subramaniyan at Purdue Univ
(asubrama at purdue.edu).

.. _matlabhome: http://www.mathworks.com




----------


.. _micelle:

micelle2d tool
----------------------------

The file micelle2d.f creates a LAMMPS data file containing short lipid
chains in a monomer solution.  It uses a text file containing lipid
definition parameters as an input.  The created molecules and solvent
atoms can strongly overlap, so LAMMPS needs to run the system
initially with a "soft" pair potential to un-overlap it.  The syntax
for running the tool is


.. parsed-literal::

   micelle2d < def.micelle2d > data.file

See the def.micelle2d file in the tools directory for an example of a
definition file.  This tool was used to create the system for the
:doc:`micelle example <Examples>`.


----------


.. _moltemplate:

moltemplate tool
----------------------------------

The moltemplate sub-directory contains instructions for installing
moltemplate, a Python-based tool for building molecular systems based
on a text-file description, and creating LAMMPS data files that encode
their molecular topology as lists of bonds, angles, dihedrals, etc.
See the README.txt file for more information.

This tool was written by Andrew Jewett (jewett.aij at gmail.com), who
supports it.  It has its own WWW page at
`http://moltemplate.org <http://moltemplate.org>`_.
The latest sources can be found `on its GitHub page <https://github.com/jewettaij/moltemplate/releases>`_


----------


.. _msi:

msi2lmp tool
----------------------

The msi2lmp sub-directory contains a tool for creating LAMMPS template
input and data files from BIOVIA's Materias Studio files (formerly
Accelrys' Insight MD code, formerly MSI/Biosym and its Discover MD code).

This tool was written by John Carpenter (Cray), Michael Peachey
(Cray), and Steve Lustig (Dupont). Several people contributed changes
to remove bugs and adapt its output to changes in LAMMPS.

This tool has several known limitations and is no longer under active
development, so there are no changes except for the occasional bug fix.

See the README file in the tools/msi2lmp folder for more information.


----------


.. _phonon:

phonon tool
------------------------

The phonon sub-directory contains a post-processing tool useful for
analyzing the output of the :doc:`fix phonon <fix_phonon>` command in
the USER-PHONON package.

See the README file for instruction on building the tool and what
library it needs.  And see the examples/USER/phonon directory
for example problems that can be post-processed with this tool.

This tool was written by Ling-Ti Kong at Shanghai Jiao Tong
University.


----------


.. _polybond:

polybond tool
----------------------------

The polybond sub-directory contains a Python-based tool useful for
performing "programmable polymer bonding".  The Python file
lmpsdata.py provides a "Lmpsdata" class with various methods which can
be invoked by a user-written Python script to create data files with
complex bonding topologies.

See the Manual.pdf for details and example scripts.

This tool was written by Zachary Kraus at Georgia Tech.


----------


.. _pymol:

pymol\_asphere tool
-------------------------------

The pymol\_asphere sub-directory contains a tool for converting a
LAMMPS dump file that contains orientation info for ellipsoidal
particles into an input file for the `PyMol visualization package <pymolhome_>`_ or its `open source variant <pymolopen_>`_.

.. _pymolhome: http://www.pymol.org



.. _pymolopen: http://sourceforge.net/scm/?type=svn&group\_id=4546



Specifically, the tool triangulates the ellipsoids so they can be
viewed as true ellipsoidal particles within PyMol.  See the README and
examples directory within pymol\_asphere for more information.

This tool was written by Mike Brown at Sandia.


----------


.. _pythontools:

python tool
-----------------------------

The python sub-directory contains several Python scripts
that perform common LAMMPS post-processing tasks, such as:

* extract thermodynamic info from a log file as columns of numbers
* plot two columns of thermodynamic info from a log file using GnuPlot
* sort the snapshots in a dump file by atom ID
* convert multiple :doc:`NEB <neb>` dump files into one dump file for viz
* convert dump files into XYZ, CFG, or PDB format for viz by other packages

These are simple scripts built on `Pizza.py <pizza_>`_ modules.  See the
README for more info on Pizza.py and how to use these scripts.


----------


.. _replica:

replica tool
--------------------------

The tools/replica directory contains the reorder\_remd\_traj python script which
can be used to reorder the replica trajectories (resulting from the use of the 
temper command) according to temperature. This will produce discontinuous
trajectories with all frames at the same temperature in each trajectory.
Additional options can be used to calculate the canonical configurational
log-weight for each frame at each temperature using the pymbar package. See
the README.md file for further details. Try out the peptide example provided.

This tool was written by (and is maintained by) Tanmoy Sanyal, 
while at the Shell lab at UC Santa Barbara. (tanmoy dot 7989 at gmail.com)


----------


.. _reax\_tool:

reax tool
--------------------------

The reax sub-directory contains stand-alone codes that can
post-process the output of the :doc:`fix reax/c/bonds <fix_reaxc_bonds>`
command from a LAMMPS simulation using :doc:`ReaxFF <pair_reaxc>`.  See
the README.txt file for more info.

These tools were written by Aidan Thompson at Sandia.


----------


.. _smd:

smd tool
------------------

The smd sub-directory contains a C++ file dump2vtk\_tris.cpp and
Makefile which can be compiled and used to convert triangle output
files created by the Smooth-Mach Dynamics (USER-SMD) package into a
VTK-compatible unstructured grid file.  It could then be read in and
visualized by VTK.

See the header of dump2vtk.cpp for more details.

This tool was written by the USER-SMD package author, Georg
Ganzenmuller at the Fraunhofer-Institute for High-Speed Dynamics,
Ernst Mach Institute in Germany (georg.ganzenmueller at emi.fhg.de).


----------


.. _spin:

spin tool
--------------------

The spin sub-directory contains a C file interpolate.c which can
be compiled and used to perform a cubic polynomial interpolation of
the MEP following a GNEB calculation.

See the README file in tools/spin/interpolate\_gneb for more details.

This tool was written by the SPIN package author, Julien
Tranchida at Sandia National Labs (jtranch at sandia.gov, and by Aleksei
Ivanov, at University of Iceland (ali5 at hi.is).


----------


.. _singularity\_tool:

singularity tool
----------------------------------------

The singularity sub-directory contains container definitions files
that can be used to build container images for building and testing
LAMMPS on specific OS variants using the `Singularity <https://sylabs.io>`_
container software. Contributions for additional variants are welcome.


----------


.. _vim:

vim tool
------------------

The files in the tools/vim directory are add-ons to the VIM editor
that allow easier editing of LAMMPS input scripts.  See the README.txt
file for details.

These files were provided by Gerolf Ziegenhain (gerolf at
ziegenhain.com)


----------


.. _xmgrace:

xmgrace tool
--------------------------

The files in the tools/xmgrace directory can be used to plot the
thermodynamic data in LAMMPS log files via the xmgrace plotting
package.  There are several tools in the directory that can be used in
post-processing mode.  The lammpsplot.cpp file can be compiled and
used to create plots from the current state of a running LAMMPS
simulation.

See the README file for details.

These files were provided by Vikas Varshney (vv0210 at gmail.com)
