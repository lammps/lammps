Visualize LAMMPS snapshots
==========================

LAMMPS itself does not do visualization, but snapshots from LAMMPS
simulations can be visualized (and analyzed) in a variety of ways.

Mention dump image and dump movie.

LAMMPS snapshots are created by the :doc:`dump <dump>` command which can
create files in several formats. The native LAMMPS dump format is a
text file (see "dump atom" or "dump custom") which can be visualized
by several popular visualization tools. The :doc:`dump image <dump_image>` and :doc:`dump movie <dump_image>` styles can
output internally rendered images and convert a sequence of them to a
movie during the MD run.  Several programs included with LAMMPS as
auxiliary tools can convert between LAMMPS format files and other
formats.  See the :doc:`Tools <Tools>` doc page for details.

A Python-based toolkit distributed by our group can read native LAMMPS
dump files, including custom dump files with additional columns of
user-specified atom information, and convert them to various formats
or pipe them into visualization software directly.  See the `Pizza.py WWW site <pizza_>`_ for details.  Specifically, Pizza.py can convert
LAMMPS dump files into PDB, XYZ, `Ensight <ensight_>`_, and VTK formats.
Pizza.py can pipe LAMMPS dump files directly into the Raster3d and
RasMol visualization programs.  Pizza.py has tools that do interactive
3d OpenGL visualization and one that creates SVG images of dump file
snapshots.

.. _pizza: https://pizza.sandia.gov

.. _ensight: http://www.ensight.com

.. _atomeye: http://li.mit.edu/Archive/Graphics/A/
