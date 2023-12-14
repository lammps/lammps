Visualize LAMMPS snapshots
==========================

Snapshots from LAMMPS simulations can be viewed, visualized, and
analyzed in a variety of ways.

LAMMPS snapshots are created by the :doc:`dump <dump>` command, which
can create files in several formats. The native LAMMPS dump format is a
text file (see "dump atom" or "dump custom") which can be visualized by
`several visualization tools <https://www.lammps.org/viz.html>`_ for MD
simulation trajectories.  `OVITO <https://www.ovito.org>`_ and `VMD
<https://www.ks.uiuc.edu/Research/vmd>`_ seem to be the most popular
choices among them.

The :doc:`dump image <dump_image>` and :doc:`dump movie <dump_image>`
styles can output internally rendered images or convert them to a movie
during the MD run.

Programs included with LAMMPS as auxiliary tools can convert
between LAMMPS format files and other formats.  See the :doc:`Tools
<Tools>` page for details.  These are rarely needed these days.
