LAMMPS non-features
-------------------

LAMMPS is designed to be a fast, parallel engine for molecular
dynamics (MD) simulations.  It provides only a modest amount of
functionality for setting up simulations and analyzing their output.

Originally, LAMMPS was not conceived and designed for:

* being run through a GUI
* building molecular systems, or building molecular topologies
* assign force-field coefficients automagically
* perform sophisticated analysis of your MD simulation
* visualize your MD simulation interactively
* plot your output data

Over the years many of these limitations have been reduced or
removed. In part through features added to LAMMPS and in part
through external tools that either closely interface with LAMMPS
or extend LAMMPS.

Here are suggestions on how to perform these tasks:

* **GUI:** LAMMPS can be built as a library and a Python module that
  wraps the library interface is provided.  Thus, GUI interfaces can be
  written in Python or C/C++ that run LAMMPS and visualize or plot its
  output.  Examples of this are provided in the python directory and
  described on the :doc:`Python <Python_head>` doc page.  As of version
  2 August 2023 :ref:`a GUI tool <lammps_gui>` is included in LAMMPS.
  Also, there are several external wrappers or GUI front ends.
* **Builder:** Several pre-processing tools are packaged with LAMMPS.
  Some of them convert input files in formats produced by other MD codes
  such as CHARMM, AMBER, or Insight into LAMMPS input formats.  Some of
  them are simple programs that will build simple molecular systems,
  such as linear bead-spring polymer chains.  The moltemplate program is
  a true molecular builder that will generate complex molecular models.
  See the :doc:`Tools <Tools>` page for details on tools packaged with
  LAMMPS.  The `Pre-/post-processing page
  <https://www.lammps.org/prepost.html>`_ of the LAMMPS homepage
  describes a variety of third party tools for this task.  Furthermore,
  some internal LAMMPS commands allow reconstructing, or selectively adding
  topology information, as well as provide the option to insert molecule
  templates instead of atoms for building bulk molecular systems.
* **Force-field assignment:** The conversion tools described in the previous
  bullet for CHARMM, AMBER, and Insight will also assign force field
  coefficients in the LAMMPS format, assuming you provide CHARMM, AMBER,
  or BIOVIA (formerly Accelrys) force field files. The tools
  `ParmEd <https://parmed.github.io/ParmEd/html/index.html>`_ and
  `InterMol <https://github.com/shirtsgroup/InterMol>`_ are particularly
  powerful and flexible in converting force field and topology data
  between various MD simulation programs.
* **Simulation analysis:** If you want to perform analysis on-the-fly as
  your simulation runs, see the :doc:`compute <compute>` and :doc:`fix
  <fix>` doc pages, which list commands that can be used in a LAMMPS
  input script.  Also see the :doc:`Modify <Modify>` page for info on
  how to add your own analysis code or algorithms to LAMMPS.  For
  post-processing, LAMMPS output such as :doc:`dump file snapshots
  <dump>` can be converted into formats used by other MD or
  post-processing codes.  To some degree, that conversion can be done
  directly inside LAMMPS by interfacing to the VMD molfile plugins.  The
  :doc:`rerun <rerun>` command also allows post-processing of existing
  trajectories, and through being able to read a variety of file
  formats, this can also be used for analyzing trajectories from other
  MD codes.  Some post-processing tools packaged with LAMMPS will do
  these conversions.  Scripts provided in the tools/python directory can
  extract and massage data in dump files to make it easier to import
  into other programs.  See the :doc:`Tools <Tools>` page for details on
  these various options.
* **Visualization:** LAMMPS can produce NETPBM, JPG, or PNG format
  snapshot images on-the-fly via its :doc:`dump image <dump_image>`
  command and pass them to an external program, `FFmpeg  <https://ffmpeg.org/>`_,
  to generate movies from them.  For high-quality, interactive visualization,
  there are many excellent and free tools available.  See the `Visualization Tools
  <https://www.lammps.org/viz.html>`_ page of the LAMMPS website for
  visualization packages that can process LAMMPS output data.
* **Plotting:** See the next bullet about Pizza.py as well as the
  :doc:`Python <Python_head>` page for examples of plotting LAMMPS
  output.  Scripts provided with the *python* tool in the ``tools``
  directory will extract and process data in log and dump files to make
  it easier to analyze and plot.  See the :doc:`Tools <Tools>` doc page
  for more discussion of the various tools.
* **Pizza.py:** Our group has also written a separate toolkit called
  `Pizza.py <https://lammps.github.io/pizza>`_ which can do certain kinds of
  setup, analysis, plotting, and visualization (via OpenGL) for LAMMPS
  simulations.  It thus provides some functionality for several of the
  above bullets.  Pizza.py is written in `Python <https://www.python.org>`_
  and is available for download from `this page <https://sjplimp.github.io/download.html>`_.
