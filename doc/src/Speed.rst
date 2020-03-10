Accelerate performance
**********************

This section describes various methods for improving LAMMPS
performance for different classes of problems running on different
kinds of machines.

There are two thrusts to the discussion that follows.  The first is
using code options that implement alternate algorithms that can
speed-up a simulation.  The second is to use one of the several
accelerator packages provided with LAMMPS that contain code optimized
for certain kinds of hardware, including multi-core CPUs, GPUs, and
Intel Xeon Phi co-processors.

The `Benchmark page <http://lammps.sandia.gov/bench.html>`_ of the LAMMPS
web site gives performance results for the various accelerator
packages discussed on the :doc:`Speed packages <Speed_packages>` doc
page, for several of the standard LAMMPS benchmark problems, as a
function of problem size and number of compute nodes, on different
hardware platforms.

.. toctree::
   :maxdepth: 1

   Speed_bench
   Speed_measure
   Speed_tips
   Speed_packages
   Speed_compare
