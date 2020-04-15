Calculate diffusion coefficients
================================

The diffusion coefficient D of a material can be measured in at least
2 ways using various options in LAMMPS.  See the examples/DIFFUSE
directory for scripts that implement the 2 methods discussed here for
a simple Lennard-Jones fluid model.

The first method is to measure the mean-squared displacement (MSD) of
the system, via the :doc:`compute msd <compute_msd>` command.  The slope
of the MSD versus time is proportional to the diffusion coefficient.
The instantaneous MSD values can be accumulated in a vector via the
:doc:`fix vector <fix_vector>` command, and a line fit to the vector to
compute its slope via the :doc:`variable slope <variable>` function, and
thus extract D.

The second method is to measure the velocity auto-correlation function
(VACF) of the system, via the :doc:`compute vacf <compute_vacf>`
command.  The time-integral of the VACF is proportional to the
diffusion coefficient.  The instantaneous VACF values can be
accumulated in a vector via the :doc:`fix vector <fix_vector>` command,
and time integrated via the :doc:`variable trap <variable>` function,
and thus extract D.
