Measuring performance
=====================

Before trying to make your simulation run faster, you should
understand how it currently performs and where the bottlenecks are.

The best way to do this is run the your system (actual number of
atoms) for a modest number of timesteps (say 100 steps) on several
different processor counts, including a single processor if possible.
Do this for an equilibrium version of your system, so that the
100-step timings are representative of a much longer run.  There is
typically no need to run for 1000s of timesteps to get accurate
timings; you can simply extrapolate from short runs.

For the set of runs, look at the timing data printed to the screen and
log file at the end of each LAMMPS run.  The
:doc:`Run_output <Run_output>` doc page gives an overview.

Running on one (or a few processors) should give a good estimate of
the serial performance and what portions of the timestep are taking
the most time.  Running the same problem on a few different processor
counts should give an estimate of parallel scalability.  I.e. if the
simulation runs 16x faster on 16 processors, its 100% parallel
efficient; if it runs 8x faster on 16 processors, it's 50% efficient.

The most important data to look at in the timing info is the timing
breakdown and relative percentages.  For example, trying different
options for speeding up the long-range solvers will have little impact
if they only consume 10% of the run time.  If the pairwise time is
dominating, you may want to look at GPU or OMP versions of the pair
style, as discussed below.  Comparing how the percentages change as
you increase the processor count gives you a sense of how different
operations within the timestep are scaling.  Note that if you are
running with a Kspace solver, there is additional output on the
breakdown of the Kspace time.  For PPPM, this includes the fraction
spent on FFTs, which can be communication intensive.

Another important detail in the timing info are the histograms of
atoms counts and neighbor counts.  If these vary widely across
processors, you have a load-imbalance issue.  This often results in
inaccurate relative timing data, because processors have to wait when
communication occurs for other processors to catch up.  Thus the
reported times for "Communication" or "Other" may be higher than they
really are, due to load-imbalance.  If this is an issue, you can
uncomment the MPI_Barrier() lines in src/timer.cpp, and re-compile
LAMMPS, to obtain synchronized timings.
