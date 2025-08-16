Common issues that are often regarded as bugs
=============================================

The list below are some random notes on behavior of LAMMPS that is
sometimes unexpected or even considered a bug.  Most of the time, these
are just issues of understanding how LAMMPS is implemented and
parallelized.  Please also have a look at the :doc:`Error details
discussions page <Errors_details>` that contains recommendations for
tracking down issues and explanations for error messages that may
sometimes be confusing or need additional explanations.

- A LAMMPS simulation typically has two stages, 1) issuing commands
  and 2) run or minimize.  Most LAMMPS errors are detected in stage 1),
  others at the beginning of stage 2), and finally others like a bond
  stretching too far may or lost atoms or bonds may not occur until the
  middle of a run.

- If two LAMMPS runs do not produce the exact same answer on different
  machines or different numbers of processors, this is typically not a
  bug.  In theory you should get identical answers on any number of
  processors and on any machine.  In practice, numerical round-off can
  cause slight differences and eventual divergence of molecular dynamics
  phase space trajectories within a few 100s or few 1000s of timesteps.
  This can be triggered by different ordering of atoms due to different
  domain decompositions, but also through different CPU architectures,
  different operating systems, different compilers or compiler versions,
  different compiler optimization levels, different FFT libraries.
  However, the statistical properties of the two runs (e.g. average
  energy or temperature) should still be the same.

- If the :doc:`velocity <velocity>` command is used to set initial atom
  velocities, a particular atom can be assigned a different velocity
  when the problem is run on a different number of processors or on
  different machines.  If this happens, the phase space trajectories of
  the two simulations will rapidly diverge.  See the discussion of the
  *loop* option in the :doc:`velocity <velocity>` command for details
  and options that avoid this issue.

- Similarly, the :doc:`create_atoms <create_atoms>` command generates a
  lattice of atoms.  For the same physical system, the ordering and
  numbering of atoms by atom ID may be different depending on the number
  of processors.

- Some commands use random number generators which may be setup to
  produce different random number streams on each processor and hence
  will produce different effects when run on different numbers of
  processors.  A commonly-used example is the :doc:`fix langevin
  <fix_langevin>` command for thermostatting.

- LAMMPS tries to flag errors and print informative error messages so
  you can fix the problem.  For most errors it will also print the last
  input script command that it was processing or even point to the
  keyword that is causing troubles.  Of course, LAMMPS cannot figure out
  your physics or numerical mistakes, like choosing too big a timestep,
  specifying erroneous force field coefficients, or putting 2 atoms on
  top of each other!  Also, LAMMPS does not know what you *intend* to
  do, but very strictly applies the syntax as described in the
  documentation.  If you run into errors that LAMMPS does not catch that
  you think it should flag, please send an email to the `developers
  <https://www.lammps.org/authors.html>`_ or create an new topic on the
  dedicated `MatSci forum section <https://matsci.org/lammps/>`_.

- If you get an error message about an invalid command in your input
  script, you can determine what command is causing the problem by
  looking in the log.lammps file or using the :doc:`echo command <echo>`
  to see it on the screen.  If you get an error like "Invalid ...
  style", with ... being fix, compute, pair, etc, it means that you
  mistyped the style name or that the command is part of an optional
  package which was not compiled into your executable.  The list of
  available styles in your executable can be listed by using
  :doc:`the -h command-line switch <Run_options>`.  The installation and
  compilation of optional packages is explained on the :doc:`Build
  packages <Build_package>` doc page.

- For a given command, LAMMPS expects certain arguments in a specified
  order.  If you mess this up, LAMMPS will often flag the error, but it
  may also simply read a bogus argument and assign a value that is
  valid, but not what you wanted.  E.g. trying to read the string "abc"
  as an integer value of 0.  Careful reading of the associated doc page
  for the command should allow you to fix these problems. In most cases,
  where LAMMPS expects to read a number, either integer or floating
  point, it performs a stringent test on whether the provided input
  actually is an integer or floating-point number, respectively, and
  reject the input with an error message (for instance, when an integer
  is required, but a floating-point number 1.0 is provided):

  .. parsed-literal::

     ERROR: Expected integer parameter instead of '1.0' in input script or data file

- Some commands allow for using variable references in place of numeric
  constants so that the value can be evaluated and may change over the
  course of a run.  This is typically done with the syntax *v_name* for
  a parameter, where name is the name of the variable. On the other
  hand, immediate variable expansion with the syntax ${name} is
  performed while reading the input and before parsing commands,

  .. note::

     Using a variable reference (i.e. *v_name*) is only allowed if
     the documentation of the corresponding command explicitly says it is.
     Otherwise, you will receive an error message of this kind:

  .. parsed-literal::

     ERROR: Expected floating point parameter instead of 'v_name' in input script or data file

- Generally, LAMMPS will print a message to the screen and logfile and
  exit gracefully when it encounters a fatal error.  When running in
  parallel this message may be stuck in an I/O buffer and LAMMPS will be
  terminated before that buffer is printed.  In that case you can try
  adding the ``-nonblock`` or ``-nb`` command-line flag to turn off that
  buffering.  Please note that this should not be used for production
  runs, since turning off buffering usually has a significant negative
  impact on performance (even worse than :doc:`thermo_modify flush yes
  <thermo_modify>`).  Sometimes LAMMPS will print a WARNING to the
  screen and logfile and continue on; you can decide if the WARNING is
  important or not, but as a general rule do not ignore warnings that
  you not understand.  A WARNING message that is generated in the middle
  of a run is only printed to the screen, not to the logfile, to avoid
  cluttering up thermodynamic output.  If LAMMPS crashes or hangs
  without generating an error message first then it could be a bug
  (see :doc:`this section <Errors_bugs>`).

- LAMMPS runs in the available memory a processor allows to be
  allocated.  Most reasonable MD runs are compute limited, not memory
  limited, so this should not be a bottleneck on most platforms.  Almost
  all large memory allocations in the code are done via C-style malloc's
  which will generate an error message if you run out of memory.
  Smaller chunks of memory are allocated via C++ "new" statements.  If
  you are unlucky you could run out of memory just when one of these
  small requests is made, in which case the code will crash or hang (in
  parallel).

- Illegal arithmetic can cause LAMMPS to run slow or crash.  This is
  typically due to invalid physics and numerics that your simulation is
  computing.  If you see wild thermodynamic values or NaN values in your
  LAMMPS output, something is wrong with your simulation.  If you
  suspect this is happening, it is a good idea to print out
  thermodynamic info frequently (e.g. every timestep) via the
  :doc:`thermo <thermo>` so you can monitor what is happening.
  Visualizing the atom movement is also a good idea to ensure your model
  is behaving as you expect.

- When running in parallel with MPI, one way LAMMPS can hang is because
  LAMMPS has come across an error condition, but only on one or a few
  MPI processes and not all of them.  LAMMPS has two different "stop
  with an error message" functions and the correct one has to be called
  or else it will hang.
