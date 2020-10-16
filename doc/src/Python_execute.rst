Executing commands
==================

Once an instance of the :py:class:`lammps <lammps.lammps>`,
:py:class:`PyLammps <lammps.PyLammps>`, or
:py:class:`IPyLammps <lammps.IPyLammps>` class is created, there are
multiple ways to "feed" it commands. In a way that is not very different from
running a LAMMPS input script, except that Python has many more facilities
for structured programming than the LAMMPS input script syntax. Furthermore
it is possible to "compute" what the next LAMMPS command should be.

.. tabs::

   .. tab:: lammps API

      Same as in the equivalent
      :doc:`C library functions <Library_execute>`, commands can be read from a file, a
      single string, a list of strings and a block of commands in a single
      multi-line string. They are processed under the same boundary conditions
      as the C library counterparts.  The example below demonstrates the use
      of :py:func:`lammps.file()`, :py:func:`lammps.command()`,
      :py:func:`lammps.commands_list()`, and :py:func:`lammps.commands_string()`:

      .. code-block:: python

         from lammps import lammps
         lmp = lammps()

         # read commands from file 'in.melt'
         lmp.file('in.melt')

         # issue a single command
         lmp.command('variable zpos index 1.0')

         # create 10 groups with 10 atoms each
         cmds = ["group g{} id {}:{}".format(i,10*i+1,10*(i+1)) for i in range(10)]
         lmp.commands_list(cmds)

         # run commands from a multi-line string
         block = """
         clear
         region  box block 0 2 0 2 0 2
         create_box 1 box
         create_atoms 1 single 1.0 1.0 ${zpos}
         """
         lmp.commands_string(block)

   .. tab:: PyLammps/IPyLammps API

      Unlike the lammps API, the PyLammps/IPyLammps APIs allow running LAMMPS
      commands by calling equivalent member functions of :py:class:`PyLammps <lammps.PyLammps>`
      and :py:class:`IPyLammps <lammps.IPyLammps>` instances.

      For instance, the following LAMMPS command

      .. code-block:: LAMMPS

         region box block 0 10 0 5 -0.5 0.5

      can be executed using with the lammps API with the following Python code if ``lmp`` is an
      instance of :py:class:`lammps <lammps.lammps>`:

      .. code-block:: Python

         from lammps import lammps

         lmp = lammps()
         lmp.command("region box block 0 10 0 5 -0.5 0.5")

      With the PyLammps interface, any LAMMPS command can be split up into arbitrary parts.
      These parts are then passed to a member function with the name of the :doc:`command <Commands_all>`.
      For the :doc:`region <region>` command that means the :code:`region()` method can be called.
      The arguments of the command can be passed as one string, or
      individually.

      .. code-block:: Python

         from lammps import PyLammps

         L = PyLammps()

         # pass command parameters as one string
         L.region("box block 0 10 0 5 -0.5 0.5")

         # OR pass them individually
         L.region("box block", 0, 10, 0, 5, -0.5, 0.5)

      In the latter example, all parameters except the first are Python floating-point literals. The
      member function takes the entire parameter list and transparently merges it to a single command
      string.

      The benefit of this approach is avoiding redundant command calls and easier
      parameterization. In the lammps API parameterization needed to be done
      manually by creating formatted command strings.

      .. code-block:: Python

         lmp.command("region box block %f %f %f %f %f %f" % (xlo, xhi, ylo, yhi, zlo, zhi))

      In contrast, methods of PyLammps accept parameters directly and will convert
      them automatically to a final command string.

      .. code-block:: Python

         L.region("box block", xlo, xhi, ylo, yhi, zlo, zhi)

      Using these facilities, the example shown for the lammps API can be rewritten as follows:

      .. code-block:: python

         from lammps import PyLammps
         L = PyLammps()

         # read commands from file 'in.melt'
         L.file('in.melt')

         # issue a single command
         L.variable('zpos', 'index', 1.0)

         # create 10 groups with 10 atoms each
         for i in range(10):
            L.group(f"g{i}", "id", f"{10*i+1}:{10*(i+1)}")

         L.clear()
         L.region("box block", 0, 2, 0, 2, 0, 2)
         L.create_box(1, "box")
         L.create_atoms(1, "single", 1.0, 1.0, "${zpos}")
