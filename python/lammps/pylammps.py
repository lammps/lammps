# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

################################################################################
# Alternative Python Wrapper
# Written by Richard Berger <richard.berger@temple.edu>
################################################################################

# for python2/3 compatibility

from __future__ import print_function

import os
import re
import select
from collections import namedtuple

from .core import lammps


class OutputCapture(object):
  """ Utility class to capture LAMMPS library output """

  def __init__(self):
    self.stdout_pipe_read, self.stdout_pipe_write = os.pipe()
    self.stdout_fd = 1

  def __enter__(self):
    self.stdout = os.dup(self.stdout_fd)
    os.dup2(self.stdout_pipe_write, self.stdout_fd)
    return self

  def __exit__(self, exc_type, exc_value, traceback):
    os.dup2(self.stdout, self.stdout_fd)
    os.close(self.stdout)
    os.close(self.stdout_pipe_read)
    os.close(self.stdout_pipe_write)

  # check if we have more to read from the pipe
  def more_data(self, pipe):
    r, _, _ = select.select([pipe], [], [], 0)
    return bool(r)

  # read the whole pipe
  def read_pipe(self, pipe):
    out = ""
    while self.more_data(pipe):
      out += os.read(pipe, 1024).decode()
    return out

  @property
  def output(self):
    return self.read_pipe(self.stdout_pipe_read)

# -------------------------------------------------------------------------

class Variable(object):
  def __init__(self, pylammps_instance, name, style, definition):
    self._pylmp = pylammps_instance
    self.name = name
    self.style = style
    self.definition = definition.split()

  @property
  def value(self):
    if self.style == 'atom':
      return list(self._pylmp.lmp.extract_variable(self.name, "all", 1))
    else:
      value = self._pylmp.lmp_print('"${%s}"' % self.name).strip()
      try:
        return float(value)
      except ValueError:
        return value

# -------------------------------------------------------------------------

class AtomList(object):
  """
  A dynamic list of atoms that returns either an :py:class:`Atom` or
  :py:class:`Atom2D` instance for each atom. Instances are only allocated
  when accessed.

  :ivar natoms: total number of atoms
  :ivar dimensions: number of dimensions in system
  """
  def __init__(self, pylammps_instance):
    self._pylmp = pylammps_instance
    self.natoms = self._pylmp.system.natoms
    self.dimensions = self._pylmp.system.dimensions
    self._loaded = {}

  def __getitem__(self, index):
    """
    Return Atom with given local index

    :param index: Local index of atom
    :type index: int
    :rtype: Atom or Atom2D
    """
    if index not in self._loaded:
        if self.dimensions == 2:
            atom = Atom2D(self._pylmp, index + 1)
        else:
            atom = Atom(self._pylmp, index + 1)
        self._loaded[index] = atom
    return self._loaded[index]

  def __len__(self):
    return self.natoms


# -------------------------------------------------------------------------

class Atom(object):
  """
  A wrapper class then represents a single atom inside of LAMMPS

  It provides access to properties of the atom and allows you to change some of them.
  """
  def __init__(self, pylammps_instance, index):
    self._pylmp = pylammps_instance
    self.index = index

  def __dir__(self):
    return [k for k in super().__dir__() if not k.startswith('_')]

  @property
  def id(self):
    """
    Return the atom ID

    :type: int
    """
    return int(self._pylmp.eval("id[%d]" % self.index))

  @property
  def type(self):
    """
    Return the atom type

    :type: int
    """
    return int(self._pylmp.eval("type[%d]" % self.index))

  @property
  def mol(self):
    """
    Return the atom molecule index

    :type: int
    """
    return self._pylmp.eval("mol[%d]" % self.index)

  @property
  def mass(self):
    """
    Return the atom mass

    :type: float
    """
    return self._pylmp.eval("mass[%d]" % self.index)

  @property
  def position(self):
    """
    :getter: Return position of atom
    :setter: Set position of atom
    :type: tuple (float, float, float)
    """
    return (self._pylmp.eval("x[%d]" % self.index),
            self._pylmp.eval("y[%d]" % self.index),
            self._pylmp.eval("z[%d]" % self.index))

  @position.setter
  def position(self, value):
    """
    :getter: Return velocity of atom
    :setter: Set velocity of atom
    :type: tuple (float, float, float)
    """
    self._pylmp.set("atom", self.index, "x", value[0])
    self._pylmp.set("atom", self.index, "y", value[1])
    self._pylmp.set("atom", self.index, "z", value[2])

  @property
  def velocity(self):
    return (self._pylmp.eval("vx[%d]" % self.index),
            self._pylmp.eval("vy[%d]" % self.index),
            self._pylmp.eval("vz[%d]" % self.index))

  @velocity.setter
  def velocity(self, value):
     self._pylmp.set("atom", self.index, "vx", value[0])
     self._pylmp.set("atom", self.index, "vy", value[1])
     self._pylmp.set("atom", self.index, "vz", value[2])

  @property
  def force(self):
    """
    Return the total force acting on the atom

    :type: tuple (float, float, float)
    """
    return (self._pylmp.eval("fx[%d]" % self.index),
            self._pylmp.eval("fy[%d]" % self.index),
            self._pylmp.eval("fz[%d]" % self.index))

  @property
  def charge(self):
    """
    Return the atom charge

    :type: float
    """
    return self._pylmp.eval("q[%d]" % self.index)

# -------------------------------------------------------------------------

class Atom2D(Atom):
  """
  A wrapper class then represents a single 2D atom inside of LAMMPS

  Inherits all properties from the :py:class:`Atom` class, but returns 2D versions
  of position, velocity, and force.

  It provides access to properties of the atom and allows you to change some of them.
  """
  def __init__(self, pylammps_instance, index):
    super(Atom2D, self).__init__(pylammps_instance, index)

  @property
  def position(self):
    """Access to coordinates of an atom

    :getter: Return position of atom
    :setter: Set position of atom
    :type: tuple (float, float)
    """
    return (self._pylmp.eval("x[%d]" % self.index),
            self._pylmp.eval("y[%d]" % self.index))

  @position.setter
  def position(self, value):
     self._pylmp.set("atom", self.index, "x", value[0])
     self._pylmp.set("atom", self.index, "y", value[1])

  @property
  def velocity(self):
    """Access to velocity of an atom
    :getter: Return velocity of atom
    :setter: Set velocity of atom
    :type: tuple (float, float)
    """
    return (self._pylmp.eval("vx[%d]" % self.index),
            self._pylmp.eval("vy[%d]" % self.index))

  @velocity.setter
  def velocity(self, value):
     self._pylmp.set("atom", self.index, "vx", value[0])
     self._pylmp.set("atom", self.index, "vy", value[1])

  @property
  def force(self):
    """Access to force of an atom

    :type: tuple (float, float)
    """
    return (self._pylmp.eval("fx[%d]" % self.index),
            self._pylmp.eval("fy[%d]" % self.index))

# -------------------------------------------------------------------------

class variable_set:
    def __init__(self, name, variable_dict):
        self._name = name
        array_pattern = re.compile(r"(?P<arr>.+)\[(?P<index>[0-9]+)\]")

        for key, value in variable_dict.items():
            m = array_pattern.match(key)
            if m:
                g = m.groupdict()
                varname = g['arr']
                idx = int(g['index'])
                if varname not in self.__dict__:
                    self.__dict__[varname] = {}
                self.__dict__[varname][idx] = value
            else:
                self.__dict__[key] = value

    def __str__(self):
        return "{}({})".format(self._name, ','.join(["{}={}".format(k, self.__dict__[k]) for k in self.__dict__.keys() if not k.startswith('_')]))

    def __dir__(self):
        return [k for k in self.__dict__.keys() if not k.startswith('_')]

    def __repr__(self):
        return self.__str__()

# -------------------------------------------------------------------------

def get_thermo_data(output):
    """ traverse output of runs and extract thermo data columns """
    if isinstance(output, str):
        lines = output.splitlines()
    else:
        lines = output

    runs = []
    columns = []
    in_run = False
    current_run = {}

    for line in lines:
        if line.startswith("Per MPI rank memory allocation"):
            in_run = True
        elif in_run and len(columns) == 0:
            # first line after memory usage are column names
            columns = line.split()

            current_run = {}

            for col in columns:
                current_run[col] = []

        elif line.startswith("Loop time of "):
            in_run = False
            columns = []
            thermo_data = variable_set('ThermoData', current_run)
            r = {'thermo' : thermo_data }
            runs.append(namedtuple('Run', list(r.keys()))(*list(r.values())))
        elif in_run and len(columns) > 0:
            items = line.split()
            # Convert thermo output and store it.
            # It must have the same number of columns and
            # all of them must be convertible to floats.
            # Otherwise we ignore the line
            if len(items) == len(columns):
                try:
                    values = [float(x) for x in items]
                    for i, col in enumerate(columns):
                        current_run[col].append(values[i])
                except ValueError:
                  # cannot convert. must be a non-thermo output. ignore.
                  pass

    return runs

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

class PyLammps(object):
  """
  This is a Python wrapper class around the lower-level
  :py:class:`lammps` class, exposing a more Python-like,
  object-oriented interface for prototyping system inside of IPython and
  Jupyter notebooks.

  It either creates its own instance of :py:class:`lammps` or can be
  initialized with an existing instance. The arguments are the same of the
  lower-level interface. The original interface can still be accessed via
  :py:attr:`PyLammps.lmp`.

  :param name: "machine" name of the shared LAMMPS library ("mpi" loads ``liblammps_mpi.so``, "" loads ``liblammps.so``)
  :type  name: string
  :param cmdargs: list of command line arguments to be passed to the :cpp:func:`lammps_open` function.  The executable name is automatically added.
  :type  cmdargs: list
  :param ptr: pointer to a LAMMPS C++ class instance when called from an embedded Python interpreter.  None means load symbols from shared library.
  :type  ptr: pointer
  :param comm: MPI communicator (as provided by `mpi4py <mpi4py_docs_>`_). ``None`` means use ``MPI_COMM_WORLD`` implicitly.
  :type  comm: MPI_Comm

  :ivar lmp:  instance of original LAMMPS Python interface
  :vartype lmp: :py:class:`lammps`

  :ivar runs:  list of completed runs, each storing the thermo output
  :vartype run: list
  """

  def __init__(self, name="", cmdargs=None, ptr=None, comm=None):
    self.has_echo = False

    if cmdargs:
      if '-echo' in cmdargs:
        idx = cmdargs.index('-echo')
        # ensures that echo line is ignored during output capture
        self.has_echo = idx+1 < len(cmdargs) and cmdargs[idx+1] in ('screen', 'both')

    if ptr:
      if isinstance(ptr,PyLammps):
        self.lmp = ptr.lmp
      elif isinstance(ptr,lammps):
        self.lmp = ptr
      else:
        self.lmp = lammps(name=name,cmdargs=cmdargs,ptr=ptr,comm=comm)
    else:
      self.lmp = lammps(name=name,cmdargs=cmdargs,ptr=None,comm=comm)
    print("LAMMPS output is captured by PyLammps wrapper")
    self._cmd_history = []
    self._enable_cmd_history = False
    self.runs = []

  def __enter__(self):
    return self

  def __exit__(self, ex_type, ex_value, ex_traceback):
    self.close()

  def __del__(self):
    if self.lmp: self.lmp.close()
    self.lmp = None

  def close(self):
    """Explicitly delete a LAMMPS instance

    This is a wrapper around the :py:meth:`lammps.close` of the Python interface.
    """
    if self.lmp: self.lmp.close()
    self.lmp = None

  def version(self):
    """Return a numerical representation of the LAMMPS version in use.

    This is a wrapper around the :py:meth:`lammps.version` function of the Python interface.

    :return: version number
    :rtype:  int
    """
    return self.lmp.version()

  def file(self, file):
    """Read LAMMPS commands from a file.

    This is a wrapper around the :py:meth:`lammps.file` function of the Python interface.

    :param path: Name of the file/path with LAMMPS commands
    :type path:  string
    """
    self.lmp.file(file)

  @property
  def enable_cmd_history(self):
    """
    :getter: Return whether command history is saved
    :setter: Set if command history should be saved
    :type: bool
    """
    return self._enable_cmd_history

  @enable_cmd_history.setter
  def enable_cmd_history(self, value):
    """
    :getter: Return whether command history is saved
    :setter: Set if command history should be saved
    :type: bool
    """
    self._enable_cmd_history = (value == True)

  def write_script(self, filepath):
    """
    Write LAMMPS script file containing all commands executed up until now

    :param filepath: path to script file that should be written
    :type filepath: string
    """
    with open(filepath, "w") as f:
      for cmd in self._cmd_history:
        print(cmd, file=f)

  def clear_cmd_history(self):
    """
    Clear LAMMPS command history up to this point
    """
    self._cmd_history = []

  def command(self, cmd):
    """
    Execute LAMMPS command

    If :py:attr:`PyLammps.enable_cmd_history` is set to ``True``, commands executed
    will be recorded. The entire command history can be written to a file using
    :py:meth:`PyLammps.write_script()`. To clear the command history, use
    :py:meth:`PyLammps.clear_cmd_history()`.

    :param cmd: command string that should be executed
    :type: cmd: string
    """
    self.lmp.command(cmd)

    if self.enable_cmd_history:
      self._cmd_history.append(cmd)

  def run(self, *args, **kwargs):
    """
    Execute LAMMPS run command with given arguments

    All thermo output during the run is captured and saved as new entry in
    :py:attr:`PyLammps.runs`. The latest run can be retrieved by
    :py:attr:`PyLammps.last_run`.
    """
    output = self.__getattr__('run')(*args, **kwargs)

    comm = self.lmp.get_mpi_comm()
    if comm:
      output = self.lmp.comm.bcast(output, root=0)

    self.runs += get_thermo_data(output)
    return output

  @property
  def last_run(self):
    """
    Return data produced of last completed run command

    :getter: Returns an object containing information about the last run command
    :type: dict
    """
    if len(self.runs) > 0:
        return self.runs[-1]
    return None

  @property
  def atoms(self):
    """
    All atoms of this LAMMPS instance

    :getter: Returns a list of atoms currently in the system
    :type: AtomList
    """
    return AtomList(self)

  @property
  def system(self):
    """
    The system state of this LAMMPS instance

    :getter: Returns an object with properties storing the current system state
    :type: namedtuple
    """
    output = self.info("system")
    d = self._parse_info_system(output)
    return namedtuple('System', d.keys())(*d.values())

  @property
  def communication(self):
    """
    The communication state of this LAMMPS instance

    :getter: Returns an object with properties storing the current communication state
    :type: namedtuple
    """
    output = self.info("communication")
    d = self._parse_info_communication(output)
    return namedtuple('Communication', d.keys())(*d.values())

  @property
  def computes(self):
    """
    The list of active computes of this LAMMPS instance

    :getter: Returns a list of computes that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.info("computes")
    return self._parse_element_list(output)

  @property
  def dumps(self):
    """
    The list of active dumps of this LAMMPS instance

    :getter: Returns a list of dumps that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.info("dumps")
    return self._parse_element_list(output)

  @property
  def fixes(self):
    """
    The list of active fixes of this LAMMPS instance

    :getter: Returns a list of fixes that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.info("fixes")
    return self._parse_element_list(output)

  @property
  def groups(self):
    """
    The list of active atom groups of this LAMMPS instance

    :getter: Returns a list of atom groups that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.info("groups")
    return self._parse_groups(output)

  @property
  def variables(self):
    """
    Returns a dictionary of all variables defined in the current LAMMPS instance

    :getter: Returns a dictionary of all variables that are defined in this LAMMPS instance
    :type: dict
    """
    output = self.info("variables")
    vars = {}
    for v in self._parse_element_list(output):
      vars[v['name']] = Variable(self, v['name'], v['style'], v['def'])
    return vars

  def eval(self, expr):
    """
    Evaluate expression

    :param expr: the expression string that should be evaluated inside of LAMMPS
    :type expr: string

    :return: the value of the evaluated expression
    :rtype: float if numeric, string otherwise
    """
    value = self.lmp_print('"$(%s)"' % expr).strip()
    try:
      return float(value)
    except ValueError:
      return value

  def _split_values(self, line):
    return [x.strip() for x in line.split(',')]

  def _get_pair(self, value):
    return [x.strip() for x in value.split('=')]

  def _parse_info_system(self, output):
    lines = output[6:-2]
    system = {}

    for line in lines:
      if line.startswith("Units"):
        system['units'] = self._get_pair(line)[1]
      elif line.startswith("Atom style"):
        system['atom_style'] = self._get_pair(line)[1]
      elif line.startswith("Atom map"):
        system['atom_map'] = self._get_pair(line)[1]
      elif line.startswith("Atoms"):
        parts = self._split_values(line)
        system['natoms'] = int(self._get_pair(parts[0])[1])
        system['ntypes'] = int(self._get_pair(parts[1])[1])
        system['style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Kspace style"):
        system['kspace_style'] = self._get_pair(line)[1]
      elif line.startswith("Dimensions"):
        system['dimensions'] = int(self._get_pair(line)[1])
      elif line.startswith("Orthogonal box"):
        system['orthogonal_box'] = [float(x) for x in self._get_pair(line)[1].split('x')]
      elif line.startswith("Boundaries"):
        system['boundaries'] = self._get_pair(line)[1]
      elif line.startswith("xlo"):
        keys, values = [self._split_values(x) for x in self._get_pair(line)]
        for key, value in zip(keys, values):
          system[key] = float(value)
      elif line.startswith("ylo"):
        keys, values = [self._split_values(x) for x in self._get_pair(line)]
        for key, value in zip(keys, values):
          system[key] = float(value)
      elif line.startswith("zlo"):
        keys, values = [self._split_values(x) for x in self._get_pair(line)]
        for key, value in zip(keys, values):
          system[key] = float(value)
      elif line.startswith("Molecule type"):
        system['molecule_type'] = self._get_pair(line)[1]
      elif line.startswith("Bonds"):
        parts = self._split_values(line)
        system['nbonds'] = int(self._get_pair(parts[0])[1])
        system['nbondtypes'] = int(self._get_pair(parts[1])[1])
        system['bond_style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Angles"):
        parts = self._split_values(line)
        system['nangles'] = int(self._get_pair(parts[0])[1])
        system['nangletypes'] = int(self._get_pair(parts[1])[1])
        system['angle_style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Dihedrals"):
        parts = self._split_values(line)
        system['ndihedrals'] = int(self._get_pair(parts[0])[1])
        system['ndihedraltypes'] = int(self._get_pair(parts[1])[1])
        system['dihedral_style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Impropers"):
        parts = self._split_values(line)
        system['nimpropers'] = int(self._get_pair(parts[0])[1])
        system['nimpropertypes'] = int(self._get_pair(parts[1])[1])
        system['improper_style'] = self._get_pair(parts[2])[1]

    return system

  def _parse_info_communication(self, output):
    lines = output[6:-3]
    comm = {}

    for line in lines:
      if line.startswith("MPI library"):
        comm['mpi_version'] = line.split(':')[1].strip()
      elif line.startswith("Comm style"):
        parts = self._split_values(line)
        comm['comm_style'] = self._get_pair(parts[0])[1]
        comm['comm_layout'] = self._get_pair(parts[1])[1]
      elif line.startswith("Processor grid"):
        comm['proc_grid'] = [int(x) for x in self._get_pair(line)[1].split('x')]
      elif line.startswith("Communicate velocities for ghost atoms"):
        comm['ghost_velocity'] = (self._get_pair(line)[1] == "yes")
      elif line.startswith("Nprocs"):
        parts = self._split_values(line)
        comm['nprocs'] = int(self._get_pair(parts[0])[1])
        comm['nthreads'] = int(self._get_pair(parts[1])[1])
    return comm

  def _parse_element_list(self, output):
    lines = output[6:-3]
    elements = []

    for line in lines:
      element_info = self._split_values(line.split(':')[1].strip())
      element = {'name': element_info[0]}
      for key, value in [self._get_pair(x) for x in element_info[1:]]:
        element[key] = value
      elements.append(element)
    return elements

  def _parse_groups(self, output):
    lines = output[6:-3]
    groups = []
    group_pattern = re.compile(r"(?P<name>.+) \((?P<type>.+)\)")

    for line in lines:
      m = group_pattern.match(line.split(':')[1].strip())
      group = {'name': m.group('name'), 'type': m.group('type')}
      groups.append(group)
    return groups

  def lmp_print(self, s):
    """ needed for Python2 compatibility, since print is a reserved keyword """
    return self.__getattr__("print")(s)

  def __dir__(self):
    return ['angle_coeff', 'angle_style', 'atom_modify', 'atom_style', 'atom_style',
    'bond_coeff', 'bond_style', 'boundary', 'change_box', 'communicate', 'compute',
    'create_atoms', 'create_box', 'delete_atoms', 'delete_bonds', 'dielectric',
    'dihedral_coeff', 'dihedral_style', 'dimension', 'dump', 'fix', 'fix_modify',
    'group', 'improper_coeff', 'improper_style', 'include', 'kspace_modify',
    'kspace_style', 'lattice', 'mass', 'minimize', 'min_style', 'neighbor',
    'neigh_modify', 'newton', 'nthreads', 'pair_coeff', 'pair_modify',
    'pair_style', 'processors', 'read', 'read_data', 'read_restart', 'region',
    'replicate', 'reset_timestep', 'restart', 'run', 'run_style', 'thermo',
    'thermo_modify', 'thermo_style', 'timestep', 'undump', 'unfix', 'units',
    'variable', 'velocity', 'write_restart']

  def __getattr__(self, name):
    """
    This method is where the Python 'magic' happens. If a method is not
    defined by the class PyLammps, it assumes it is a LAMMPS command. It takes
    all the arguments, concatinates them to a single string, and executes it using
    :py:meth:`lammps.PyLammps.command()`.

    :param verbose: Print output of command
    :type verbose:  bool
    :return: line or list of lines of output, None if no output
    :rtype: list or string
    """
    def handler(*args, **kwargs):
      cmd_args = [name] + [str(x) for x in args]

      with OutputCapture() as capture:
        cmd = ' '.join(cmd_args)
        self.command(cmd)
        output = capture.output

      if 'verbose' in kwargs and kwargs['verbose']:
        print(output)

      lines = output.splitlines()

      if self.has_echo:
        lines = lines[1:]

      if len(lines) > 1:
        return lines
      elif len(lines) == 1:
        return lines[0]
      return None

    return handler


class IPyLammps(PyLammps):
  """
  IPython wrapper for LAMMPS which adds embedded graphics capabilities to PyLammmps interface

  It either creates its own instance of :py:class:`lammps` or can be
  initialized with an existing instance. The arguments are the same of the
  lower-level interface. The original interface can still be accessed via
  :py:attr:`PyLammps.lmp`.

  :param name: "machine" name of the shared LAMMPS library ("mpi" loads ``liblammps_mpi.so``, "" loads ``liblammps.so``)
  :type  name: string
  :param cmdargs: list of command line arguments to be passed to the :cpp:func:`lammps_open` function.  The executable name is automatically added.
  :type  cmdargs: list
  :param ptr: pointer to a LAMMPS C++ class instance when called from an embedded Python interpreter.  None means load symbols from shared library.
  :type  ptr: pointer
  :param comm: MPI communicator (as provided by `mpi4py <mpi4py_docs_>`_). ``None`` means use ``MPI_COMM_WORLD`` implicitly.
  :type  comm: MPI_Comm
  """

  def __init__(self,name="",cmdargs=None,ptr=None,comm=None):
    super(IPyLammps, self).__init__(name=name,cmdargs=cmdargs,ptr=ptr,comm=comm)

  def image(self, filename="snapshot.png", group="all", color="type", diameter="type",
            size=None, view=None, center=None, up=None, zoom=1.0, background_color="white"):
    """ Generate image using write_dump command and display it

    See :doc:`dump image <dump_image>` for more information.

    :param filename: Name of the image file that should be generated. The extension determines whether it is PNG or JPEG
    :type filename: string
    :param group: the group of atoms write_image should use
    :type group: string
    :param color: name of property used to determine color
    :type color: string
    :param diameter: name of property used to determine atom diameter
    :type diameter: string
    :param size: dimensions of image
    :type size: tuple (width, height)
    :param view: view parameters
    :type view: tuple (theta, phi)
    :param center: center parameters
    :type center: tuple (flag, center_x, center_y, center_z)
    :param up: vector pointing to up direction
    :type up: tuple (up_x, up_y, up_z)
    :param zoom: zoom factor
    :type zoom: float
    :param background_color: background color of scene
    :type background_color: string

    :return: Image instance used to display image in notebook
    :rtype: :py:class:`IPython.core.display.Image`
    """
    cmd_args = [group, "image", filename, color, diameter]

    if size:
      width = size[0]
      height = size[1]
      cmd_args += ["size", width, height]

    if view:
      theta = view[0]
      phi = view[1]
      cmd_args += ["view", theta, phi]

    if center:
      flag = center[0]
      Cx = center[1]
      Cy = center[2]
      Cz = center[3]
      cmd_args += ["center", flag, Cx, Cy, Cz]

    if up:
      Ux = up[0]
      Uy = up[1]
      Uz = up[2]
      cmd_args += ["up", Ux, Uy, Uz]

    if zoom:
      cmd_args += ["zoom", zoom]

    cmd_args.append("modify backcolor " + background_color)

    self.write_dump(*cmd_args)
    from IPython.core.display import Image
    return Image(filename)

  def video(self, filename):
    """
    Load video from file

    Can be used to visualize videos from :doc:`dump movie <dump_image>`.

    :param filename: Path to video file
    :type filename: string
    :return: HTML Video Tag used by notebook to embed a video
    :rtype: :py:class:`IPython.display.HTML`
    """
    from IPython.display import HTML
    return HTML("<video controls><source src=\"" + filename + "\"></video>")
