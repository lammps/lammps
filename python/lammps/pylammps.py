# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   LAMMPS Development team: developers@lammps.org
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
# Written by Richard Berger <richard.berger@outlook.com>
################################################################################

# for python2/3 compatibility

from __future__ import print_function

import io
import os
import re
import sys
import tempfile
from collections import namedtuple

from .core import lammps
from .constants import *                # lgtm [py/polluting-import]

# -------------------------------------------------------------------------

class OutputCapture(object):
  """ Utility class to capture LAMMPS library output """
  def __init__(self):
    self.stdout_fd = 1
    self.captured_output = ""

  def __enter__(self):
    self.tmpfile = tempfile.TemporaryFile(mode='w+b')

    sys.stdout.flush()

    # make copy of original stdout
    self.stdout_orig = os.dup(self.stdout_fd)

    # replace stdout and redirect to temp file
    os.dup2(self.tmpfile.fileno(), self.stdout_fd)
    return self

  def __exit__(self, exc_type, exc_value, traceback):
    os.dup2(self.stdout_orig, self.stdout_fd)
    os.close(self.stdout_orig)
    self.tmpfile.close()

  @property
  def output(self):
    sys.stdout.flush()
    self.tmpfile.flush()
    self.tmpfile.seek(0, io.SEEK_SET)
    self.captured_output = self.tmpfile.read().decode('utf-8')
    return self.captured_output

# -------------------------------------------------------------------------

class Variable(object):
  def __init__(self, pylammps_instance, name):
    self._pylmp = pylammps_instance
    self.name = name

  @property
  def style(self):
    vartype = self._pylmp.lmp.lib.lammps_extract_variable_datatype(self._pylmp.lmp.lmp, self.name.encode())
    if vartype == LMP_VAR_EQUAL:
      return "equal"
    elif vartype == LMP_VAR_ATOM:
      return "atom"
    elif vartype == LMP_VAR_VECTOR:
      return "vector"
    elif vartype == LMP_VAR_STRING:
      return "string"
    return None

  @property
  def value(self):
    return self._pylmp.lmp.extract_variable(self.name)

  @value.setter
  def value(self, newvalue):
    style = self.style
    if style == "equal" or style == "string":
      self._pylmp.variable("{} {} {}".format(self.name, style, newvalue))
    else:
      raise Exception("Setter not implemented for {} style variables.".format(style))

  def __str__(self):
    value = self.value
    if isinstance(value, str):
      value = "\"{}\"".format(value)
    return "Variable(name=\"{}\", value={})".format(self.name, value)

  def __repr__(self):
    return self.__str__()

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
            atom = Atom2D(self._pylmp, index)
        else:
            atom = Atom(self._pylmp, index)
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

  def get(self, name, index):
    prop = self._pylmp.lmp.numpy.extract_atom(name)
    if prop is not None:
      return prop[index]
    return None

  @property
  def id(self):
    """
    Return the atom ID

    :type: int
    """
    return self.get("id", self.index)

  @property
  def type(self):
    """
    Return the atom type

    :type: int
    """
    return self.get("type", self.index)

  @property
  def mol(self):
    """
    Return the atom molecule index

    :type: int
    """
    return self.get("mol", self.index)

  @property
  def mass(self):
    """
    Return the atom mass

    :type: float
    """
    return self.get("mass", self.index)

  @property
  def radius(self):
    """
    Return the particle radius

    :type: float
    """
    return self.get("radius", self.index)

  @property
  def position(self):
    """
    :getter: Return position of atom
    :setter: Set position of atom
    :type: numpy.array (float, float, float)
    """
    return self.get("x", self.index)

  @position.setter
  def position(self, value):
    current = self.position
    current[:] = value

  @property
  def velocity(self):
    """
    :getter: Return velocity of atom
    :setter: Set velocity of atom
    :type: numpy.array (float, float, float)
    """
    return self.get("v", self.index)

  @velocity.setter
  def velocity(self, value):
    current = self.velocity
    current[:] = value

  @property
  def force(self):
    """
    Return the total force acting on the atom

    :type: numpy.array (float, float, float)
    """
    return self.get("f", self.index)

  @force.setter
  def force(self, value):
    current = self.force
    current[:] = value

  @property
  def torque(self):
    """
    Return the total torque acting on the atom

    :type: numpy.array (float, float, float)
    """
    return self.get("torque", self.index)

  @force.setter
  def torque(self, value):
    current = self.torque
    current[:] = value

  @property
  def omega(self):
    """
    Return the rotational velocity of the particle

    :type: numpy.array (float, float, float)
    """
    return self.get("torque", self.index)

  @omega.setter
  def omega(self, value):
    current = self.torque
    current[:] = value

  @property
  def torque(self):
    """
    Return the total torque acting on the particle

    :type: numpy.array (float, float, float)
    """
    return self.get("torque", self.index)

  @torque.setter
  def torque(self, value):
    current = self.torque
    current[:] = value

  @property
  def angular_momentum(self):
    """
    Return the angular momentum of the particle

    :type: numpy.array (float, float, float)
    """
    return self.get("angmom", self.index)

  @angular_momentum.setter
  def angular_momentum(self, value):
    current = self.angular_momentum
    current[:] = value

  @property
  def charge(self):
    """
    Return the atom charge

    :type: float
    """
    return self.get("q", self.index)

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
    :type: numpy.array (float, float)
    """
    return super(Atom2D, self).position[0:2]

  @position.setter
  def position(self, value):
    current = self.position
    current[:] = value

  @property
  def velocity(self):
    """Access to velocity of an atom
    :getter: Return velocity of atom
    :setter: Set velocity of atom
    :type: numpy.array (float, float)
    """
    return super(Atom2D, self).velocity[0:2]

  @velocity.setter
  def velocity(self, value):
    current = self.velocity
    current[:] = value

  @property
  def force(self):
    """Access to force of an atom
    :getter: Return force of atom
    :setter: Set force of atom
    :type: numpy.array (float, float)
    """
    return super(Atom2D, self).force[0:2]

  @force.setter
  def force(self, value):
    current = self.force
    current[:] = value

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
  :param verbose: print all LAMMPS output to stdout
  :type  verbose: bool

  :ivar lmp:  instance of original LAMMPS Python interface
  :vartype lmp: :py:class:`lammps`

  :ivar runs:  list of completed runs, each storing the thermo output
  :vartype run: list
  """

  def __init__(self, name="", cmdargs=None, ptr=None, comm=None, verbose=False):
    self.has_echo = False
    self.verbose = verbose

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

    if not self.lmp.has_package("PYTHON"):
      print("WARNING: run thermo data not captured since PYTHON LAMMPS package is not enabled")

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

  def _append_run_thermo(self, thermo):
    for k, v in thermo.items():
      if k in self._current_run:
        self._current_run[k].append(v)
      else:
        self._current_run[k] = [v]

  def run(self, *args, **kwargs):
    """
    Execute LAMMPS run command with given arguments

    Thermo data of the run is recorded and saved as new entry in
    :py:attr:`PyLammps.runs`. The latest run can be retrieved by
    :py:attr:`PyLammps.last_run`.

    Note, for recording of all thermo steps during a run, the PYTHON package
    needs to be enabled in LAMMPS. Otherwise, it will only capture the final
    timestep.
    """
    self._current_run = {}
    self._last_thermo_step = -1
    def end_of_step_callback(lmp):
      if self.lmp.last_thermo_step == self._last_thermo_step: return
      thermo = self.lmp.last_thermo()
      self._append_run_thermo(thermo)
      self._last_thermo_step = thermo['Step']

    import __main__
    __main__._PyLammps_end_of_step_callback = end_of_step_callback
    capture_thermo = self.lmp.has_package("PYTHON")

    if capture_thermo:
        self.fix("__pylammps_internal_run_callback", "all", "python/invoke", "1", "end_of_step", "_PyLammps_end_of_step_callback")

    output = self.__getattr__('run')(*args, **kwargs)

    if capture_thermo:
        self.unfix("__pylammps_internal_run_callback")
    self._append_run_thermo(self.lmp.last_thermo())

    thermo_data = variable_set('ThermoData', self._current_run)
    r = {'thermo' : thermo_data }
    self.runs.append(namedtuple('Run', list(r.keys()))(*list(r.values())))
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
    output = self.lmp_info("system")
    output = output[output.index("System information:")+1:]
    d = self._parse_info_system(output)
    return namedtuple('System', d.keys())(*d.values())

  @property
  def communication(self):
    """
    The communication state of this LAMMPS instance

    :getter: Returns an object with properties storing the current communication state
    :type: namedtuple
    """
    output = self.lmp_info("communication")
    output = output[output.index("Communication information:")+1:]
    d = self._parse_info_communication(output)
    return namedtuple('Communication', d.keys())(*d.values())

  @property
  def computes(self):
    """
    The list of active computes of this LAMMPS instance

    :getter: Returns a list of computes that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.lmp_info("computes")
    output = output[output.index("Compute information:")+1:]
    return self._parse_element_list(output)

  @property
  def dumps(self):
    """
    The list of active dumps of this LAMMPS instance

    :getter: Returns a list of dumps that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.lmp_info("dumps")
    output = output[output.index("Dump information:")+1:]
    return self._parse_element_list(output)

  @property
  def fixes(self):
    """
    The list of active fixes of this LAMMPS instance

    :getter: Returns a list of fixes that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.lmp_info("fixes")
    output = output[output.index("Fix information:")+1:]
    return self._parse_element_list(output)

  @property
  def groups(self):
    """
    The list of active atom groups of this LAMMPS instance

    :getter: Returns a list of atom groups that are currently active in this LAMMPS instance
    :type: list
    """
    return self.lmp.available_ids("group")

  @property
  def variables(self):
    """
    Returns a dictionary of all variables defined in the current LAMMPS instance

    :getter: Returns a dictionary of all variables that are defined in this LAMMPS instance
    :type: dict
    """
    variables = {}
    for name in self.lmp.available_ids("variable"):
      variables[name] = Variable(self, name)
    return variables

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
    system = {}
    system['dimensions'] = self.lmp.extract_setting("dimension")
    system['xlo'] = self.lmp.extract_global("boxxlo")
    system['ylo'] = self.lmp.extract_global("boxylo")
    system['zlo'] = self.lmp.extract_global("boxzlo")
    system['xhi'] = self.lmp.extract_global("boxxhi")
    system['yhi'] = self.lmp.extract_global("boxyhi")
    system['zhi'] = self.lmp.extract_global("boxzhi")
    xprd = system["xhi"] - system["xlo"]
    yprd = system["yhi"] - system["ylo"]
    zprd = system["zhi"] - system["zlo"]
    if self.lmp.extract_setting("triclinic") == 1:
      system['triclinic_box'] = (xprd, yprd, zprd)
    else:
      system['orthogonal_box'] = (xprd, yprd, zprd)
    system['nangles'] = self.lmp.extract_global("nbonds")
    system['nangletypes'] = self.lmp.extract_setting("nbondtypes")
    system['angle_style'] = self.lmp.extract_global("angle_style")
    system['nbonds'] = self.lmp.extract_global("nbonds")
    system['nbondtypes'] = self.lmp.extract_setting("nbondtypes")
    system['bond_style'] = self.lmp.extract_global("bond_style")
    system['ndihedrals'] = self.lmp.extract_global("ndihedrals")
    system['ndihedraltypes'] = self.lmp.extract_setting("ndihedraltypes")
    system['dihedral_style'] = self.lmp.extract_global("dihedral_style")
    system['nimpropers'] = self.lmp.extract_global("nimpropers")
    system['nimpropertypes'] = self.lmp.extract_setting("nimpropertypes")
    system['improper_style'] = self.lmp.extract_global("improper_style")
    system['kspace_style'] = self.lmp.extract_global("kspace_style")
    system['natoms'] = self.lmp.extract_global("natoms")
    system['ntypes'] = self.lmp.extract_global("ntypes")
    system['pair_style'] = self.lmp.extract_global("pair_style")
    system['atom_style'] = self.lmp.extract_global("atom_style")
    system['units'] = self.lmp.extract_global("units")

    for line in output:
      if line.startswith("Atom map"):
        system['atom_map'] = self._get_pair(line)[1]
      elif line.startswith("Boundaries"):
        system['boundaries'] = self._get_pair(line)[1]
      elif line.startswith("Molecule type"):
        system['molecule_type'] = self._get_pair(line)[1]

    return system

  def _parse_info_communication(self, output):
    comm = {}
    comm['nprocs'] = self.lmp.extract_setting("world_size")
    comm['nthreads'] = self.lmp.extract_setting("nthreads")

    for line in output:
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
    return comm

  def _parse_element_list(self, output):
    elements = []

    for line in output:
      if not line or (":" not in line): continue
      element_info = self._split_values(line.split(':')[1].strip())
      element = {'name': element_info[0]}
      for key, value in [self._get_pair(x) for x in element_info[1:]]:
        element[key] = value
      elements.append(element)
    return elements

  def lmp_print(self, s):
    """ needed for Python2 compatibility, since print is a reserved keyword """
    return self.__getattr__("print")(s)

  def __dir__(self):
    return sorted(set(['angle_coeff', 'angle_style', 'atom_modify', 'atom_style', 'atom_style',
    'bond_coeff', 'bond_style', 'boundary', 'change_box', 'communicate', 'compute',
    'create_atoms', 'create_box', 'delete_atoms', 'delete_bonds', 'dielectric',
    'dihedral_coeff', 'dihedral_style', 'dimension', 'dump', 'fix', 'fix_modify',
    'group', 'improper_coeff', 'improper_style', 'include', 'kspace_modify',
    'kspace_style', 'lattice', 'mass', 'minimize', 'min_style', 'neighbor',
    'neigh_modify', 'newton', 'nthreads', 'pair_coeff', 'pair_modify',
    'pair_style', 'processors', 'read', 'read_data', 'read_restart', 'region',
    'replicate', 'reset_timestep', 'restart', 'run', 'run_style', 'thermo',
    'thermo_modify', 'thermo_style', 'timestep', 'undump', 'unfix', 'units',
    'variable', 'velocity', 'write_restart'] + self.lmp.available_styles("command")))

  def lmp_info(self, s):
      # skip anything before and after Info-Info-Info
      # also skip timestamp line
      output = self.__getattr__("info")(s)
      indices = [index for index, line in enumerate(output) if line.startswith("Info-Info-Info-Info")]
      start = indices[0]
      end = indices[1]
      return [line for line in output[start+2:end] if line]

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
      self.lmp.flush_buffers()

      with OutputCapture() as capture:
        cmd = ' '.join(cmd_args)
        self.command(cmd)
        self.lmp.flush_buffers()
        output = capture.output

      comm = self.lmp.get_mpi_comm()
      if comm:
        output = self.lmp.comm.bcast(output, root=0)

      if self.verbose or ('verbose' in kwargs and kwargs['verbose']):
        print(output, end = '')

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

    if size is not None:
      width = size[0]
      height = size[1]
      cmd_args += ["size", width, height]

    if view is not None:
      theta = view[0]
      phi = view[1]
      cmd_args += ["view", theta, phi]

    if center is not None:
      flag = center[0]
      Cx = center[1]
      Cy = center[2]
      Cz = center[3]
      cmd_args += ["center", flag, Cx, Cy, Cz]

    if up is not None:
      Ux = up[0]
      Uy = up[1]
      Uz = up[2]
      cmd_args += ["up", Ux, Uy, Uz]

    if zoom is not None:
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
