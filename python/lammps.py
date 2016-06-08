# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   http://lammps.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

# Python wrapper on LAMMPS library via ctypes

import sys, traceback, types
from ctypes import *
from os.path import dirname, abspath, join
from inspect import getsourcefile
from collections import namedtuple
import os
import select
import re


class lammps:
    # detect if Python is using version of mpi4py that can pass a communicator

    has_mpi4py_v2 = False
    try:
        from mpi4py import MPI
        from mpi4py import __version__ as mpi4py_version
        if mpi4py_version.split('.')[0] == '2':
            has_mpi4py_v2 = True
    except:
        pass

    # create instance of LAMMPS

    def __init__(self, name="", cmdargs=None, ptr=None, comm=None):

        # determine module location

        modpath = dirname(abspath(getsourcefile(lambda: 0)))

        # load liblammps.so unless name is given.
        # e.g. if name = "g++", load liblammps_g++.so
        # try loading the LAMMPS shared object from the location
        # of lammps.py with an absolute path (so that LD_LIBRARY_PATH
        # does not need to be set for regular installations.
        # fall back to loading with a relative path, which typically
        # requires LD_LIBRARY_PATH to be set appropriately.

        try:
            if not name:
                self.lib = CDLL(join(modpath, "liblammps.so"), RTLD_GLOBAL)
            else:
                self.lib = CDLL(join(modpath, "liblammps_%s.so" % name), RTLD_GLOBAL)
        except:
            if not name:
                self.lib = CDLL("liblammps.so", RTLD_GLOBAL)
            else:
                self.lib = CDLL("liblammps_%s.so" % name, RTLD_GLOBAL)

        # if no ptr provided, create an instance of LAMMPS
        #   don't know how to pass an MPI communicator from PyPar
        #   but we can pass an MPI communicator from mpi4py v2.0.0 and later
        #   no_mpi call lets LAMMPS use MPI_COMM_WORLD
        #   cargs = array of C strings from args
        # if ptr, then are embedding Python in LAMMPS input script
        #   ptr is the desired instance of LAMMPS
        #   just convert it to ctypes ptr and store in self.lmp

        if not ptr:
            # with mpi4py v2, can pass MPI communicator to LAMMPS
            # need to adjust for type of MPI communicator object
            # allow for int (like MPICH) or void* (like OpenMPI)

            if lammps.has_mpi4py_v2 and comm != None:
                if lammps.MPI._sizeof(lammps.MPI.Comm) == sizeof(c_int):
                    MPI_Comm = c_int
                else:
                    MPI_Comm = c_void_p

                narg = 0
                cargs = 0
                if cmdargs:
                    cmdargs.insert(0, "lammps.py")
                    narg = len(cmdargs)
                    for i in range(narg):
                        if type(cmdargs[i]) is str:
                            cmdargs[i] = cmdargs[i].encode()
                    cargs = (c_char_p * narg)(*cmdargs)
                    self.lib.lammps_open.argtypes = [c_int, c_char_p * narg, \
                                                     MPI_Comm, c_void_p()]
                else:
                    self.lib.lammps_open.argtypes = [c_int, c_int, \
                                                     MPI_Comm, c_void_p()]

                self.lib.lammps_open.restype = None
                self.opened = 1
                self.lmp = c_void_p()
                comm_ptr = lammps.MPI._addressof(comm)
                comm_val = MPI_Comm.from_address(comm_ptr)
                self.lib.lammps_open(narg, cargs, comm_val, byref(self.lmp))

            else:
                self.opened = 1
                if cmdargs:
                    cmdargs.insert(0, "lammps.py")
                    narg = len(cmdargs)
                    for i in range(narg):
                        if type(cmdargs[i]) is str:
                            cmdargs[i] = cmdargs[i].encode()
                    cargs = (c_char_p * narg)(*cmdargs)
                    self.lmp = c_void_p()
                    self.lib.lammps_open_no_mpi(narg, cargs, byref(self.lmp))
                else:
                    self.lmp = c_void_p()
                    self.lib.lammps_open_no_mpi(0, None, byref(self.lmp))
                    # could use just this if LAMMPS lib interface supported it
                    # self.lmp = self.lib.lammps_open_no_mpi(0,None)

        else:
            self.opened = 0
            # magic to convert ptr to ctypes ptr
            pythonapi.PyCObject_AsVoidPtr.restype = c_void_p
            pythonapi.PyCObject_AsVoidPtr.argtypes = [py_object]
            self.lmp = c_void_p(pythonapi.PyCObject_AsVoidPtr(ptr))

    def __del__(self):
        if self.lmp and self.opened: self.lib.lammps_close(self.lmp)

    def close(self):
        if self.opened: self.lib.lammps_close(self.lmp)
        self.lmp = None

    def version(self):
        return self.lib.lammps_version(self.lmp)

    def file(self, file):
        file = file.encode()
        self.lib.lammps_file(self.lmp, file)

    def command(self, cmd):
        cmd = cmd.encode()
        self.lib.lammps_command(self.lmp, cmd)

    def extract_global(self, name, type):
        name = name.encode()
        if type == 0:
            self.lib.lammps_extract_global.restype = POINTER(c_int)
        elif type == 1:
            self.lib.lammps_extract_global.restype = POINTER(c_double)
        else:
            return None
        ptr = self.lib.lammps_extract_global(self.lmp, name)
        return ptr[0]

    def extract_atom(self, name, type):
        name = name.encode()
        if type == 0:
            self.lib.lammps_extract_atom.restype = POINTER(c_int)
        elif type == 1:
            self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_int))
        elif type == 2:
            self.lib.lammps_extract_atom.restype = POINTER(c_double)
        elif type == 3:
            self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_double))
        else:
            return None
        ptr = self.lib.lammps_extract_atom(self.lmp, name)
        return ptr

    def extract_compute(self, id, style, type):
        id = id.encode()
        if type == 0:
            if style > 0: return None
            self.lib.lammps_extract_compute.restype = POINTER(c_double)
            ptr = self.lib.lammps_extract_compute(self.lmp, id, style, type)
            return ptr[0]
        if type == 1:
            self.lib.lammps_extract_compute.restype = POINTER(c_double)
            ptr = self.lib.lammps_extract_compute(self.lmp, id, style, type)
            return ptr
        if type == 2:
            self.lib.lammps_extract_compute.restype = POINTER(POINTER(c_double))
            ptr = self.lib.lammps_extract_compute(self.lmp, id, style, type)
            return ptr
        return None

    # in case of global datum, free memory for 1 double via lammps_free()
    # double was allocated by library interface function

    def extract_fix(self, id, style, type, i=0, j=0):
        id = id.encode()
        if style == 0:
            self.lib.lammps_extract_fix.restype = POINTER(c_double)
            ptr = self.lib.lammps_extract_fix(self.lmp, id, style, type, i, j)
            result = ptr[0]
            self.lib.lammps_free(ptr)
            return result
        elif (style == 1) or (style == 2):
            if type == 1:
                self.lib.lammps_extract_fix.restype = POINTER(c_double)
            elif type == 2:
                self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
            else:
                return None
            ptr = self.lib.lammps_extract_fix(self.lmp, id, style, type, i, j)
            return ptr
        else:
            return None

    # free memory for 1 double or 1 vector of doubles via lammps_free()
    # for vector, must copy nlocal returned values to local c_double vector
    # memory was allocated by library interface function

    def extract_variable(self, name, group, type):
        name = name.encode()
        group = group.encode()
        if type == 0:
            self.lib.lammps_extract_variable.restype = POINTER(c_double)
            ptr = self.lib.lammps_extract_variable(self.lmp, name, group)
            result = ptr[0]
            self.lib.lammps_free(ptr)
            return result
        if type == 1:
            self.lib.lammps_extract_global.restype = POINTER(c_int)
            nlocalptr = self.lib.lammps_extract_global(self.lmp, "nlocal")
            nlocal = nlocalptr[0]
            result = (c_double * nlocal)()
            self.lib.lammps_extract_variable.restype = POINTER(c_double)
            ptr = self.lib.lammps_extract_variable(self.lmp, name, group)
            for i in range(nlocal): result[i] = ptr[i]
            self.lib.lammps_free(ptr)
            return result
        return None

    # set variable value
    # value is converted to string
    # returns 0 for success, -1 if failed

    def set_variable(self, name, value):
        name = name.encode()
        value = str(value).encode()
        return self.lib.lammps_set_variable(self.lmp, name, str(value))

    # return total number of atoms in system

    def get_natoms(self):
        return self.lib.lammps_get_natoms(self.lmp)

    # return vector of atom properties gathered across procs, ordered by atom ID

    def gather_atoms(self, name, type, count):
        name = name.encode()
        natoms = self.lib.lammps_get_natoms(self.lmp)
        if type == 0:
            data = ((count * natoms) * c_int)()
            self.lib.lammps_gather_atoms(self.lmp, name, type, count, data)
        elif type == 1:
            data = ((count * natoms) * c_double)()
            self.lib.lammps_gather_atoms(self.lmp, name, type, count, data)
        else:
            return None
        return data

    # scatter vector of atom properties across procs, ordered by atom ID
    # assume vector is of correct type and length, as created by gather_atoms()

    def scatter_atoms(self, name, type, count, data):
        name = name.encode()
        self.lib.lammps_scatter_atoms(self.lmp, name, type, count, data)


class OutputCapture(object):
    def __init__(self):
        self.stdout_pipe_read, self.stdout_pipe_write = os.pipe()
        self.stdout_fd = 1

    def __enter__(self):
        self.stdout = os.dup(self.stdout_fd)
        os.dup2(self.stdout_pipe_write, self.stdout_fd)
        return self

    def __exit__(self, type, value, tracebac):
        os.dup2(self.stdout, self.stdout_fd)

    # check if we have more to read from the pipe
    def more_data(self, pipe):
        r, _, _ = select.select([pipe], [], [], 0)
        return bool(r)

    # read the whole pipe
    def read_pipe(self, pipe):
        out = ""
        while self.more_data(pipe):
            out += str(os.read(pipe, 1024), 'utf-8')
        return out

    @property
    def output(self):
        return self.read_pipe(self.stdout_pipe_read)


class Variable(object):
    def __init__(self, lammps_wrapper_instance, name, style, definition):
        self.lmp = lammps_wrapper_instance
        self.name = name
        self.style = style
        self.definition = definition.split()

    @property
    def value(self):
        value = self.lmp.print('"${%s}"' % self.name).strip()
        try:
            return float(value)
        except ValueError:
            return value


class AtomList(object):
    def __init__(self, lammps_wrapper_instance):
        self.lmp = lammps_wrapper_instance
        self.natoms = self.lmp.system.natoms

    def __getitem__(self, index):
        return Atom(self.lmp, index+1)


class Atom(object):
    def __init__(self, lammps_wrapper_instance, index):
        self.lmp = lammps_wrapper_instance
        self.index = index

    @property
    def id(self):
        return int(self.lmp.eval("id[%d]" % self.index))

    @property
    def type(self):
        return int(self.lmp.eval("type[%d]" % self.index))

    @property
    def mol(self):
        return self.lmp.eval("mol[%d]" % self.index)

    @property
    def mass(self):
        return self.lmp.eval("mass[%d]" % self.index)

    @property
    def position(self):
        return (self.lmp.eval("x[%d]" % self.index),
                self.lmp.eval("y[%d]" % self.index),
                self.lmp.eval("z[%d]" % self.index))

    @property
    def velocity(self):
        return (self.lmp.eval("vx[%d]" % self.index),
                self.lmp.eval("vy[%d]" % self.index),
                self.lmp.eval("vz[%d]" % self.index))

    @property
    def force(self):
        return (self.lmp.eval("fx[%d]" % self.index),
                self.lmp.eval("fy[%d]" % self.index),
                self.lmp.eval("fz[%d]" % self.index))

    @property
    def charge(self):
        return self.lmp.eval("q[%d]" % self.index)


class LammpsWrapper(object):
    def __init__(self, lmp):
        self.lmp = lmp

    @property
    def atoms(self):
        return AtomList(self)

    @property
    def system(self):
        output = self.info("system")
        d = self._parse_info_system(output)
        return namedtuple('System', d.keys())(*d.values())

    @property
    def communication(self):
        output = self.info("communication")
        d = self._parse_info_communication(output)
        return namedtuple('Communication', d.keys())(*d.values())

    @property
    def computes(self):
        output = self.info("computes")
        return self._parse_element_list(output)

    @property
    def dumps(self):
        output = self.info("dumps")
        return self._parse_element_list(output)

    @property
    def fixes(self):
        output = self.info("fixes")
        return self._parse_element_list(output)

    @property
    def groups(self):
        output = self.info("groups")
        return self._parse_groups(output)

    @property
    def variables(self):
        output = self.info("variables")
        vars = {}
        for v in self._parse_element_list(output):
            vars[v['name']] = Variable(self, v['name'], v['style'], v['def'])
        return vars

    def eval(self, expr):
        value = self.print('"$(%s)"' % expr).strip()
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

    def __getattr__(self, name):
        def handler(*args, **kwargs):
            cmd_args = [name] + [str(x) for x in args]

            with OutputCapture() as capture:
                self.lmp.command(' '.join(cmd_args))
                output = capture.output
            lines = output.splitlines()

            if len(lines) > 1:
                return lines
            elif len(lines) == 1:
                return lines[0]
            return None

        return handler


class LammpsIPythonWrapper(LammpsWrapper):
    def __init__(self, lmp):
        super(LammpsIPythonWrapper, self).__init__(lmp)

    def image(self, filename="snapshot.png", group="all", color="type", diameter="type",
              size=None, view=None, center=None, up=None, zoom=1.0):
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

        cmd_args.append("modify backcolor white")

        self.write_dump(*cmd_args)
        from IPython.core.display import Image
        return Image('snapshot.png')
