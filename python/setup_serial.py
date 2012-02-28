#!/usr/local/bin/python

"""
setup_serial.py file for LAMMPS with dummy serial MPI library
"""

from distutils.core import setup, Extension

import os, glob
path = os.path.dirname(os.getcwd())

# list of src files for LAMMPS and MPI STUBS

libfiles = glob.glob("%s/src/*.cpp" % path) + \
           glob.glob("%s/src/STUBS/*.c" % path)

lammps_library = Extension("_lammps_serial",
                           sources = libfiles,
                           define_macros = [("MPICH_IGNORE_CXX_SEEK",1),
                                            ("LAMMPS_GZIP",1),
                                            ("FFT_NONE",1),],
                           # src files for LAMMPS and MPI STUBS
                           include_dirs = ["../src", "../src/STUBS"]
                           )

setup(name = "lammps_serial",
      version = "28Nov11",
      author = "Steve Plimpton",
      author_email = "sjplimp@sandia.gov",
      url = "http://lammps.sandia.gov",
      description = """LAMMPS molecular dynamics library - serial""",
      py_modules = ["lammps"],
      ext_modules = [lammps_library]
      )
