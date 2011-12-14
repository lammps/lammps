#!/usr/local/bin/python

"""
setup.py file for LAMMPS with system MPICH library
"""

from distutils.core import setup, Extension

import os, glob
path = os.path.dirname(os.getcwd())

# list of src files for LAMMPS

libfiles = glob.glob("%s/src/*.cpp" % path)

lammps_library = Extension("_lammps",
                           sources = libfiles,
                           define_macros = [("MPICH_IGNORE_CXX_SEEK",1),
                                            ("LAMMPS_GZIP",1),
                                            ("FFT_NONE",1),],
                           # src files for LAMMPS
                           include_dirs = ["../src"],
                           # additional libs for MPICH on Linux
                           libraries = ["mpich","mpl","pthread"],
                           # where to find the MPICH lib on Linux
                           library_dirs = ["/usr/local/lib"],
                           # additional libs for MPI on Mac
                           # libraries = ["mpi"],
                           )

setup(name = "lammps",
      version = "28Nov11",
      author = "Steve Plimpton",
      author_email = "sjplimp@sandia.gov",
      url = "http://lammps.sandia.gov",
      description = """LAMMPS molecular dynamics library - parallel""",
      py_modules = ["lammps"],
      ext_modules = [lammps_library]
      )
