# this only installs the LAMMPS python package
# it assumes the LAMMPS shared library is already installed
from distutils.core import setup
import os

LAMMPS_PYTHON_DIR = os.path.dirname(os.path.realpath(__file__))
LAMMPS_DIR = os.path.dirname(LAMMPS_PYTHON_DIR)
LAMMPS_SOURCE_DIR = os.path.join(LAMMPS_DIR, 'src')

def get_lammps_version():
    with open(os.path.join(LAMMPS_SOURCE_DIR, 'version.h'), 'r') as f:
        line = f.readline()
        start_pos = line.find('"')+1
        end_pos = line.find('"', start_pos)
        return "".join(line[start_pos:end_pos].split())

setup(
    name = "lammps",
    version = get_lammps_version(),
    author = "Steve Plimpton",
    author_email = "sjplimp@sandia.gov",
    url = "https://lammps.sandia.gov",
    description = "LAMMPS Molecular Dynamics Python package",
    license = "GPL",
    packages=["lammps","lammps.mliap"],
)
