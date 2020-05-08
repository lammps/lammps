
import sys
import unittest
from lammps import lammps
from io import StringIO
from contextlib import redirect_stdout

class PythonOpen(unittest.TestCase):

    def setUp(self):
        self.has_mpi=False
        try:
            from mpi4py import MPI
            self.has_mpi=True
        except:
            pass

    def testNoArgs(self):
        """Create LAMMPS instance without any arguments"""

        lmp=lammps()
        lmp.__del__()

    def testWithArgs(self):
        """Create LAMMPS instance with a few arguments"""
        lmp=lammps(cmdargs=['-nocite','-sf','opt','-log','none'])
        lmp.__del__()

    def testWithMPI(self):
        if self.has_mpi:
            from mpi4py import MPI
            mycomm=MPI.Comm.Split(MPI.COMM_WORLD, 0, 1)
            lmp=lammps(comm=mycomm)
            lmp.__del__()

if __name__ == "__main__":
    unittest.main()
