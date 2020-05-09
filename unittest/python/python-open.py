
import sys,os,unittest
from lammps import lammps

class PythonOpen(unittest.TestCase):

    def setUp(self):
        self.has_mpi=False
        try:
            from mpi4py import MPI
            self.has_mpi=True
        except:
            pass
        self.machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            self.machine=os.environ['LAMMPS_MACHINE_NAME']

    def testNoArgs(self):
        """Create LAMMPS instance without any arguments"""

        lmp=lammps(name=self.machine)
        lmp.__del__()

    def testWithArgs(self):
        """Create LAMMPS instance with a few arguments"""
        lmp=lammps(name=self.machine,
                   cmdargs=['-nocite','-sf','opt','-log','none'])
        lmp.close()
        lmp.__del__()

    def testWithMPI(self):
        if self.has_mpi:
            from mpi4py import MPI
            mycomm=MPI.Comm.Split(MPI.COMM_WORLD, 0, 1)
            lmp=lammps(name=self.machine,comm=mycomm)
            lmp.__del__()

if __name__ == "__main__":
    unittest.main()
