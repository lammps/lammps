
import sys,os,unittest
from lammps import lammps

has_mpi=False
has_mpi4py=False

try:
    from mpi4py import __version__ as mpi4py_version
    # tested to work with mpi4py versions 2 and 3
    has_mpi4py = mpi4py_version.split('.')[0] in ['2','3']
except:
    pass

try:
    if 'LAMMPS_MACHINE_NAME' in os.environ:
        machine = os.environ['LAMMPS_MACHINE_NAME']
    else:
        machine = ""
    lmp = lammps(name=machine)
    has_mpi = lmp.has_mpi_support
    lmp.close()
except:
    pass

class PythonOpen(unittest.TestCase):

    def setUp(self):
        self.machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            self.machine=os.environ['LAMMPS_MACHINE_NAME']

    def testNoArgs(self):
        """Create LAMMPS instance without any arguments"""

        lmp=lammps(name=self.machine)
        self.assertIsNot(lmp.lmp,None)
        self.assertEqual(lmp.opened,1)
        self.assertEqual(has_mpi and has_mpi4py,lmp.has_mpi4py)
        self.assertEqual(has_mpi,lmp.has_mpi_support)
        lmp.close()
        self.assertIsNone(lmp.lmp,None)
        self.assertEqual(lmp.opened,0)

    def testWithArgs(self):
        """Create LAMMPS instance with a few arguments"""
        lmp=lammps(name=self.machine,cmdargs=['-nocite','-sf','opt','-log','none'])
        self.assertIsNot(lmp.lmp,None)
        self.assertEqual(lmp.opened,1)

    def testError(self):
        """Print warning message through LAMMPS Error class"""
        lmp=lammps(name=self.machine,cmdargs=['-nocite','-log','none','-screen','tmp.error.output'])
        lmp.error(0,'test_warning')
        lmp.close()
        with open('tmp.error.output','r') as f:
            output = f.read()
        self.assertTrue('LAMMPS' in output)
        self.assertTrue('Total wall time' in output)
        self.assertTrue('WARNING: test_warning' in output)

    def testContextManager(self):
        """Automatically clean up LAMMPS instance"""
        with lammps(name=self.machine) as lmp:
            self.assertIsNot(lmp.lmp,None)
            self.assertEqual(lmp.opened,1)
            self.assertEqual(has_mpi and has_mpi4py,lmp.has_mpi4py)
            self.assertEqual(has_mpi,lmp.has_mpi_support)
        self.assertIsNone(lmp.lmp,None)
        self.assertEqual(lmp.opened,0)

    @unittest.skipIf(not (has_mpi and has_mpi4py),"Skipping MPI test since LAMMPS is not parallel or mpi4py is not found")
    def testWithMPI(self):
        from mpi4py import MPI
        mycomm=MPI.Comm.Split(MPI.COMM_WORLD, 0, 1)
        lmp=lammps(name=self.machine,comm=mycomm)
        self.assertIsNot(lmp.lmp,None)
        self.assertEqual(lmp.opened,1)
        lmp.close()

    def testUnknownCommand(self):
        lmp = lammps(name=self.machine)

        with self.assertRaisesRegex(Exception, "ERROR: Unknown command: write_paper"):
            lmp.command("write_paper")

        lmp.close()

    def testUnknownCommandInList(self):
        lmp = lammps(name=self.machine)

        with self.assertRaisesRegex(Exception, "ERROR: Unknown command: write_paper"):
            lmp.commands_list(["write_paper"])

        lmp.close()

    def testUnknownCommandInString(self):
        lmp = lammps(name=self.machine)

        with self.assertRaisesRegex(Exception, "ERROR: Unknown command: write_paper"):
            lmp.commands_string("write_paper")

        lmp.close()

if __name__ == "__main__":
    unittest.main()
