
import sys,os,unittest
from lammps import lammps

has_full=False
try:
    machine=None
    if 'LAMMPS_MACHINE_NAME' in os.environ:
        machine=os.environ['LAMMPS_MACHINE_NAME']
    lmp=lammps(name=machine)
    has_full = lmp.has_style("atom","full")
    lmp.close()
except:
    pass

class PythonGather(unittest.TestCase):

    def setUp(self):
        machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.lmp=lammps(name=machine,
                        cmdargs=['-nocite',
                                 '-log','none',
                                 '-echo','screen'])
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])

    # clean up temporary files
    def tearDown(self):
        self.lmp.close()

    # bond data comparison
    def checkBond(self, vals, btype, batom1, batom2):
        if ((vals[1] == batom1 and vals[2] == batom2)
            or (vals[1] == batom2 and vals[2] == batom1)):
            self.assertEqual(vals[0], btype)
            return 1
        else:
            return 0

    ##############################
    @unittest.skipIf(not has_full, "Gather_bonds test")
    def testGatherBond_newton_on(self):
        """Test gather_bonds() with newton on"""
        self.lmp.command("newton on on")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        nbonds, bonds = self.lmp.gather_bonds()
        self.assertEqual(nbonds, 24)
        self.assertEqual(len(bonds), 3*24)
        count = 0;
        for i in range(0,nbonds):
            count += self.checkBond(bonds[3*i:3*i+3], 5, 1, 2)
            count += self.checkBond(bonds[3*i:3*i+3], 3, 1, 3)
            count += self.checkBond(bonds[3*i:3*i+3], 2, 3, 4)
            count += self.checkBond(bonds[3*i:3*i+3], 2, 3, 5)
            count += self.checkBond(bonds[3*i:3*i+3], 1, 6, 3)
            count += self.checkBond(bonds[3*i:3*i+3], 3, 6, 8)
            count += self.checkBond(bonds[3*i:3*i+3], 4, 6, 7)
            count += self.checkBond(bonds[3*i:3*i+3], 5, 8, 9)
            count += self.checkBond(bonds[3*i:3*i+3], 5, 27, 28)
            count += self.checkBond(bonds[3*i:3*i+3], 5, 29, 27)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full, "Gather_bonds test")
    def testGatherBond_newton_off(self):
        """Test gather_bonds() with newton off"""
        self.lmp.command("newton off off")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        nbonds, bonds = self.lmp.gather_bonds()
        self.assertEqual(nbonds, 24)
        self.assertEqual(len(bonds), 3*24)
        count = 0;
        for i in range(0,nbonds):
            count += self.checkBond(bonds[3*i:3*i+3], 5, 1, 2)
            count += self.checkBond(bonds[3*i:3*i+3], 3, 1, 3)
            count += self.checkBond(bonds[3*i:3*i+3], 2, 3, 4)
            count += self.checkBond(bonds[3*i:3*i+3], 2, 3, 5)
            count += self.checkBond(bonds[3*i:3*i+3], 1, 3, 6)
            count += self.checkBond(bonds[3*i:3*i+3], 3, 6, 8)
            count += self.checkBond(bonds[3*i:3*i+3], 4, 6, 7)
            count += self.checkBond(bonds[3*i:3*i+3], 5, 8, 9)
            count += self.checkBond(bonds[3*i:3*i+3], 5, 27, 28)
            count += self.checkBond(bonds[3*i:3*i+3], 5, 27, 29)
        self.assertEqual(count,10)

##############################
if __name__ == "__main__":
    unittest.main()
