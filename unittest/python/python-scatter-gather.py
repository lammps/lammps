
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

    # angle data comparison
    def checkAngle(self, vals, atype, aatom1, aatom2, aatom3):
        set1 = {aatom1, aatom2, aatom3}
        set2 = {vals[1], vals[2], vals[3]}
        if len(set1.intersection(set2)) == 3:
            self.assertEqual(vals[0], atype)
            return 1
        else:
            return 0

    # dihedral data comparison
    def checkDihedral(self, vals, dtype, datom1, datom2, datom3, datom4):
        set1 = {datom1, datom2, datom3, datom4}
        set2 = {vals[1], vals[2], vals[3], vals[4]}
        if len(set1.intersection(set2)) == 4:
            self.assertEqual(vals[0], dtype)
            return 1
        else:
            return 0

    # improper data comparison
    def checkImproper(self, vals, itype, iatom1, iatom2, iatom3, iatom4):
        set1 = {iatom1, iatom2, iatom3, iatom4}
        set2 = {vals[1], vals[2], vals[3], vals[4]}
        if len(set1.intersection(set2)) == 4:
            self.assertEqual(vals[0], itype)
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

    @unittest.skipIf(not has_full, "Gather_angles test")
    def testGatherAngle_newton_on(self):
        """Test gather_angles() with newton on"""
        self.lmp.command("newton on on")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        nangles, angles = self.lmp.gather_angles()
        self.assertEqual(nangles, 30)
        self.assertEqual(len(angles), 120)
        count = 0;
        for i in range(0,nangles):
            count += self.checkAngle(angles[4*i:4*i+4], 4, 2, 1, 3)
            count += self.checkAngle(angles[4*i:4*i+4], 4, 1, 3, 6)
            count += self.checkAngle(angles[4*i:4*i+4], 2, 4, 3, 6)
            count += self.checkAngle(angles[4*i:4*i+4], 3, 7, 6, 8)
            count += self.checkAngle(angles[4*i:4*i+4], 3, 8, 10, 16)
            count += self.checkAngle(angles[4*i:4*i+4], 2, 8, 10, 11)
            count += self.checkAngle(angles[4*i:4*i+4], 1, 12, 10, 16)
            count += self.checkAngle(angles[4*i:4*i+4], 2, 10, 12, 14)
            count += self.checkAngle(angles[4*i:4*i+4], 1, 19, 18, 20)
            count += self.checkAngle(angles[4*i:4*i+4], 1, 28, 27, 29)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full, "Gather_angles test")
    def testGatherAngle_newton_off(self):
        """Test gather_angles() with newton off"""
        self.lmp.command("newton off off")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        nangles, angles = self.lmp.gather_angles()
        self.assertEqual(nangles, 30)
        self.assertEqual(len(angles), 120)
        count = 0;
        for i in range(0,nangles):
            count += self.checkAngle(angles[4*i:4*i+4], 4, 2, 1, 3)
            count += self.checkAngle(angles[4*i:4*i+4], 4, 1, 3, 6)
            count += self.checkAngle(angles[4*i:4*i+4], 2, 4, 3, 6)
            count += self.checkAngle(angles[4*i:4*i+4], 3, 7, 6, 8)
            count += self.checkAngle(angles[4*i:4*i+4], 3, 8, 10, 16)
            count += self.checkAngle(angles[4*i:4*i+4], 2, 8, 10, 11)
            count += self.checkAngle(angles[4*i:4*i+4], 1, 12, 10, 16)
            count += self.checkAngle(angles[4*i:4*i+4], 2, 10, 12, 14)
            count += self.checkAngle(angles[4*i:4*i+4], 1, 19, 18, 20)
            count += self.checkAngle(angles[4*i:4*i+4], 1, 28, 27, 29)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full, "Gather_dihedrals test")
    def testGatherDihedral_newton_on(self):
        """Test gather_dihedrals() with newton on"""
        self.lmp.command("newton on on")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        ndihedrals, dihedrals = self.lmp.gather_dihedrals()
        self.assertEqual(ndihedrals, 31)
        self.assertEqual(len(dihedrals), 155)
        count = 0;
        for i in range(0,ndihedrals):
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 2, 2, 1, 3, 6)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 3, 2, 1, 3, 5)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 1, 1, 3, 6, 7)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 5, 4, 3, 6, 7)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 4, 3, 6, 8, 9)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 2, 6, 8, 10, 16)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 5, 8, 10, 12, 13)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 1, 8, 10, 12, 14)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 2, 11, 10, 16, 17)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 5, 16, 10, 12, 14)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full, "Gather_dihedrals test")
    def testGatherDihedral_newton_off(self):
        """Test gather_dihedrals() with newton off"""
        self.lmp.command("newton off off")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        ndihedrals, dihedrals = self.lmp.gather_dihedrals()
        self.assertEqual(ndihedrals, 31)
        self.assertEqual(len(dihedrals), 155)
        count = 0;
        for i in range(0,ndihedrals):
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 2, 2, 1, 3, 6)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 3, 2, 1, 3, 5)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 1, 1, 3, 6, 7)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 5, 4, 3, 6, 7)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 4, 3, 6, 8, 9)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 2, 6, 8, 10, 16)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 5, 8, 10, 12, 13)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 1, 8, 10, 12, 14)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 2, 11, 10, 16, 17)
            count += self.checkDihedral(dihedrals[5*i:5*i+5], 5, 16, 10, 12, 14)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full, "Gather_impropers test")
    def testGatherImproper_newton_on(self):
        """Test gather_impropers() with newton on"""
        self.lmp.command("newton on on")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        nimpropers, impropers = self.lmp.gather_impropers()
        self.assertEqual(nimpropers, 2)
        self.assertEqual(len(impropers), 10)
        count = 0;
        for i in range(0,nimpropers):
            count += self.checkImproper(impropers[5*i:5*i+5], 1, 6, 3, 8, 7)
            count += self.checkImproper(impropers[5*i:5*i+5], 2, 8, 6, 10, 9)
        self.assertEqual(count,2)

    @unittest.skipIf(not has_full, "Gather_impropers test")
    def testGatherImproper_newton_off(self):
        """Test gather_impropers() with newton off"""
        self.lmp.command("newton off off")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        nimpropers, impropers = self.lmp.gather_impropers()
        self.assertEqual(nimpropers, 2)
        self.assertEqual(len(impropers), 10)
        count = 0;
        for i in range(0,nimpropers):
            count += self.checkImproper(impropers[5*i:5*i+5], 1, 6, 3, 8, 7)
            count += self.checkImproper(impropers[5*i:5*i+5], 2, 8, 6, 10, 9)
        self.assertEqual(count,2)

##############################
if __name__ == "__main__":
    unittest.main()
