import sys,os,unittest
from lammps import lammps, LAMMPS_INT, LMP_STYLE_GLOBAL, LMP_STYLE_LOCAL, \
                   LMP_STYLE_ATOM, LMP_TYPE_VECTOR, LMP_TYPE_SCALAR, LMP_TYPE_ARRAY, \
                   LMP_VAR_ATOM
from ctypes import c_void_p

has_manybody=False
try:
    machine=None
    if 'LAMMPS_MACHINE_NAME' in os.environ:
        machine=os.environ['LAMMPS_MACHINE_NAME']
    lmp=lammps(name=machine)
    has_manybody = lmp.has_style("pair","sw")
    lmp.close()
except:
    pass

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

try:
    import numpy
    NUMPY_INSTALLED = True
except ImportError:
    NUMPY_INSTALLED = False

@unittest.skipIf(not NUMPY_INSTALLED, "numpy is not available")
class PythonNumpy(unittest.TestCase):
    def setUp(self):
        machine = None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.lmp = lammps(name=machine,  cmdargs=['-nocite', '-log','none', '-echo','screen'])

    def tearDown(self):
        del self.lmp

    def checkBond(self, vals, btype, batom1, batom2):
        if ((vals[1] == batom1 and vals[2] == batom2)
            or (vals[1] == batom2 and vals[2] == batom1)):
            self.assertEqual(vals[0], btype)
            return 1
        else:
            return 0

    def checkAngle(self, vals, atype, aatom1, aatom2, aatom3):
        set1 = {aatom1, aatom2, aatom3}
        set2 = {vals[1], vals[2], vals[3]}
        if len(set1.intersection(set2)) == 3:
            self.assertEqual(vals[0], atype)
            return 1
        else:
            return 0

    def checkDihedral(self, vals, dtype, datom1, datom2, datom3, datom4):
        set1 = {datom1, datom2, datom3, datom4}
        set2 = {vals[1], vals[2], vals[3], vals[4]}
        if len(set1.intersection(set2)) == 4:
            self.assertEqual(vals[0], dtype)
            return 1
        else:
            return 0

    def checkImproper(self, vals, itype, iatom1, iatom2, iatom3, iatom4):
        set1 = {iatom1, iatom2, iatom3, iatom4}
        set2 = {vals[1], vals[2], vals[3], vals[4]}
        if len(set1.intersection(set2)) == 4:
            self.assertEqual(vals[0], itype)
            return 1
        else:
            return 0

    def testLammpsPointer(self):
        self.assertEqual(type(self.lmp.lmp), c_void_p)

    def testExtractComputeGlobalScalar(self):
        # TODO
        pass

    def testExtractComputeGlobalVector(self):
        self.lmp.command("region       box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")
        self.lmp.command("create_atoms 1 single 1.0 1.0 1.0")
        self.lmp.command("create_atoms 1 single 1.0 1.0 1.5")
        self.lmp.command("compute coordsum all reduce sum x y z")
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,2)
        values = self.lmp.numpy.extract_compute("coordsum", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR)
        self.assertEqual(len(values), 3)
        self.assertEqual(values[0], 2.0)
        self.assertEqual(values[1], 2.0)
        self.assertEqual(values[2], 2.5)

    def testExtractComputeGlobalArray(self):
        # TODO
        pass

    def testExtractComputePerAtomVector(self):
        self.lmp.command("region       box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")
        self.lmp.command("create_atoms 1 single 1.0 1.0 1.0")
        self.lmp.command("create_atoms 1 single 1.0 1.0 1.5")
        self.lmp.command("mass 1 1.0")
        self.lmp.command("compute ke all ke/atom")
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,2)
        values = self.lmp.numpy.extract_compute("ke", LMP_STYLE_ATOM, LMP_TYPE_VECTOR)
        self.assertEqual(len(values), 2)
        self.assertEqual(values[0], 0.0)
        self.assertEqual(values[1], 0.0)

    def testExtractComputePerAtomArray(self):
        # TODO
        pass

    def testExtractComputeLocalVector(self):
        self.lmp.command("region       box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")
        self.lmp.command("create_atoms 1 single 1.0 1.0 1.0")
        self.lmp.command("create_atoms 1 single 1.0 1.0 1.5")
        self.lmp.command("mass 1 1.0")
        self.lmp.command("pair_style lj/cut 1.9")
        self.lmp.command("pair_coeff 1 1 1.0 1.0")
        self.lmp.command("compute r0 all pair/local dist")
        self.lmp.command("run 0 post no")
        values = self.lmp.numpy.extract_compute("r0", LMP_STYLE_LOCAL, LMP_TYPE_VECTOR)
        self.assertEqual(values.ndim, 1)
        self.assertEqual(values.size, 2)
        self.assertEqual(values[0], 0.5)
        self.assertEqual(values[1], 1.5)

    def testExtractComputeLocalArray(self):
        self.lmp.command("region       box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")
        self.lmp.command("create_atoms 1 single 1.0 1.0 1.0")
        self.lmp.command("create_atoms 1 single 1.0 1.0 1.5")
        self.lmp.command("mass 1 1.0")
        self.lmp.command("pair_style lj/cut 1.9")
        self.lmp.command("pair_coeff 1 1 1.0 1.0")
        self.lmp.command("compute r0 all pair/local dist dx dy dz")
        self.lmp.command("run 0 post no")
        values = self.lmp.numpy.extract_compute("r0", LMP_STYLE_LOCAL, LMP_TYPE_ARRAY)
        self.assertEqual(values.ndim, 2)
        self.assertEqual(values.size, 8)
        self.assertEqual(values[0,0], 0.5)
        self.assertEqual(values[0,3], -0.5)
        self.assertEqual(values[1,0], 1.5)
        self.assertEqual(values[1,3], 1.5)

    def testExtractAtom(self):
        self.lmp.command("units lj")
        self.lmp.command("atom_style atomic")
        self.lmp.command("atom_modify map array")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 2 box")

        x = [ 1.0, 1.0, 1.0,  1.0, 1.0, 1.5,  1.5, 1.0, 1.0 ]
        types = [1, 2, 1]
        ids = [1, 2, 3]
        self.assertEqual(self.lmp.create_atoms(3, id=ids, type=types, x=x), 3)
        self.lmp.command("mass * 2.0")
        self.lmp.command("pair_style zero 1.1")
        self.lmp.command("pair_coeff * *")
        self.lmp.command("fix props all property/atom i_one i2_two 2 d_three d2_four 2");
        self.lmp.command("fix rmass all property/atom mol q rmass ghost yes");
        self.lmp.command("fix 1 all nve")
        self.lmp.command("run 0 post no")
        ntypes = self.lmp.extract_setting("ntypes");
        nlocal = self.lmp.extract_setting("nlocal");
        nall   = self.lmp.extract_setting("nall");
        self.assertEqual(nlocal, 3)
        self.assertEqual(ntypes, 2)
        self.assertEqual(nall, 63)

        self.lmp.command("set atom 1 charge -1");
        self.lmp.command("set atom 2 charge  1");
        self.lmp.command("set atom 3 charge  0");
        self.lmp.command("set atom * mol 2");
        self.lmp.command("set atom 2 mol 1");
        self.lmp.command("set atom 1 i_one -3");
        self.lmp.command("set atom 2 i_one  3");
        self.lmp.command("set atom 2 d_three -1.3");
        self.lmp.command("set atom 3 d_three  3.5");
        self.lmp.command("set atom 1 i_two[1] -3");
        self.lmp.command("set atom 2 i_two[2]  3");
        self.lmp.command("set atom * d_four[1] -1.3");
        self.lmp.command("set atom * d_four[2]  3.5");
        self.lmp.command("run 0 post no")

        mass = self.lmp.numpy.extract_atom("mass")
        self.assertEqual(len(mass), ntypes + 1)
        self.assertTrue((mass == (0.0, 2.0, 2.0)).all())

        rmass = self.lmp.numpy.extract_atom("rmass")
        self.assertEqual(len(rmass), nall)
        self.assertTrue((rmass[0:3] == (0.0, 0.0, 0.0)).all())

        charge = self.lmp.numpy.extract_atom("q")
        self.assertEqual(len(charge), nall)
        self.assertTrue((charge[0:3] == (-1.0, 1.0, 0.0)).all())

        molecule = self.lmp.numpy.extract_atom("molecule")
        self.assertEqual(len(molecule), nall)
        self.assertTrue((molecule[0:3] == (2, 1, 2)).all())

        ident = self.lmp.numpy.extract_atom("id")
        self.assertEqual(len(ident), nall)
        self.assertTrue((ident[0:3] == (1, 2, 3)).all())

        atype = self.lmp.numpy.extract_atom("type")
        self.assertEqual(len(atype), nall)
        self.assertTrue((atype[0:3] == (1, 2, 1)).all())

        x = self.lmp.numpy.extract_atom("x")
        v = self.lmp.numpy.extract_atom("v")
        self.assertEqual(len(x), nall)
        self.assertEqual(len(x[0]), 3)
        self.assertTrue((x[0] == (1.0, 1.0, 1.0)).all())
        self.assertTrue((x[1] == (1.0, 1.0, 1.5)).all())
        self.assertTrue((x[2] == (1.5, 1.0, 1.0)).all())
        self.assertEqual(len(v), nlocal)
        self.assertEqual(len(v[0]), 3)

        self.lmp.command("comm_modify vel yes");
        self.lmp.command("run 0 post no")

        v = self.lmp.numpy.extract_atom("v")
        self.assertEqual(len(v), nall)

        one = self.lmp.numpy.extract_atom("i_one")
        two = self.lmp.numpy.extract_atom("i2_two")
        three = self.lmp.numpy.extract_atom("d_three")
        four = self.lmp.numpy.extract_atom("d2_four")
        self.assertEqual(len(one), nlocal)
        self.assertTrue((one == (-3, 3, 0)).all())
        self.assertEqual(len(two), nlocal)
        self.assertEqual(len(two[0]), 2)
        self.assertTrue((two[0] == (-3, 0)).all())
        self.assertTrue((two[1] == (0, 3)).all())
        self.assertTrue((two[2] == (0, 0)).all())
        self.assertEqual(len(three), nlocal)
        self.assertTrue((three == (0.0, -1.3, 3.5)).all())
        self.assertEqual(len(four), nlocal)
        self.assertEqual(len(four[0]), 2)
        self.assertTrue((four[0] == (-1.3, 3.5)).all())
        self.assertTrue((four[1] == (-1.3, 3.5)).all())
        self.assertTrue((four[2] == (-1.3, 3.5)).all())

    @unittest.skipIf(not has_full,"Gather bonds test")
    def testGatherBond_newton_on(self):
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])
        self.lmp.command("newton on on")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        bonds = self.lmp.numpy.gather_bonds()
        self.assertEqual(len(bonds),24)
        count = 0
        for bond in bonds:
            count += self.checkBond(bond, 5, 1, 2)
            count += self.checkBond(bond, 3, 1, 3)
            count += self.checkBond(bond, 2, 3, 4)
            count += self.checkBond(bond, 2, 3, 5)
            count += self.checkBond(bond, 1, 6, 3)
            count += self.checkBond(bond, 3, 6, 8)
            count += self.checkBond(bond, 4, 6, 7)
            count += self.checkBond(bond, 5, 8, 9)
            count += self.checkBond(bond, 5, 27, 28)
            count += self.checkBond(bond, 5, 29, 27)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full,"Gather bonds test")
    def testGatherBond_newton_off(self):
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])
        self.lmp.command("newton off off")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        bonds = self.lmp.numpy.gather_bonds()
        self.assertEqual(len(bonds),24)
        count = 0
        for bond in bonds:
            count += self.checkBond(bond, 5, 1, 2)
            count += self.checkBond(bond, 3, 1, 3)
            count += self.checkBond(bond, 2, 3, 4)
            count += self.checkBond(bond, 2, 3, 5)
            count += self.checkBond(bond, 1, 6, 3)
            count += self.checkBond(bond, 3, 6, 8)
            count += self.checkBond(bond, 4, 6, 7)
            count += self.checkBond(bond, 5, 8, 9)
            count += self.checkBond(bond, 5, 27, 28)
            count += self.checkBond(bond, 5, 29, 27)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full,"Gather angles test")
    def testGatherAngle_newton_on(self):
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])
        self.lmp.command("newton on on")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        angles = self.lmp.numpy.gather_angles()
        self.assertEqual(len(angles),30)
        count = 0
        for angle in angles:
            count += self.checkAngle(angle, 4, 2, 1, 3)
            count += self.checkAngle(angle, 4, 1, 3, 6)
            count += self.checkAngle(angle, 2, 4, 3, 6)
            count += self.checkAngle(angle, 3, 7, 6, 8)
            count += self.checkAngle(angle, 3, 8, 10, 16)
            count += self.checkAngle(angle, 2, 8, 10, 11)
            count += self.checkAngle(angle, 1, 12, 10, 16)
            count += self.checkAngle(angle, 2, 10, 12, 14)
            count += self.checkAngle(angle, 1, 19, 18, 20)
            count += self.checkAngle(angle, 1, 28, 27, 29)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full,"Gather angles test")
    def testGatherAngle_newton_off(self):
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])
        self.lmp.command("newton off off")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        angles = self.lmp.numpy.gather_angles()
        self.assertEqual(len(angles),30)
        count = 0
        for angle in angles:
            count += self.checkAngle(angle, 4, 2, 1, 3)
            count += self.checkAngle(angle, 4, 1, 3, 6)
            count += self.checkAngle(angle, 2, 4, 3, 6)
            count += self.checkAngle(angle, 3, 7, 6, 8)
            count += self.checkAngle(angle, 3, 8, 10, 16)
            count += self.checkAngle(angle, 2, 8, 10, 11)
            count += self.checkAngle(angle, 1, 12, 10, 16)
            count += self.checkAngle(angle, 2, 10, 12, 14)
            count += self.checkAngle(angle, 1, 19, 18, 20)
            count += self.checkAngle(angle, 1, 28, 27, 29)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full,"Gather dihedrals test")
    def testGatherDihedral_newton_on(self):
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])
        self.lmp.command("newton on on")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        dihedrals = self.lmp.numpy.gather_dihedrals()
        self.assertEqual(len(dihedrals),31)
        count = 0
        for dihedral in dihedrals:
            count += self.checkDihedral(dihedral, 2, 2, 1, 3, 6)
            count += self.checkDihedral(dihedral, 3, 2, 1, 3, 5)
            count += self.checkDihedral(dihedral, 1, 1, 3, 6, 7)
            count += self.checkDihedral(dihedral, 5, 4, 3, 6, 7)
            count += self.checkDihedral(dihedral, 4, 3, 6, 8, 9)
            count += self.checkDihedral(dihedral, 2, 6, 8, 10, 16)
            count += self.checkDihedral(dihedral, 5, 8, 10, 12, 13)
            count += self.checkDihedral(dihedral, 1, 8, 10, 12, 14)
            count += self.checkDihedral(dihedral, 2, 11, 10, 16, 17)
            count += self.checkDihedral(dihedral, 5, 16, 10, 12, 14)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full,"Gather dihedrals test")
    def testGatherDihedral_newton_off(self):
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])
        self.lmp.command("newton off off")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        dihedrals = self.lmp.numpy.gather_dihedrals()
        self.assertEqual(len(dihedrals),31)
        count = 0
        for dihedral in dihedrals:
            count += self.checkDihedral(dihedral, 2, 2, 1, 3, 6)
            count += self.checkDihedral(dihedral, 3, 2, 1, 3, 5)
            count += self.checkDihedral(dihedral, 1, 1, 3, 6, 7)
            count += self.checkDihedral(dihedral, 5, 4, 3, 6, 7)
            count += self.checkDihedral(dihedral, 4, 3, 6, 8, 9)
            count += self.checkDihedral(dihedral, 2, 6, 8, 10, 16)
            count += self.checkDihedral(dihedral, 5, 8, 10, 12, 13)
            count += self.checkDihedral(dihedral, 1, 8, 10, 12, 14)
            count += self.checkDihedral(dihedral, 2, 11, 10, 16, 17)
            count += self.checkDihedral(dihedral, 5, 16, 10, 12, 14)
        self.assertEqual(count,10)

    @unittest.skipIf(not has_full,"Gather impropers test")
    def testGatherImproper_newton_on(self):
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])
        self.lmp.command("newton on on")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        impropers = self.lmp.numpy.gather_impropers()
        self.assertEqual(len(impropers),2)
        count = 0
        for improper in impropers:
            count += self.checkImproper(improper, 1, 6, 3, 8, 7)
            count += self.checkImproper(improper, 2, 8, 6, 10, 9)
        self.assertEqual(count,2)

    @unittest.skipIf(not has_full,"Gather impropers test")
    def testGatherImproper_newton_off(self):
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])
        self.lmp.command("newton off off")
        self.lmp.file("in.fourmol")
        self.lmp.command("run 0 post no")
        impropers = self.lmp.numpy.gather_impropers()
        self.assertEqual(len(impropers),2)
        count = 0
        for improper in impropers:
            count += self.checkImproper(improper, 1, 6, 3, 8, 7)
            count += self.checkImproper(improper, 2, 8, 6, 10, 9)
        self.assertEqual(count,2)

    def testNeighborListSimple(self):
        self.lmp.commands_string("""
        units lj
        atom_style atomic
        atom_modify map array
        boundary f f f
        region box block 0 2 0 2 0 2
        create_box 1 box""")

        x = [ 1.0, 1.0, 1.0,  1.0, 1.0, 1.5 ]
        types = [1, 1]

        self.assertEqual(self.lmp.create_atoms(2, id=None, type=types, x=x), 2)
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 2)

        self.lmp.commands_string("""
        mass 1 1.0
        velocity all create 3.0 87287
        pair_style lj/cut 2.5
        pair_coeff 1 1 1.0 1.0 2.5
        neighbor 0.1 bin
        neigh_modify every 20 delay 0 check no
        run 0 post no""")

        idx = self.lmp.find_pair_neighlist("lj/cut")
        self.assertNotEqual(idx, -1)
        nlist = self.lmp.numpy.get_neighlist(idx)
        self.assertEqual(len(nlist), 2)
        atom_i, neighbors_i = nlist[0]
        atom_j, neighbors_j = nlist[1]

        self.assertEqual(atom_i, 0)
        self.assertEqual(atom_j, 1)

        self.assertEqual(len(neighbors_i), 1)
        self.assertEqual(len(neighbors_j), 0)

        self.assertIn(1, neighbors_i)
        self.assertNotIn(0, neighbors_j)

    def testNeighborListHalf(self):
        self.lmp.commands_string("""
        boundary f f f
        units real
        region box block -5 5 -5 5 -5 5
        create_box 1 box
        mass 1 1.0
        pair_style lj/cut 4.0
        pair_coeff 1 1 0.2 2.0
        """)
        x = [ 0.0,  0.0,  0.0,  -1.1,  0.0,  0.0,  1.0,  0.0,  0.0,
              0.0, -1.1,  0.0,   0.0,  1.0,  0.0,  0.0,  0.0, -1.1,
              0.0,  0.0,  1.0 ]
        tags = [1, 2, 3, 4, 5, 6, 7]
        types = [1, 1, 1, 1, 1, 1, 1]

        self.assertEqual(self.lmp.create_atoms(7, id=tags, type=types, x=x), 7)
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 7)

        self.lmp.command("run 0 post no")

        self.assertEqual(self.lmp.find_pair_neighlist("lj/cut"),0)
        nlist = self.lmp.numpy.get_neighlist(0)
        self.assertEqual(nlist.size, 7)
        for i in range(0,nlist.size):
            idx, neighs = nlist.get(i)
            self.assertEqual(idx,i)
            self.assertEqual(neighs.size,nlocal-1-i)

        # look up neighbor list by atom index
        neighs = nlist.find(2)
        self.assertEqual(neighs.size,4)
        self.assertIsNotNone(neighs,None)
        # this one will fail
        neighs = nlist.find(10)
        self.assertIsNone(neighs,None)

    @unittest.skipIf(not has_manybody,"Full neighbor list test for manybody potential")
    def testNeighborListFull(self):
        self.lmp.commands_string("""
        boundary f f f
        units metal
        region box block -5 5 -5 5 -5 5
        create_box 1 box
        mass 1 1.0
        pair_style sw
        pair_coeff * * Si.sw Si
        """)
        x = [ 0.0,  0.0,  0.0,  -1.1,  0.0,  0.0,  1.0,  0.0,  0.0,
              0.0, -1.1,  0.0,   0.0,  1.0,  0.0,  0.0,  0.0, -1.1,
              0.0,  0.0,  1.0 ]
        tags = [1, 2, 3, 4, 5, 6, 7]
        types = [1, 1, 1, 1, 1, 1, 1]

        self.assertEqual(self.lmp.create_atoms(7, id=tags, type=types, x=x), 7)
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 7)

        self.lmp.command("run 0 post no")

        self.assertEqual(self.lmp.find_pair_neighlist("sw"),0)
        nlist = self.lmp.numpy.get_neighlist(0)
        self.assertEqual(nlist.size, 7)
        for i in range(0,nlist.size):
            idx, neighs = nlist.get(i)
            self.assertEqual(idx,i)
            self.assertEqual(neighs.size,nlocal-1)

    @unittest.skipIf(not has_manybody,"Hybrid neighbor list test for manybody potential")
    def testNeighborListHybrid(self):
        self.lmp.commands_string("""
        boundary f f f
        units metal
        region box block -5 5 -5 5 -5 5
        create_box 2 box
        mass * 1.0
        pair_style hybrid/overlay morse 4.0 lj/cut 4.0 lj/cut 4.0 sw
        pair_coeff * * sw Si.sw Si NULL
        pair_coeff 1 2 morse 0.2 2.0 2.0
        pair_coeff 2 2 lj/cut 1 0.1 2.0
        pair_coeff * * lj/cut 2 0.01 2.0
        """)
        x = [ 0.0,  0.0,  0.0,  -1.1,  0.0,  0.0,  1.0,  0.0,  0.0,
              0.0, -1.1,  0.0,   0.0,  1.0,  0.0,  0.0,  0.0, -1.1,
              0.0,  0.0,  1.0 ]
        tags = [1, 2, 3, 4, 5, 6, 7]
        types = [1, 1, 1, 1, 2, 2, 2]

        self.assertEqual(self.lmp.create_atoms(7, id=tags, type=types, x=x), 7)
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 7)

        self.lmp.command("run 0 post no")

        # valid and invalid lookups
        self.assertNotEqual(self.lmp.find_pair_neighlist("sw"),-1)
        self.assertNotEqual(self.lmp.find_pair_neighlist("morse"),-1)
        self.assertNotEqual(self.lmp.find_pair_neighlist("lj/cut",nsub=1),-1)
        self.assertNotEqual(self.lmp.find_pair_neighlist("lj/cut",nsub=2),-1)
        self.assertEqual(self.lmp.find_pair_neighlist("lj/cut"),-1)
        self.assertEqual(self.lmp.find_pair_neighlist("hybrid/overlay"),-1)
        self.assertNotEqual(self.lmp.numpy.get_neighlist(4).size,0)
        self.assertEqual(self.lmp.numpy.get_neighlist(5).size,-1)

        # full neighbor list for 4 type 1 atoms
        # all have 3 type 1 atom neighbors
        nlist = self.lmp.numpy.get_neighlist(self.lmp.find_pair_neighlist("sw"))
        self.assertEqual(nlist.size, 4)
        for i in range(0,nlist.size):
            idx, neighs = nlist.get(i)
            self.assertEqual(idx,i)
            self.assertEqual(neighs.size,3)

        # half neighbor list for all pairs between type 1 and type 2
        # 4 type 1 atoms with 3 type 2 neighbors and 3 type 2 atoms without neighbors
        nlist = self.lmp.numpy.get_neighlist(self.lmp.find_pair_neighlist("morse"))
        self.assertEqual(nlist.size, 7)
        for i in range(0,nlist.size):
            idx, neighs = nlist.get(i)
            if (i < 4): self.assertEqual(neighs.size,3)
            else: self.assertEqual(neighs.size,0)

        # half neighbor list between type 2 atoms only
        # 3 pairs with 2, 1, 0 neighbors
        nlist = self.lmp.numpy.get_neighlist(self.lmp.find_pair_neighlist("lj/cut",nsub=1))
        self.assertEqual(nlist.size, 3)
        for i in range(0,nlist.size):
            idx, neighs = nlist.get(i)
            self.assertEqual(neighs.size,2-i)

        # half neighbor list between all pairs. same as simple lj/cut case
        nlist = self.lmp.numpy.get_neighlist(self.lmp.find_pair_neighlist("lj/cut",nsub=2))
        self.assertEqual(nlist.size, 7)
        for i in range(0,nlist.size):
            idx, neighs = nlist.get(i)
            self.assertEqual(neighs.size,nlocal-1-i)

    def testNeighborListZeroHalf(self):
        self.lmp.commands_string("""
        boundary f f f
        units real
        region box block -5 5 -5 5 -5 5
        create_box 1 box
        mass 1 1.0
        pair_style zero 4.0
        pair_coeff 1 1
        """)
        x = [ 0.0,  0.0,  0.0,  -1.1,  0.0,  0.0,  1.0,  0.0,  0.0,
              0.0, -1.1,  0.0,   0.0,  1.0,  0.0,  0.0,  0.0, -1.1,
              0.0,  0.0,  1.0 ]
        tags = [1, 2, 3, 4, 5, 6, 7]
        types = [1, 1, 1, 1, 1, 1, 1]

        self.assertEqual(self.lmp.create_atoms(7, id=tags, type=types, x=x), 7)
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 7)

        self.lmp.command("run 0 post no")

        self.assertEqual(self.lmp.find_pair_neighlist("zero"),0)
        nlist = self.lmp.numpy.get_neighlist(0)
        self.assertEqual(nlist.size, 7)
        for i in range(0,nlist.size):
            idx, neighs = nlist.get(i)
            self.assertEqual(idx,i)
            self.assertEqual(neighs.size,nlocal-1-i)

        # look up neighbor list by atom index
        neighs = nlist.find(2)
        self.assertEqual(neighs.size,4)
        self.assertIsNotNone(neighs,None)
        # this one will fail
        neighs = nlist.find(10)
        self.assertIsNone(neighs,None)

    def testNeighborListZeroFull(self):
        self.lmp.commands_string("""
        boundary f f f
        units metal
        region box block -5 5 -5 5 -5 5
        create_box 1 box
        mass 1 1.0
        pair_style zero 4.0 full
        pair_coeff * *
        """)
        x = [ 0.0,  0.0,  0.0,  -1.1,  0.0,  0.0,  1.0,  0.0,  0.0,
              0.0, -1.1,  0.0,   0.0,  1.0,  0.0,  0.0,  0.0, -1.1,
              0.0,  0.0,  1.0 ]
        tags = [1, 2, 3, 4, 5, 6, 7]
        types = [1, 1, 1, 1, 1, 1, 1]

        self.assertEqual(self.lmp.create_atoms(7, id=tags, type=types, x=x), 7)
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 7)

        self.lmp.command("run 0 post no")

        self.assertEqual(self.lmp.find_pair_neighlist("zero"),0)
        nlist = self.lmp.numpy.get_neighlist(0)
        self.assertEqual(nlist.size, 7)
        for i in range(0,nlist.size):
            idx, neighs = nlist.get(i)
            self.assertEqual(idx,i)
            self.assertEqual(neighs.size,nlocal-1)

    def testNeighborListCompute(self):
        self.lmp.commands_string("""
        boundary f f f
        units real
        region box block -5 5 -5 5 -5 5
        create_box 1 box
        mass 1 1.0
        pair_style lj/cut 4.0
        pair_coeff 1 1 0.2 2.0
        compute dist all pair/local dist
        fix dist all ave/histo 1 1 1 0.0 3.0 4 c_dist mode vector
        thermo_style custom f_dist[*]
        """)
        x = [ 0.0,  0.0,  0.0,  -1.1,  0.0,  0.0,  1.0,  0.0,  0.0,
              0.0, -1.1,  0.0,   0.0,  1.0,  0.0,  0.0,  0.0, -1.1,
              0.0,  0.0,  1.0 ]
        tags = [1, 2, 3, 4, 5, 6, 7]
        types = [1, 1, 1, 1, 1, 1, 1]

        self.assertEqual(self.lmp.create_atoms(7, id=tags, type=types, x=x), 7)
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 7)

        self.lmp.command("run 0 post no")
        # check compute data from histogram summary
        nhisto = self.lmp.extract_fix("dist",LMP_STYLE_GLOBAL,LMP_TYPE_VECTOR,nrow=0)
        nskip = self.lmp.extract_fix("dist",LMP_STYLE_GLOBAL,LMP_TYPE_VECTOR,nrow=1)
        minval = self.lmp.extract_fix("dist",LMP_STYLE_GLOBAL,LMP_TYPE_VECTOR,nrow=2)
        maxval = self.lmp.extract_fix("dist",LMP_STYLE_GLOBAL,LMP_TYPE_VECTOR,nrow=3)
        # 21 pair distances counted, none skipped, smallest 1.0, largest 2.1
        self.assertEqual(nhisto,21)
        self.assertEqual(nskip,0)
        self.assertEqual(minval,1.0)
        self.assertEqual(maxval,2.1)

        self.assertNotEqual(self.lmp.find_pair_neighlist("lj/cut"),-1)
        self.assertNotEqual(self.lmp.find_compute_neighlist("dist"),-1)
        self.assertEqual(self.lmp.find_compute_neighlist("xxx"),-1)
        self.assertEqual(self.lmp.find_fix_neighlist("dist"),-1)

        # the compute has a half neighbor list
        nlist = self.lmp.numpy.get_neighlist(self.lmp.find_compute_neighlist("dist"))
        self.assertEqual(nlist.size, 7)
        for i in range(0,nlist.size):
            idx, neighs = nlist.get(i)
            self.assertEqual(idx,i)
            self.assertEqual(neighs.size,nlocal-1-i)

    def test_extract_variable_equalstyle(self):
        self.lmp.command("variable a equal 100")
        a = self.lmp.numpy.extract_variable("a")
        self.assertEqual(a, 100)

        self.lmp.command("variable a equal 3.14")
        a = self.lmp.numpy.extract_variable("a")
        self.assertEqual(a, 3.14)

    def test_extract_variable_atomstyle(self):
        self.lmp.command("units lj")
        self.lmp.command("atom_style atomic")
        self.lmp.command("atom_modify map array")
        self.lmp.command("boundary f f f")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")

        x = [
          1.0, 1.0, 1.0,
          1.0, 1.0, 1.5
        ]

        types = [1, 1]

        self.assertEqual(self.lmp.create_atoms(2, id=None, type=types, x=x), 2)
        self.lmp.command("variable a atom x*x+y*y+z*z")
        a = self.lmp.numpy.extract_variable("a", "all", LMP_VAR_ATOM)
        self.assertIs(type(a), numpy.ndarray)
        self.assertEqual(a[0], x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
        self.assertEqual(a[1], x[3]*x[3]+x[4]*x[4]+x[5]*x[5])

if __name__ == "__main__":
    unittest.main()
