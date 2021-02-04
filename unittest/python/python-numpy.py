import sys,os,unittest
from lammps import lammps, LAMMPS_INT, LMP_STYLE_GLOBAL, LMP_STYLE_LOCAL, \
                   LMP_STYLE_ATOM, LMP_TYPE_VECTOR, LMP_TYPE_SCALAR, LMP_TYPE_ARRAY, \
                   LMP_VAR_ATOM
from ctypes import c_void_p

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

    def testExtractComputeLocalScalar(self):
        # TODO
        pass

    def testExtractComputeLocalVector(self):
        # TODO
        pass

    def testExtractComputeLocalArray(self):
        # TODO
        pass

    def testExtractAtomDeprecated(self):
        self.lmp.command("units lj")
        self.lmp.command("atom_style atomic")
        self.lmp.command("atom_modify map array")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")

        x = [
          1.0, 1.0, 1.0,
          1.0, 1.0, 1.5
        ]

        types = [1, 1]

        self.assertEqual(self.lmp.create_atoms(2, id=None, type=types, x=x), 2)
        nlocal = self.lmp.extract_global("nlocal", LAMMPS_INT)
        self.assertEqual(nlocal, 2)

        ident = self.lmp.numpy.extract_atom_iarray("id", nlocal, dim=1)
        self.assertEqual(len(ident), 2)

        ntypes = self.lmp.extract_global("ntypes", LAMMPS_INT)
        self.assertEqual(ntypes, 1)

        x = self.lmp.numpy.extract_atom_darray("x", nlocal, dim=3)
        v = self.lmp.numpy.extract_atom_darray("v", nlocal, dim=3)
        self.assertEqual(len(x), 2)
        self.assertTrue((x[0] == (1.0, 1.0, 1.0)).all())
        self.assertTrue((x[1] == (1.0, 1.0, 1.5)).all())
        self.assertEqual(len(v), 2)

    def testExtractAtom(self):
        self.lmp.command("units lj")
        self.lmp.command("atom_style atomic")
        self.lmp.command("atom_modify map array")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")

        x = [
          1.0, 1.0, 1.0,
          1.0, 1.0, 1.5
        ]

        types = [1, 1]

        self.assertEqual(self.lmp.create_atoms(2, id=None, type=types, x=x), 2)
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 2)

        ident = self.lmp.numpy.extract_atom("id")
        self.assertEqual(len(ident), 2)

        ntypes = self.lmp.extract_global("ntypes")
        self.assertEqual(ntypes, 1)

        x = self.lmp.numpy.extract_atom("x")
        v = self.lmp.numpy.extract_atom("v")
        self.assertEqual(len(x), 2)
        self.assertTrue((x[0] == (1.0, 1.0, 1.0)).all())
        self.assertTrue((x[1] == (1.0, 1.0, 1.5)).all())
        self.assertEqual(len(v), 2)

    def testNeighborList(self):
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
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 2)

        self.lmp.command("mass 1 1.0")
        self.lmp.command("velocity all create 3.0 87287")
        self.lmp.command("pair_style lj/cut 2.5")
        self.lmp.command("pair_coeff 1 1 1.0 1.0 2.5")
        self.lmp.command("neighbor 0.1 bin")
        self.lmp.command("neigh_modify every 20 delay 0 check no")

        self.lmp.command("run 0")

        self.assertEqual(self.lmp.find_pair_neighlist("lj/cut"), 0)
        nlist = self.lmp.numpy.get_neighlist(0)
        self.assertEqual(len(nlist), 2)
        atom_i, neighbors_i = nlist[0]
        atom_j, neighbors_j = nlist[1]

        self.assertEqual(atom_i, 0)
        self.assertEqual(atom_j, 1)

        self.assertEqual(len(neighbors_i), 1)
        self.assertEqual(len(neighbors_j), 0)

        self.assertIn(1, neighbors_i)
        self.assertNotIn(0, neighbors_j)

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
