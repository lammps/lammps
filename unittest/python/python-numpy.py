import sys,os,unittest
from lammps import lammps, LAMMPS_INT, LMP_STYLE_GLOBAL, LMP_STYLE_LOCAL, LMP_STYLE_ATOM, LMP_TYPE_VECTOR, LMP_TYPE_SCALAR, LMP_TYPE_ARRAY
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

if __name__ == "__main__":
    unittest.main()
