import os,unittest
from lammps import lammps

try:
    import numpy
    NUMPY_INSTALLED = True
except ImportError:
    NUMPY_INSTALLED = False

@unittest.skipIf(not NUMPY_INSTALLED, "numpy is not available")
class PythonCmdWrapper(unittest.TestCase):
    def setUp(self):
        machine = None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.lmp = lammps(name=machine,  cmdargs=['-nocite', '-log','none', '-echo', 'screen'])
        self.lmp.cmd.units("lj")
        self.lmp.cmd.atom_style("atomic")
        self.lmp.cmd.atom_modify("map array")

        if 'LAMMPS_CMAKE_CACHE' in os.environ:
            self.cmake_cache = {}

            with open(os.environ['LAMMPS_CMAKE_CACHE'], 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#') or line.startswith('//'): continue
                    parts = line.split('=')
                    key, value_type = parts[0].split(':')
                    if len(parts) > 1:
                        value = parts[1]
                        if value_type == "BOOL":
                            value = (value.upper() == "ON")
                    else:
                        value = None
                    self.cmake_cache[key] = value

    def tearDown(self):
        self.lmp.close()
        del self.lmp

    def test_version(self):
        self.assertGreaterEqual(self.lmp.version(), 20200824)

    def test_create_atoms(self):
        self.lmp.cmd.region("box block", 0, 2, 0, 2, 0, 2)
        self.lmp.cmd.create_box(1, "box")

        x = [
          1.0, 1.0, 1.0,
          1.0, 1.0, 1.5
        ]

        types = [1, 1]

        self.assertEqual(self.lmp.create_atoms(2, id=None, type=types, x=x), 2)
        self.assertEqual(self.lmp.extract_global("natoms"), 2)
        pos = self.lmp.numpy.extract_atom("x")
        self.assertEqual(pos.shape[0], 2)
        numpy.testing.assert_array_equal(pos[0], tuple(x[0:3]))
        numpy.testing.assert_array_equal(pos[1], tuple(x[3:6]))

    def test_thermo_capture(self):
        self.lmp.cmd.lattice("fcc", 0.8442),
        self.lmp.cmd.region("box block", 0, 4, 0, 4, 0, 4)
        self.lmp.cmd.create_box(1, "box")
        self.lmp.cmd.create_atoms(1, "box")
        self.lmp.cmd.mass(1, 1.0)
        self.lmp.cmd.velocity("all create", 1.44, 87287, "loop geom")
        self.lmp.cmd.pair_style("lj/cut", 2.5)
        self.lmp.cmd.pair_coeff(1, 1, 1.0, 1.0, 2.5)
        self.lmp.cmd.neighbor(0.3, "bin")
        self.lmp.cmd.neigh_modify("delay 0 every 20 check no")
        self.lmp.cmd.fix("1 all nve")

        current_run = {}

        def append_thermo_data(lmp):
          for k, v in lmp.last_thermo().items():
            current_run.setdefault(k, []).append(v)

        # thermo data is only captured during a run if PYTHON package is enabled
        # without it, it will only capture the final thermo at completion
        nvalues = 1
        if self.lmp.has_package("PYTHON"):
            self.lmp.cmd.fix("myfix", "all", "python/invoke", 10, "end_of_step", append_thermo_data)
            nvalues = 2

        self.lmp.cmd.run(10)
        append_thermo_data(self.lmp)

        for k in ('Step', 'Temp', 'E_pair', 'E_mol', 'TotEng', 'Press'):
            self.assertIn(k, current_run)
            self.assertEqual(len(current_run[k]), nvalues)

if __name__ == "__main__":
    unittest.main()
