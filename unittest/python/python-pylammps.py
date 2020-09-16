import sys,os,unittest
from lammps import PyLammps

class PythonPyLammps(unittest.TestCase):
    def setUp(self):
        machine = None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.pylmp = PyLammps(name=machine,  cmdargs=['-nocite', '-log','none', '-echo', 'screen'])
        self.pylmp.units("lj")
        self.pylmp.atom_style("atomic")
        self.pylmp.atom_modify("map array")

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
        self.pylmp.close()
        del self.pylmp

    def test_version(self):
        self.assertGreaterEqual(self.pylmp.version(), 20200824)

    def test_create_atoms(self):
        self.pylmp.region("box block", 0, 2, 0, 2, 0, 2)
        self.pylmp.create_box(1, "box")

        x = [
          1.0, 1.0, 1.0,
          1.0, 1.0, 1.5
        ]

        types = [1, 1]

        self.assertEqual(self.pylmp.lmp.create_atoms(2, id=None, type=types, x=x), 2)
        self.assertEqual(self.pylmp.system.natoms, 2)
        self.assertEqual(len(self.pylmp.atoms), 2)
        self.assertEqual(self.pylmp.atoms[0].position, tuple(x[0:3]))
        self.assertEqual(self.pylmp.atoms[1].position, tuple(x[3:6]))
        self.assertEqual(self.pylmp.last_run, None)


    def test_write_script(self):
        outfile = 'in.test_write_script'
        self.pylmp.write_script(outfile)
        self.assertTrue(os.path.exists(outfile))
        os.remove(outfile)

    def test_runs(self):
        self.pylmp.lattice("fcc", 0.8442),
        self.pylmp.region("box block", 0, 4, 0, 4, 0, 4)
        self.pylmp.create_box(1, "box")
        self.pylmp.create_atoms(1, "box")
        self.pylmp.mass(1, 1.0)
        self.pylmp.velocity("all create", 1.44, 87287, "loop geom")
        self.pylmp.pair_style("lj/cut", 2.5)
        self.pylmp.pair_coeff(1, 1, 1.0, 1.0, 2.5)
        self.pylmp.neighbor(0.3, "bin")
        self.pylmp.neigh_modify("delay 0 every 20 check no")
        self.pylmp.fix("1 all nve")
        self.pylmp.variable("fx atom fx")
        self.pylmp.run(10)

        self.assertEqual(len(self.pylmp.runs), 1)
        self.assertEqual(self.pylmp.last_run, self.pylmp.runs[0])
        self.assertEqual(len(self.pylmp.last_run.thermo.Step), 2)
        self.assertEqual(len(self.pylmp.last_run.thermo.Temp), 2)
        self.assertEqual(len(self.pylmp.last_run.thermo.E_pair), 2)
        self.assertEqual(len(self.pylmp.last_run.thermo.E_mol), 2)
        self.assertEqual(len(self.pylmp.last_run.thermo.TotEng), 2)
        self.assertEqual(len(self.pylmp.last_run.thermo.Press), 2)

if __name__ == "__main__":
    unittest.main()
