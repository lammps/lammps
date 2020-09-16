import sys,os,unittest
from lammps import PyLammps

class PythonPyLammps(unittest.TestCase):
    def setUp(self):
        machine = None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.pylmp = PyLammps(name=machine,  cmdargs=['-nocite', '-log','none', '-echo','screen'])

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

    def test_create_box(self):
        self.pylmp.region("box block", 0, 2, 0, 2, 0, 2)
        self.pylmp.create_box(1, "box")

        x = [
          1.0, 1.0, 1.0,
          1.0, 1.0, 1.5
        ]

        types = [1, 1]

        self.assertEqual(self.pylmp.lmp.create_atoms(2, id=None, type=types, x=x), 2)
        self.assertEqual(self.pylmp.system.natoms, 2)

if __name__ == "__main__":
    unittest.main()
