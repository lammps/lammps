import sys,os,unittest
from lammps import lammps

class PythonCapabilities(unittest.TestCase):
    def setUp(self):
        machine = None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.lmp = lammps(name=machine,  cmdargs=['-nocite', '-log','none', '-echo','screen'])

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
                            value = (value == "ON")
                    else:
                        value = None
                    self.cmake_cache[key] = value

    def tearDown(self):
        del self.lmp

    def test_has_gzip_support(self):
        self.assertEqual(self.lmp.has_gzip_support, self.cmake_cache['WITH_GZIP'])

    def test_has_png_support(self):
        self.assertEqual(self.lmp.has_png_support, self.cmake_cache['WITH_PNG'])

    def test_has_jpeg_support(self):
        self.assertEqual(self.lmp.has_jpeg_support, self.cmake_cache['WITH_JPEG'])

if __name__ == "__main__":
    unittest.main()
