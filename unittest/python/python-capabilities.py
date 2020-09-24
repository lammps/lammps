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
                    if value_type == "UNINITIALIZED": continue
                    if value_type == "INTERNAL": continue
                    if len(parts) > 1:
                        value = parts[1]
                        if value_type == "BOOL":
                            value = (value.upper() == "ON") or (value.upper() == "YES") or (value == "1")
                    else:
                        value = None
                    self.cmake_cache[key] = value

    def tearDown(self):
        del self.lmp

    def test_version(self):
        self.assertGreaterEqual(self.lmp.version(), 20200824)

    def test_has_gzip_support(self):
        self.assertEqual(self.lmp.has_gzip_support, self.cmake_cache['WITH_GZIP'])

    def test_has_png_support(self):
        self.assertEqual(self.lmp.has_png_support, self.cmake_cache['WITH_PNG'])

    def test_has_jpeg_support(self):
        self.assertEqual(self.lmp.has_jpeg_support, self.cmake_cache['WITH_JPEG'])

    def test_has_ffmpeg_support(self):
        self.assertEqual(self.lmp.has_ffmpeg_support, self.cmake_cache['WITH_FFMPEG'])

    def test_installed_packages(self):
        installed_packages = self.lmp.installed_packages
        selected_packages = [key[4:] for key in self.cmake_cache.keys() if not key.startswith('PKG_CONFIG') and key.startswith('PKG_') and self.cmake_cache[key]]

        for pkg in selected_packages:
            self.assertIn(pkg, installed_packages)

    def test_has_style(self):
        self.assertTrue(self.lmp.has_style('pair', 'lj/cut'))
        self.assertFalse(self.lmp.has_style('pair', 'lennard_jones'))

    def test_available_styles(self):
        pairs = self.lmp.available_styles('pair')
        self.assertIn('lj/cut', pairs)

if __name__ == "__main__":
    unittest.main()
