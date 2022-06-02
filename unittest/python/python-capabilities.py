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

    def test_os_info(self):
        import platform

        system = platform.system()
        osinfo = self.lmp.get_os_info()
        print("System: %s   LAMMPS OS Info: %s" % (system, osinfo))
        self.assertEqual(osinfo.find(system),0)

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

    def test_has_id(self):
        self.lmp.command('fix charge all property/atom q ghost yes')
        self.lmp.command('region box block 0 1 0 1 0 1')
        self.lmp.command('create_box 1 box')
        self.lmp.command('group none empty')
        self.lmp.command('variable test index 1')
        self.assertTrue(self.lmp.has_id('compute', 'thermo_temp'))
        self.assertTrue(self.lmp.has_id('compute', 'thermo_press'))
        self.assertTrue(self.lmp.has_id('compute', 'thermo_pe'))
        self.assertFalse(self.lmp.has_id('compute', 'xxx'))
        self.assertFalse(self.lmp.has_id('dump', 'xxx'))
        self.assertTrue(self.lmp.has_id('fix', 'charge'))
        self.assertFalse(self.lmp.has_id('fix', 'xxx'))
        self.assertTrue(self.lmp.has_id('group', 'all'))
        self.assertTrue(self.lmp.has_id('group', 'none'))
        self.assertFalse(self.lmp.has_id('group', 'xxx'))
        self.assertTrue(self.lmp.has_id('region', 'box'))
        self.assertFalse(self.lmp.has_id('region', 'xxx'))
        self.assertTrue(self.lmp.has_id('variable', 'test'))
        self.assertFalse(self.lmp.has_id('variable', 'xxx'))

    def test_available_id(self):
        self.lmp.command('fix charge all property/atom q ghost yes')
        self.lmp.command('region box block 0 1 0 1 0 1')
        self.lmp.command('create_box 1 box')
        self.lmp.command('group none empty')
        self.lmp.command('variable test index 1')
        ids = self.lmp.available_ids('compute')
        self.assertIn('thermo_pe', ids)
        self.assertEqual(len(ids),3)
        ids = self.lmp.available_ids('dump')
        self.assertEqual(len(ids),0)
        ids = self.lmp.available_ids('fix')
        self.assertIn('charge', ids)
        self.assertEqual(len(ids),1)
        ids = self.lmp.available_ids('group')
        self.assertIn('none', ids)
        self.assertEqual(len(ids),2)
        ids = self.lmp.available_ids('molecule')
        self.assertEqual(len(ids),0)
        ids = self.lmp.available_ids('region')
        self.assertIn('box', ids)
        self.assertEqual(len(ids),1)
        ids = self.lmp.available_ids('variable')
        self.assertIn('test', ids)
        self.assertEqual(len(ids),1)

    def test_is_running(self):
        self.assertFalse(self.lmp.is_running)

    def test_force_timeout(self):
        self.lmp.command('region box block 0 1 0 1 0 1')
        self.lmp.command('create_box 1 box')
        self.lmp.command('mass * 1.0')
        self.lmp.command('run 10')
        self.assertEqual(self.lmp.extract_global('ntimestep'),10)
        self.lmp.force_timeout()
        self.lmp.command('run 10')
        self.assertEqual(self.lmp.extract_global('ntimestep'),10)
        self.lmp.command('timer timeout off')
        self.lmp.command('run 10')
        self.assertEqual(self.lmp.extract_global('ntimestep'),20)

    def test_accelerator_config(self):

        settings = self.lmp.accelerator_config
        if self.cmake_cache['PKG_OPENMP']:
            if self.cmake_cache['BUILD_OMP']:
                self.assertIn('openmp',settings['OPENMP']['api'])
            else:
                self.assertIn('serial',settings['OPENMP']['api'])
            self.assertIn('double',settings['OPENMP']['precision'])

        if self.cmake_cache['PKG_INTEL']:
            if 'LMP_INTEL_OFFLOAD' in self.cmake_cache.keys():
                self.assertIn('phi',settings['INTEL']['api'])
            elif self.cmake_cache['BUILD_OMP']:
                self.assertIn('openmp',settings['INTEL']['api'])
            else:
                self.assertIn('serial',settings['INTEL']['api'])
            self.assertIn('double',settings['INTEL']['precision'])
            self.assertIn('mixed',settings['INTEL']['precision'])
            self.assertIn('single',settings['INTEL']['precision'])

        if self.cmake_cache['PKG_GPU']:
            if self.cmake_cache['GPU_API'].lower() == 'opencl':
                 self.assertIn('opencl',settings['GPU']['api'])
            if self.cmake_cache['GPU_API'].lower() == 'cuda':
                 self.assertIn('cuda',settings['GPU']['api'])
            if self.cmake_cache['GPU_API'].lower() == 'hip':
                 self.assertIn('hip',settings['GPU']['api'])
            if self.cmake_cache['GPU_PREC'].lower() == 'double':
                 self.assertIn('double',settings['GPU']['precision'])
            if self.cmake_cache['GPU_PREC'].lower() == 'mixed':
                 self.assertIn('mixed',settings['GPU']['precision'])
            if self.cmake_cache['GPU_PREC'].lower() == 'single':
                 self.assertIn('single',settings['GPU']['precision'])

        if self.cmake_cache['PKG_KOKKOS']:
            if 'Kokkos_ENABLE_OPENMP' in self.cmake_cache and self.cmake_cache['Kokkos_ENABLE_OPENMP']:
                self.assertIn('openmp',settings['KOKKOS']['api'])
            if 'Kokkos_ENABLE_SERIAL' in self.cmake_cache and self.cmake_cache['Kokkos_ENABLE_SERIAL']:
                self.assertIn('serial',settings['KOKKOS']['api'])
            self.assertIn('double',settings['KOKKOS']['precision'])

    def test_gpu_device(self):

        info = self.lmp.get_gpu_device_info()
        if self.lmp.has_gpu_device:
            self.assertTrue(info)
            self.assertGreaterEqual(info.find("Device"),0)
        else:
            self.assertFalse(info)

if __name__ == "__main__":
    unittest.main()
