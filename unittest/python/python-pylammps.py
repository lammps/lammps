import os,unittest
from lammps import PyLammps

try:
    import numpy
    NUMPY_INSTALLED = True
except ImportError:
    NUMPY_INSTALLED = False

@unittest.skipIf(not NUMPY_INSTALLED, "numpy is not available")
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
        numpy.testing.assert_array_equal(self.pylmp.atoms[0].position, tuple(x[0:3]))
        numpy.testing.assert_array_equal(self.pylmp.atoms[1].position, tuple(x[3:6]))
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

        # thermo data is only captured during a run if PYTHON package is enabled
        # without it, it will only capture the final thermo at completion
        if self.pylmp.lmp.has_package("PYTHON"):
            self.assertEqual(len(self.pylmp.runs), 1)
            self.assertEqual(self.pylmp.last_run, self.pylmp.runs[0])
            self.assertEqual(len(self.pylmp.last_run.thermo.Step), 2)
            self.assertEqual(len(self.pylmp.last_run.thermo.Temp), 2)
            self.assertEqual(len(self.pylmp.last_run.thermo.E_pair), 2)
            self.assertEqual(len(self.pylmp.last_run.thermo.E_mol), 2)
            self.assertEqual(len(self.pylmp.last_run.thermo.TotEng), 2)
            self.assertEqual(len(self.pylmp.last_run.thermo.Press), 2)
        else:
            self.assertEqual(len(self.pylmp.runs), 1)
            self.assertEqual(self.pylmp.last_run, self.pylmp.runs[0])
            self.assertEqual(len(self.pylmp.last_run.thermo.Step), 1)
            self.assertEqual(len(self.pylmp.last_run.thermo.Temp), 1)
            self.assertEqual(len(self.pylmp.last_run.thermo.E_pair), 1)
            self.assertEqual(len(self.pylmp.last_run.thermo.E_mol), 1)
            self.assertEqual(len(self.pylmp.last_run.thermo.TotEng), 1)
            self.assertEqual(len(self.pylmp.last_run.thermo.Press), 1)

    def test_info_queries(self):
        self.pylmp.lattice("fcc", 0.8442),
        self.pylmp.region("box block", 0, 4, 0, 4, 0, 4)
        self.pylmp.create_box(1, "box")
        self.pylmp.variable("a equal 10.0")
        self.pylmp.variable("b string value")
        self.assertEqual(self.pylmp.variables['a'].value, 10.0)
        self.assertEqual(self.pylmp.variables['b'].value, 'value')
        self.assertEqual(len(self.pylmp.variables),2)
        self.assertEqual(self.pylmp.system.units,'lj')
        self.assertEqual(self.pylmp.system.atom_style,'atomic')
        self.assertEqual(self.pylmp.system.ntypes,1)
        self.assertEqual(self.pylmp.system.natoms,0)
        self.assertEqual(self.pylmp.communication.comm_style,'brick')
        self.assertEqual(self.pylmp.communication.comm_layout,'uniform')
        self.assertEqual(self.pylmp.communication.nprocs,1)
        self.assertEqual(self.pylmp.communication.nthreads,1)
        self.assertEqual(self.pylmp.communication.procgrid,[1,1,1])
        self.assertEqual(self.pylmp.communication.proc_grid,[1,1,1])
        self.assertEqual(self.pylmp.communication.ghost_velocity,0)
        self.assertEqual(len(self.pylmp.computes),3)
        self.assertEqual(self.pylmp.computes[0]['name'], 'thermo_temp')
        self.assertEqual(self.pylmp.computes[0]['style'], 'temp')
        self.assertEqual(self.pylmp.computes[0]['group'], 'all')
        self.assertEqual(self.pylmp.computes[1]['name'], 'thermo_press')
        self.assertEqual(self.pylmp.computes[1]['style'], 'pressure')
        self.assertEqual(self.pylmp.computes[1]['group'], 'all')
        self.assertEqual(self.pylmp.computes[2]['name'], 'thermo_pe')
        self.assertEqual(self.pylmp.computes[2]['style'], 'pe')
        self.assertEqual(self.pylmp.computes[2]['group'], 'all')
        self.assertEqual(len(self.pylmp.dumps),0)
        self.pylmp.fix('one','all','nve')
        self.assertEqual(len(self.pylmp.fixes),1)
        self.assertEqual(self.pylmp.fixes[0]['name'], 'one')
        self.assertEqual(self.pylmp.fixes[0]['style'], 'nve')
        self.assertEqual(self.pylmp.fixes[0]['group'], 'all')
        self.pylmp.group('none','empty')
        self.assertEqual(len(self.pylmp.groups),2)
        self.pylmp.comm_style('tiled')
        self.pylmp.mass('*',1.0)
        self.pylmp.run('0','post','no')
        self.pylmp.balance(0.1,'rcb')
        self.assertEqual(self.pylmp.communication.procgrid,None)

if __name__ == "__main__":
    unittest.main()
