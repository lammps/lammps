import os
import unittest
from lammps.formats import LogFile, AvgChunkFile

has_yaml = False
try:
    import yaml
    has_yaml = True
    try:
        from yaml import CSafeLoader as Loader, CSafeDumper as Dumper
    except ImportError:
        from yaml import SafeLoader, SafeDumper
except:
    pass

EXAMPLES_DIR=os.path.abspath(os.path.join(__file__, '..', '..', '..', 'examples'))

DEFAULT_STYLE_EXAMPLE_LOG="melt/log.8Apr21.melt.g++.1"
MULTI_STYLE_EXAMPLE_LOG="peptide/log.27Nov18.peptide.g++.1"
AVG_CHUNK_FILE="VISCOSITY/profile.13Oct16.nemd.2d.g++.1"
YAML_STYLE_EXAMPLE_LOG="yaml/log.7Apr22.yaml.g++.1"

class Logfiles(unittest.TestCase):
    def testLogFileNotFound(self):
        with self.assertRaises(FileNotFoundError):
            LogFile('test.log')

    def testDefaultLogFile(self):
        log = LogFile(os.path.join(EXAMPLES_DIR, DEFAULT_STYLE_EXAMPLE_LOG))
        self.assertEqual(len(log.runs), 1)
        run = log.runs[0]
        self.assertEqual(len(run.keys()), 6)
        self.assertIn("Step", run)
        self.assertIn("Temp", run)
        self.assertIn("E_pair", run)
        self.assertIn("E_mol", run)
        self.assertIn("TotEng", run)
        self.assertIn("Press", run)
        self.assertEqual(len(run["Step"]), 6)
        self.assertEqual(len(run["Temp"]), 6)
        self.assertEqual(len(run["E_pair"]), 6)
        self.assertEqual(len(run["E_mol"]), 6)
        self.assertEqual(len(run["TotEng"]), 6)
        self.assertEqual(len(run["Press"]), 6)
        self.assertEqual(log.runs[0]["Step"], [0, 50, 100, 150, 200, 250])

    def testMultiLogFile(self):
        log = LogFile(os.path.join(EXAMPLES_DIR, MULTI_STYLE_EXAMPLE_LOG))
        self.assertEqual(len(log.runs), 1)
        run0 = log.runs[0]

        self.assertEqual(len(run0.keys()), 14)
        self.assertIn("Step", run0)
        self.assertIn("CPU", run0)
        self.assertIn("TotEng", run0)
        self.assertIn("KinEng", run0)
        self.assertIn("Temp", run0)
        self.assertIn("PotEng", run0)
        self.assertIn("E_bond", run0)
        self.assertIn("E_angle", run0)
        self.assertIn("E_dihed", run0)
        self.assertIn("E_impro", run0)
        self.assertIn("E_vdwl", run0)
        self.assertIn("E_coul", run0)
        self.assertIn("E_long", run0)
        self.assertIn("Press", run0)

        for k in run0:
            self.assertEqual(len(run0[k]), 7)

        self.assertEqual(run0["Step"], list(range(0,350, 50)))

    @unittest.skipIf(not has_yaml,"Missing the PyYAML python module")
    def testYamlLogFile(self):
        log = LogFile(os.path.join(EXAMPLES_DIR, YAML_STYLE_EXAMPLE_LOG))
        self.assertEqual(len(log.runs), 2)
        run = log.runs[0]
        self.assertEqual(len(run.keys()), 12)
        self.assertIn("Step", run)
        self.assertIn("Temp", run)
        self.assertIn("E_vdwl", run)
        self.assertIn("E_coul", run)
        self.assertIn("E_bond", run)
        self.assertIn("E_angle", run)
        self.assertIn("Press", run)
        self.assertEqual(len(run["Step"]), 11)
        self.assertEqual(len(run["Temp"]), 11)
        self.assertEqual(len(run["E_vdwl"]), 11)
        self.assertEqual(len(run["E_coul"]), 11)
        self.assertEqual(len(run["E_bond"]), 11)
        self.assertEqual(len(run["E_angle"]), 11)
        self.assertEqual(len(run["Press"]), 11)
        self.assertEqual(log.runs[0]["Step"], [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])


class AvgChunkFiles(unittest.TestCase):
    def testAvgChunkFileNotFound(self):
        with self.assertRaises(FileNotFoundError):
            AvgChunkFile('test.log')

    def testRead(self):
        cfile = AvgChunkFile(os.path.join(EXAMPLES_DIR, AVG_CHUNK_FILE))
        self.assertEqual(cfile.fix_name, "4")
        self.assertEqual(cfile.group_name, "all")
        self.assertEqual(cfile.timesteps, list(range(10000, 110000, 5000)))

        ntimesteps = len(cfile.timesteps)
        ntotal_count = len(cfile.total_count)
        nchunks = len(cfile.chunks)
        self.assertEqual(ntimesteps, ntotal_count)
        self.assertEqual(nchunks, 20)

        for i in range(1, nchunks+1):
            chunk  = cfile.chunks[i-1];
            self.assertEqual(chunk['id'], i)
            self.assertEqual(len(chunk['coord']), ntimesteps)
            self.assertEqual(len(chunk['ncount']), ntimesteps)
            self.assertIn("vx", chunk)
            self.assertEqual(len(chunk['vx']), ntimesteps)

        self.assertEqual(len(chunk['coord'][0]), 1)


from lammps import lammps
has_dump_yaml = False
try:
    machine=None
    if 'LAMMPS_MACHINE_NAME' in os.environ:
        machine=os.environ['LAMMPS_MACHINE_NAME']
    lmp=lammps(name=machine)
    has_dump_yaml = lmp.has_style("atom","full") and lmp.has_style("dump", "yaml")
    lmp.close()
except:
    pass

@unittest.skipIf(not (has_dump_yaml and has_yaml), "Either atom_style full, dump_style yaml, or the python PyYAML module are not available")
class PythonDump(unittest.TestCase):
    def setUp(self):
        machine = None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.lmp = lammps(name=machine,  cmdargs=['-nocite', '-log','none', '-echo','screen'])

    def tearDown(self):
        del self.lmp

    def testDumpYaml(self):
        dumpfile = os.path.join(os.path.abspath('.'), 'dump.yaml')
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])
        self.lmp.command("newton on on")
        self.lmp.file("in.fourmol")
        self.lmp.command("dump 1 all yaml 2 " + dumpfile + " id type mol q x y z vx vy vz")
        self.lmp.command("dump_modify 1 time yes sort id units yes")
        self.lmp.command("run 4 post no")
        with open(dumpfile) as d:
            traj = tuple(yaml.load_all(d, Loader=Loader))
        self.assertEqual(len(traj), 3)
        self.assertEqual(traj[0]['timestep'], 0)
        self.assertEqual(traj[0]['time'], 0)
        self.assertEqual(traj[0]['natoms'], 29)
        self.assertEqual(traj[0]['units'], 'real')
        self.assertEqual(len(traj[0]['boundary']), 6)
        self.assertEqual(traj[0]['boundary'][0], 'p')
        self.assertEqual(traj[1]['timestep'], 2)
        self.assertEqual(traj[1]['time'], 0.2)
        self.assertEqual(traj[2]['timestep'], 4)
        self.assertEqual(traj[2]['time'], 0.4)
        self.assertEqual(traj[0]['keywords'],['id', 'type', 'mol', 'q', 'x', 'y', 'z',
                                              'vx', 'vy', 'vz'])
        self.assertEqual(traj[0]['data'][0],[1, 3, 1, -0.47, -0.279937, 2.47266, -0.172009,
                                             0.000778678, 0.000589703, -0.000221795])

if __name__ == "__main__":
    unittest.main()
