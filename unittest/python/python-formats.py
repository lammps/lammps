import os
import unittest
from lammps.formats import LogFile, AvgChunkFile

EXAMPLES_DIR=os.path.abspath(os.path.join(__file__, '..', '..', '..', 'examples'))

DEFAULT_STYLE_EXAMPLE_LOG="melt/log.8Apr21.melt.g++.1"
MULTI_STYLE_EXAMPLE_LOG="peptide/log.27Nov18.peptide.g++.1"
AVG_CHUNK_FILE="VISCOSITY/profile.13Oct16.nemd.2d.g++.1"

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


if __name__ == "__main__":
    unittest.main()
