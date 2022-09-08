import os
import sys
import unittest
from lammps import lammps

EXAMPLES_DIR=os.path.abspath(os.path.join(__file__, '..', '..', '..', 'examples'))
PIZZA_DIR=os.path.abspath(os.path.join(__file__, '..', '..', '..', 'tools', 'python', 'pizza'))
DEFAULT_STYLE_EXAMPLE_LOG=os.path.join('melt', 'log.*.melt.g++.1')
MULTI_STYLE_EXAMPLE_LOG=os.path.join('peptide', 'log.27Nov18.peptide.g++.1')
sys.path.insert(1,PIZZA_DIR)

# dump class uses NumPy, so only load and test dump if NumPy is available
has_numpy = False
try:
    import numpy
    has_numpy = True
    import dump
except:
    pass

import log

class Logfiles(unittest.TestCase):
    def testLogFileNotFound(self):
        with self.assertRaises(ValueError):
            l = log.log('test.log')

    def testDefaultLogFile(self):
        l = log.log(os.path.join(EXAMPLES_DIR, DEFAULT_STYLE_EXAMPLE_LOG))
        self.assertEqual(l.nvec, 6)
        self.assertEqual(l.nlen, 6)
        self.assertEqual(l.style, 2)
        self.assertEqual(l.increment, 0)
        n = l.names
        self.assertEqual(len(n), 6)
        self.assertIn("Step", n)
        self.assertIn("Temp", n)
        self.assertIn("E_pair", n)
        self.assertIn("E_mol", n)
        self.assertIn("TotEng", n)
        self.assertIn("Press", n)
        s = l.get("Step")
        self.assertEqual(0.0, s[0])
        self.assertEqual(50.0, s[1])
        self.assertEqual(100.0, s[2])
        self.assertEqual(150.0, s[3])
        self.assertEqual(200.0, s[4])
        self.assertEqual(250.0, s[5])
        t,m = l.get("Temp", "E_mol")
        self.assertEqual(t[0],3.0)
        self.assertEqual(m[2],0.0)
        l.write("all.txt",0)
        l.write("some.txt",1,"Step","Temp","Press")

        with self.assertRaises(Exception):
            t = l.next()

        os.remove('all.txt')
        os.remove('some.txt')

    def testIncrementalLogFile(self):
        l = log.log(os.path.join(EXAMPLES_DIR, DEFAULT_STYLE_EXAMPLE_LOG), 0)
        self.assertEqual(l.style, -1)
        self.assertEqual(l.nvec, 0)
        self.assertEqual(l.increment, 1)
        t = l.next()
        self.assertEqual(l.style, 2)
        self.assertEqual(l.nvec, 6)

    def testMultiLogFile(self):
        l = log.log(os.path.join(EXAMPLES_DIR, MULTI_STYLE_EXAMPLE_LOG))
        self.assertEqual(l.nvec, 14)
        self.assertEqual(l.nlen, 7)
        self.assertEqual(l.style, 1)
        n = l.names
        self.assertEqual(len(n), 14)
        self.assertIn("Step", n)
        self.assertIn("CPU", n)
        self.assertIn("TotEng", n)
        self.assertIn("KinEng", n)
        self.assertIn("Temp", n)
        self.assertIn("PotEng", n)
        self.assertIn("E_bond", n)
        self.assertIn("E_angle", n)
        self.assertIn("E_dihed", n)
        self.assertIn("E_impro", n)
        self.assertIn("E_vdwl", n)
        self.assertIn("E_coul", n)
        self.assertIn("E_long", n)
        self.assertIn("Press", n)
        s = l.get("Step")
        self.assertEqual(0.0, s[0])
        self.assertEqual(50.0, s[1])
        self.assertEqual(100.0, s[2])
        self.assertEqual(150.0, s[3])
        self.assertEqual(200.0, s[4])
        self.assertEqual(250.0, s[5])
        self.assertEqual(300.0, s[6])
        v,c = l.get("E_vdwl", "E_coul")
        self.assertEqual(v[0],692.8945)
        self.assertEqual(c[0],26772.2646)
        l.write("all.txt",0)
        l.write("some.txt",1,"Step","Temp","Press")

        with self.assertRaises(ValueError):
            v = l.get("Volume")
        with self.assertRaises(Exception):
            t = l.next()

        os.remove('all.txt')
        os.remove('some.txt')

class PythonDump(unittest.TestCase):
    def setUp(self):
        machine = None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.lmp = lammps(name=machine,  cmdargs=['-nocite', '-log', 'none', '-echo', 'screen'])

    def tearDown(self):
        del self.lmp

    @unittest.skipIf(not has_numpy,"Missing the NumPy python module")
    def testDumpCustom(self):
        dumpfile = os.path.join(os.path.abspath('.'), 'dump.custom')
        self.lmp.command('shell cd ' + os.environ['TEST_INPUT_DIR'])
        self.lmp.command("newton on on")
        self.lmp.file("in.fourmol")
        self.lmp.command("dump 1 all custom 2 " + dumpfile + " id type mol q x y z vx vy vz")
        self.lmp.command("dump_modify 1 time yes units yes")
        self.lmp.command("run 4 post no")
        self.lmp.command("undump 1")
        d = dump.dump(dumpfile)
        id1, id2 = d.minmax("id")
        self.assertEqual(id1,1)
        self.assertEqual(id2,29)
        t = d.time()
        self.assertEqual(len(t),3)
        d.tselect.one(2,4)
        index, time, flag = d.iterator(0)
        self.assertEqual(index,1)
        self.assertEqual(time,2)
        self.assertEqual(flag,1)
        index, time, flag = d.iterator(1)
        self.assertEqual(index,2)
        self.assertEqual(time,4)
        self.assertEqual(flag,1)
        index, time, flag = d.iterator(1)
        self.assertEqual(index,0)
        self.assertEqual(time,0)
        self.assertEqual(flag,-1)

        with self.assertRaises(Exception):
          t = d.next()

        os.remove(dumpfile)

if __name__ == "__main__":
    unittest.main()
