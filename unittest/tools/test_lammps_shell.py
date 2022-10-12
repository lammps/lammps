
import os, re, subprocess, unittest

# enable test mode
os.putenv('LAMMPS_SHELL_TESTING','1')

shell_prompt_re = r"([^>]*LAMMPS Shell> ([a-z0-9_]+) *([a-z0-9_\.]+)?.*\n)+"
cmd_group_re = r"([^>]*LAMMPS Shell> ([a-z0-9_]+) +([a-z0-9]+) +([a-z0-9]+)? *([a-z/0-9]+)?.*\n)+"

#
class LammpsShell(unittest.TestCase):

    def setUp(self):
        self.proc = subprocess.Popen('./lammps-shell',
                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE)


    def tearDown(self):
        self.proc.kill()


    def InputRunner(self,text):
        """Test tab expansions"""
        try:
            [outs,errs] = self.proc.communicate(input=text, timeout=10)
            self.timeout = 0
        except subprocess.TimeoutExpired:
            self.proc.kill()
            [outs,errs] = self.proc.communicate()
            self.timeout = 1

        return outs.decode('UTF-8')

    def testExpandClearHistory(self):
        """Test expansion of a shell specific command"""
        matches = re.findall(shell_prompt_re, self.InputRunner(b'clear_his\t\n'), re.MULTILINE)
        if self.timeout:
            self.fail("Timeout")
        else:
            self.assertEqual(matches[0][1],"clear_history")

    def testExpandDimension(self):
        """Test expansion of a LAMMPS command"""
        matches = re.findall(shell_prompt_re, self.InputRunner(b'dimens\t\n'), re.MULTILINE)
        if self.timeout:
            self.fail("Timeout")
        else:
            self.assertEqual(matches[0][1],"dimension")

    def testExpandPairStyle(self):
        """Test expansion of a pair style"""
        matches = re.findall(shell_prompt_re, self.InputRunner(b'pair_st\t zer\t\n'), re.MULTILINE)
        if self.timeout:
            self.fail("Timeout")
        else:
            self.assertEqual(matches[0][1],"pair_style")
            self.assertEqual(matches[0][2],"zero")

    def testExpandBondStyle(self):
        """Test expansion of a bond style"""
        matches = re.findall(shell_prompt_re, self.InputRunner(b'bond_st\t zer\t\n'), re.MULTILINE)
        if self.timeout:
            self.fail("Timeout")
        else:
            self.assertEqual(matches[0][1],"bond_style")
            self.assertEqual(matches[0][2],"zero")

    def testExpandAngleStyle(self):
        """Test expansion of a angle style"""
        matches = re.findall(shell_prompt_re, self.InputRunner(b'angle_st\t zer\t\n'), re.MULTILINE)
        if self.timeout:
            self.fail("Timeout")
        else:
            self.assertEqual(matches[0][1],"angle_style")
            self.assertEqual(matches[0][2],"zero")

    def testExpandDihedralStyle(self):
        """Test expansion of a dihedral style"""
        matches = re.findall(shell_prompt_re, self.InputRunner(b'dihedral_st\t zer\t\n'), re.MULTILINE)
        if self.timeout:
            self.fail("Timeout")
        else:
            self.assertEqual(matches[0][1],"dihedral_style")
            self.assertEqual(matches[0][2],"zero")

    def testExpandImproperStyle(self):
        """Test expansion of a improper style"""
        matches = re.findall(shell_prompt_re, self.InputRunner(b'improper_st\t zer\t\n'), re.MULTILINE)
        if self.timeout:
            self.fail("Timeout")
        else:
            self.assertEqual(matches[0][1],"improper_style")
            self.assertEqual(matches[0][2],"zero")

    def testExpandComputeGroup(self):
        """Test expansion of a group-ID and a compute command"""
        matches = re.findall(cmd_group_re, self.InputRunner(b'compute test al\tstress/at\t\n'), re.MULTILINE)
        if self.timeout:
            self.fail("Timeout")
        else:
            self.assertEqual(matches[0][1],"compute")
            self.assertEqual(matches[0][2],"test")
            self.assertEqual(matches[0][3],"all")
            self.assertEqual(matches[0][4],"stress/atom")

    def testExpandFixGroup(self):
        """Test expansion of a group-ID and a fix command"""
        matches = re.findall(cmd_group_re, self.InputRunner(b'fix test al\tpropert\t\n'), re.MULTILINE)
        if self.timeout:
            self.fail("Timeout")
        else:
            self.assertEqual(matches[0][1],"fix")
            self.assertEqual(matches[0][2],"test")
            self.assertEqual(matches[0][3],"all")
            self.assertEqual(matches[0][4],"property/atom")

    def testExpandSource(self):
        """Test expansion of a shell command and a file name"""
        with open('.tmp.in.source', 'w') as out:
           print('units real', file=out)
        out.close()
        matches = re.findall(shell_prompt_re, self.InputRunner(b'sour\t.tmp.in.sou\t\n'), re.MULTILINE)
        os.remove('.tmp.in.source')
        if self.timeout:
            self.fail("Timeout")
        else:
            self.assertEqual(matches[0][1],"source")
            self.assertEqual(matches[0][2],".tmp.in.source")

    def testHistory(self):
        """Test history expansion"""
        out = self.InputRunner(b'clear_history\nunits real\ndimension 2\n!!:p\n!-3:p\n!dim:p\n!uni:p\nprint !!:$\nprint !dim:1\n')
        idx = 0
        if self.timeout:
            self.fail("Timeout")
        else:
            lines = out.splitlines()
            for line in lines:
                if line.startswith('LAMMPS Shell>'): break
                idx += 1

            self.assertEqual(lines[idx+4],"dimension 2")
            self.assertEqual(lines[idx+6],"units real")
            self.assertEqual(lines[idx+8],"dimension 2")
            self.assertEqual(lines[idx+10],"units real")
            self.assertEqual(lines[idx+12],"real")
            self.assertEqual(lines[idx+14],"2")

###########################
if __name__ == "__main__":
    unittest.main()
