
import os, re, subprocess, unittest

# enable test mode
os.putenv('LAMMPS_SHELL_TESTING','1')

shell_prompt_re = r"([^>]*LAMMPS Shell> ([a-z0-9_]+) *([a-z0-9_\.]+)?.*\n)+"

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
            [outs,errs] = self.proc.communicate(input=text, timeout=1)
        except TimeoutExpired:
            proc.kill()
            [outs,errs] = proc.communicate()

        return outs.decode('UTF-8')

    def testExpandClear(self):
        """Test expansion of a shell specific command"""
        matches = re.findall(shell_prompt_re, self.InputRunner(b'clear_his\t\n'), re.MULTILINE)
        self.assertEqual(matches[0][1],"clear_history")

    def testExpandSource(self):
        """Test expansion of a shell command and a file name"""
        with open('.tmp.in.source', 'w') as out:
           print('units real', file=out)
        out.close()
        matches = re.findall(shell_prompt_re, self.InputRunner(b'sour\t.tmp.in.sou\t\n'), re.MULTILINE)
        os.remove('.tmp.in.source')
        self.assertEqual(matches[0][1],"source")
        self.assertEqual(matches[0][2],".tmp.in.source")

    def testHistory(self):
        """Test history expansion"""
        out = self.InputRunner(b'clear_history\nunits real\ndimension 2\n!!:p\n!-3:p\n!dim:p\n!uni:p\nprint !!:$\nprint !dim:1\n')
        idx = 0
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
