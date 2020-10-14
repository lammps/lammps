
import os, re, subprocess, unittest

# enable test mode
os.putenv('LAMMPS_SHELL_TESTING','1')

shell_prompt_re = r"(.*LAMMPS Shell> ([a-z0-9_]+) *([a-z0-9_\.]+)?.*)+"

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

        #print(outs.decode())
        return re.findall(shell_prompt_re, outs.decode('UTF-8'), re.MULTILINE)[0]

    def testExpandClear(self):
        """Test expansion of a shell specific command"""
        self.assertEqual(self.InputRunner(b'cle\t\n')[1],"clear")

    def testExpandSource(self):
        """Test expansion of a shell command and a file name"""
        with open('.tmp.in.source', 'w') as out:
           print('units real', file=out)
        out.close()
        matches=self.InputRunner(b'sour\t.tmp.in.sou\t\n')
        os.remove('.tmp.in.source')
        self.assertEqual(matches[1],"source")
        self.assertEqual(matches[2],".tmp.in.source")


###########################
if __name__ == "__main__":
    unittest.main()
