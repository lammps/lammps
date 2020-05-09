
import sys,os,unittest
from lammps import lammps

class PythonCommand(unittest.TestCase):

    def setUp(self):
        machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.lmp=lammps(name=machine,
                        cmdargs=['-nocite','-log','none','-echo','screen'])
        # create demo input strings and files
        self.demo_input="""
region       box block 0 2 0 2 0 2
create_box 1 box
create_atoms 1 single 1.0 1.0 1.0
"""
        self.cont_input="""
create_atoms 1 single &
            0.2 0.1 0.1
"""
        self.demo_file='in.test'
        with open(self.demo_file,'w') as f:
            f.write(self.demo_input)
        self.cont_file='in.cont'
        with open(self.cont_file,'w') as f:
            f.write(self.cont_input)

    def tearDown(self):
        if os.path.exists(self.demo_file):
            os.remove(self.demo_file)
        if os.path.exists(self.cont_file):
            os.remove(self.cont_file)

    ##############################
    def testFile(self):
        """Test reading commands from a file"""
        self.lmp.file(self.demo_file)
        self.lmp.file(self.cont_file)

    def testNoFile(self):
        """Test (not) reading commands from no file"""
        self.lmp.file(None)

    def testCommand(self):
        """Test executing individual commands"""
        cmds = self.demo_input.splitlines()
        for cmd in cmds:
            self.lmp.command(cmd)

    def testCommandsList(self):
        """Test executing list of commands"""
        cmds = self.demo_input.splitlines()+self.cont_input.splitlines()
        self.lmp.commands_list(cmds)

    def testCommandsString(self):
        """Test executing block of commands"""
        self.lmp.commands_string(self.demo_input+self.cont_input)

##############################
if __name__ == "__main__":
    unittest.main()
