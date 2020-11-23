
import sys,os,unittest
from lammps import lammps, LMP_VAR_ATOM

class PythonCommand(unittest.TestCase):

    def setUp(self):
        machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.lmp=lammps(name=machine,
                        cmdargs=['-nocite',
                                 '-log','none',
                                 '-echo','screen',
                                 '-var','zpos','1.5',
                                 '-var','x','2'])
        # create demo input strings and files
        # a few commands to set up a box with a single atom
        self.demo_input="""
region       box block 0 $x 0 2 0 2
create_box 1 box
create_atoms 1 single 1.0 1.0 ${zpos}
"""
        # another command to add an atom and use a continuation line
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

    # clean up temporary files
    def tearDown(self):
        if os.path.exists(self.demo_file):
            os.remove(self.demo_file)
        if os.path.exists(self.cont_file):
            os.remove(self.cont_file)

    ##############################
    def testFile(self):
        """Test reading commands from a file"""
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,0)
        self.lmp.file(self.demo_file)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,1)
        self.lmp.file(self.cont_file)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,2)

    def testNoFile(self):
        """Test (not) reading commands from no file"""
        self.lmp.file(None)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,0)

    def testCommand(self):
        """Test executing individual commands"""
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,0)
        cmds = self.demo_input.splitlines()
        for cmd in cmds:
            self.lmp.command(cmd)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,1)

    def testCommandsList(self):
        """Test executing commands from list of strings"""
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,0)
        cmds = self.demo_input.splitlines()+self.cont_input.splitlines()
        self.lmp.commands_list(cmds)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,2)

    def testCommandsString(self):
        """Test executing block of commands from string"""
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,0)
        self.lmp.commands_string(self.demo_input+self.cont_input)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,2)

    def testNeighborList(self):
        self.lmp.command("units lj")
        self.lmp.command("atom_style atomic")
        self.lmp.command("atom_modify map array")
        self.lmp.command("boundary f f f")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")

        x = [
          1.0, 1.0, 1.0,
          1.0, 1.0, 1.5
        ]

        types = [1, 1]

        self.assertEqual(self.lmp.create_atoms(2, id=None, type=types, x=x), 2)
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 2)

        self.lmp.command("mass 1 1.0")
        self.lmp.command("velocity all create 3.0 87287")
        self.lmp.command("pair_style lj/cut 2.5")
        self.lmp.command("pair_coeff 1 1 1.0 1.0 2.5")
        self.lmp.command("neighbor 0.1 bin")
        self.lmp.command("neigh_modify every 20 delay 0 check no")

        self.lmp.command("run 0")

        self.assertEqual(self.lmp.find_pair_neighlist("lj/cut"), 0)
        nlist = self.lmp.get_neighlist(0)
        self.assertEqual(len(nlist), 2)
        atom_i, numneigh_i, neighbors_i = nlist[0]
        atom_j, numneigh_j, _ = nlist[1]

        self.assertEqual(atom_i, 0)
        self.assertEqual(atom_j, 1)

        self.assertEqual(numneigh_i, 1)
        self.assertEqual(numneigh_j, 0)

        self.assertEqual(1, neighbors_i[0])

    def test_extract_box_non_periodic(self):
        self.lmp.command("boundary f f f")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")

        boxlo, boxhi, xy, yz, xz, periodicity, box_change = self.lmp.extract_box()

        self.assertEqual(boxlo, [0.0, 0.0, 0.0])
        self.assertEqual(boxhi, [2.0, 2.0, 2.0])
        self.assertEqual(xy, 0.0)
        self.assertEqual(yz, 0.0)
        self.assertEqual(xz, 0.0)
        self.assertEqual(periodicity, [0, 0, 0])
        self.assertEqual(box_change, 0)

    def test_extract_box_periodic(self):
        self.lmp.command("boundary p p p")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")

        boxlo, boxhi, xy, yz, xz, periodicity, box_change = self.lmp.extract_box()

        self.assertEqual(boxlo, [0.0, 0.0, 0.0])
        self.assertEqual(boxhi, [2.0, 2.0, 2.0])
        self.assertEqual(xy, 0.0)
        self.assertEqual(yz, 0.0)
        self.assertEqual(xz, 0.0)
        self.assertEqual(periodicity, [1, 1, 1])
        self.assertEqual(box_change, 0)

    def test_extract_box_triclinic(self):
        self.lmp.command("boundary p p p")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")
        self.lmp.command("change_box all triclinic")
        self.lmp.command("change_box all xy final 0.1 yz final 0.2 xz final 0.3")

        boxlo, boxhi, xy, yz, xz, periodicity, box_change = self.lmp.extract_box()

        self.assertEqual(boxlo, [0.0, 0.0, 0.0])
        self.assertEqual(boxhi, [2.0, 2.0, 2.0])
        self.assertEqual(xy, 0.1)
        self.assertEqual(yz, 0.2)
        self.assertEqual(xz, 0.3)
        self.assertEqual(periodicity, [1, 1, 1])
        self.assertEqual(box_change, 0)

    def test_reset_box(self):
        self.lmp.command("boundary p p p")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")
        self.lmp.command("change_box all triclinic")
        self.lmp.command("change_box all xy final 0.1 yz final 0.2 xz final 0.3")
        self.lmp.reset_box([0,0,0], [1,1,1], 0, 0, 0)

        boxlo, boxhi, xy, yz, xz, periodicity, box_change = self.lmp.extract_box()

        self.assertEqual(boxlo, [0.0, 0.0, 0.0])
        self.assertEqual(boxhi, [1.0, 1.0, 1.0])
        self.assertEqual(xy, 0)
        self.assertEqual(yz, 0)
        self.assertEqual(xz, 0)
        self.assertEqual(periodicity, [1, 1, 1])
        self.assertEqual(box_change, 0)

    def test_extract_variable_equalstyle(self):
        self.lmp.command("variable a equal 100")
        a = self.lmp.extract_variable("a")
        self.assertEqual(a, 100)

        self.lmp.command("variable a equal 3.14")
        a = self.lmp.extract_variable("a")
        self.assertEqual(a, 3.14)

    def test_extract_variable_atomstyle(self):
        self.lmp.command("units lj")
        self.lmp.command("atom_style atomic")
        self.lmp.command("atom_modify map array")
        self.lmp.command("boundary f f f")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")

        x = [
          1.0, 1.0, 1.0,
          1.0, 1.0, 1.5
        ]

        types = [1, 1]

        self.assertEqual(self.lmp.create_atoms(2, id=None, type=types, x=x), 2)
        self.lmp.command("variable a atom x*x+y*y+z*z")
        a = self.lmp.extract_variable("a", "all", LMP_VAR_ATOM)
        self.assertEqual(a[0], x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
        self.assertEqual(a[1], x[3]*x[3]+x[4]*x[4]+x[5]*x[5])

    def test_get_thermo(self):
        self.lmp.command("units lj")
        self.lmp.command("atom_style atomic")
        self.lmp.command("atom_modify map array")
        self.lmp.command("boundary f f f")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")

        x = [
          1.0, 1.0, 1.0,
          1.0, 1.0, 1.5
        ]

        types = [1, 1]
        self.lmp.create_atoms(2, id=None, type=types, x=x)

        state = {
            "step": 0,
            "elapsed" : 0.0,
            "elaplong": 0,
            "dt" : 0.005,
            "time" : 0.0,
            "atoms" : 2.0,
            "temp" : 0,
            "press" : 0,
            "pe" : 0.0,
            "ke" : 0.0,
            "etotal" : 0.0,
            "enthalpy" : 0.0,
            "vol" : 8.0,
            "lx" : 2.0,
            "ly" : 2.0,
            "lz" : 2.0,
            "xlo" : 0,
            "xhi" : 2.0,
            "ylo" : 0,
            "yhi" : 2.0,
            "zlo" : 0,
            "zhi" : 2.0
        }

        for key, value in state.items():
            result = self.lmp.get_thermo(key)
            self.assertEqual(value, result, key)

    def test_extract_global_double(self):
        self.lmp.command("region box block -1 1 -2 2 -3 3")
        self.lmp.command("create_box 1 box")
        self.assertEqual(self.lmp.extract_global("boxxlo"), -1.0)
        self.assertEqual(self.lmp.extract_global("boxxhi"), 1.0)
        self.assertEqual(self.lmp.extract_global("boxylo"), -2.0)
        self.assertEqual(self.lmp.extract_global("boxyhi"), 2.0)
        self.assertEqual(self.lmp.extract_global("boxzlo"), -3.0)
        self.assertEqual(self.lmp.extract_global("boxzhi"), 3.0)

##############################
if __name__ == "__main__":
    unittest.main()
