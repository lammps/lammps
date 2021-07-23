import sys,os,unittest
from ctypes import *
from lammps import lammps, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR

try:
    import numpy
    NUMPY_INSTALLED = True
except ImportError:
    NUMPY_INSTALLED = False

# add timestep dependent force
def callback_one(lmp, ntimestep, nlocal, tag, x, f):
    lmp.fix_external_set_virial_global("ext",[1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    for i in range(nlocal):
        f[i][0] = float(ntimestep)
        f[i][1] = float(ntimestep)
        f[i][2] = float(ntimestep)
    if ntimestep < 10:
        lmp.fix_external_set_energy_global("ext", 0.5)
        lmp.fix_external_set_vector("ext", 1, ntimestep)
        lmp.fix_external_set_vector("ext", 3, 1.0)
        lmp.fix_external_set_vector("ext", 4, -0.25)
    else:
        lmp.fix_external_set_energy_global("ext", 1.0)
        lmp.fix_external_set_vector("ext", 2, ntimestep)
        lmp.fix_external_set_vector("ext", 5, -1.0)
        lmp.fix_external_set_vector("ext", 6, 0.25)

    eatom = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7]
    vatom = [ [0.1,0.0,0.0,0.0,0.0,0.0],
              [0.0,0.2,0.0,0.0,0.0,0.0],
              [0.0,0.0,0.3,0.0,0.0,0.0],
              [0.0,0.0,0.0,0.4,0.0,0.0],
              [0.0,0.0,0.0,0.0,0.5,0.0],
              [0.0,0.0,0.0,0.0,0.0,0.6],
              [0.0,0.0,0.0,0.0,-7.0,0.0],
              [0.0,-8.0,0.0,0.0,0.0,0.0] ]
    if ntimestep < 5:
        lmp.fix_external_set_energy_peratom("ext",eatom)
        lmp.fix_external_set_virial_peratom("ext",vatom)
    else:
        import numpy as np
        eng = np.array(eatom)
        vir = np.array(vatom)

        lmp.numpy.fix_external_set_energy_peratom("ext",eng)
        lmp.numpy.fix_external_set_virial_peratom("ext",vir)

    # ------------------------------------------------------------------------

class PythonExternal(unittest.TestCase):
    @unittest.skipIf(not NUMPY_INSTALLED, "NumPy is not available")
    def testExternalCallback(self):
        """Test fix external from Python with pf/callback"""

        machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        lmp=lammps(name=machine, cmdargs=['-nocite', '-log','none', '-echo', 'screen'])

        # a few commands to set up simple system
        basic_system="""lattice sc 1.0
                        region box block -1 1 -1 1 -1 1
                        create_box 1 box
                        create_atoms 1 box
                        mass 1 1.0
                        pair_style zero 0.1
                        pair_coeff 1 1
                        velocity all set 0.1 0.0 -0.1
                        thermo_style custom step temp pe ke etotal press
                        thermo 5
                        fix 1 all nve
                        fix ext all external pf/callback 5 1
                        compute eatm all pe/atom fix
                        compute vatm all stress/atom NULL fix
                        compute sum all reduce sum c_eatm c_vatm[*]
                        thermo_style custom step temp pe ke etotal press c_sum[*]
                        fix_modify ext energy yes virial yes
"""
        lmp.commands_string(basic_system)
        lmp.fix_external_set_vector_length("ext",6);
        lmp.set_fix_external_callback("ext",callback_one,lmp)

        # check setting per-atom data with python lists
        lmp.command("run 0 post no")
        reduce = lmp.extract_compute("sum", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR)
        self.assertAlmostEqual(reduce[0],2.8,14)
        self.assertAlmostEqual(reduce[1],-0.1,14)
        self.assertAlmostEqual(reduce[2],7.8,14)
        self.assertAlmostEqual(reduce[3],-0.3,14)
        self.assertAlmostEqual(reduce[4],-0.4,14)
        self.assertAlmostEqual(reduce[5],6.5,14)
        self.assertAlmostEqual(reduce[6],-0.6,14)

        lmp.command("run 10 post no")
        self.assertAlmostEqual(lmp.get_thermo("temp"),1.0/30.0,14)
        self.assertAlmostEqual(lmp.get_thermo("pe"),1.0/8.0,14)
        self.assertAlmostEqual(lmp.get_thermo("press"),0.15416666666666667,14)
        # check setting per-atom data numpy arrays
        reduce = lmp.extract_compute("sum", LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR)
        self.assertAlmostEqual(reduce[0],2.8,14)
        self.assertAlmostEqual(reduce[1],-0.1,14)
        self.assertAlmostEqual(reduce[2],7.8,14)
        self.assertAlmostEqual(reduce[3],-0.3,14)
        self.assertAlmostEqual(reduce[4],-0.4,14)
        self.assertAlmostEqual(reduce[5],6.5,14)
        self.assertAlmostEqual(reduce[6],-0.6,14)
        val = 0.0
        for i in range(0,6):
            val += lmp.extract_fix("ext",LMP_STYLE_GLOBAL,LMP_TYPE_VECTOR,nrow=i)
        self.assertAlmostEqual(val,15.0,14)

    def testExternalArray(self):
        """Test fix external from Python with pf/array"""

        machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        lmp=lammps(name=machine, cmdargs=['-nocite', '-log','none', '-echo', 'screen'])

        # a few commands to set up simple system
        basic_system="""lattice sc 1.0
                        region box block -1 1 -1 1 -1 1
                        create_box 1 box
                        create_atoms 1 box
                        mass 1 1.0
                        pair_style zero 0.1
                        pair_coeff 1 1
                        velocity all set 0.1 0.0 -0.1
                        fix 1 all nve
                        fix ext all external pf/array 1
                        fix_modify ext energy yes virial yes
                        thermo_style custom step temp pe ke press
                        thermo 5
"""
        lmp.commands_string(basic_system)
        force = lmp.fix_external_get_force("ext");
        nlocal = lmp.extract_setting("nlocal");
        for i in range(nlocal):
            force[i][0] = 0.0
            force[i][1] = 0.0
            force[i][2] = 0.0
        lmp.fix_external_set_energy_global("ext", 0.5)
        lmp.fix_external_set_virial_global("ext",[0.5, 0.5, 0.5, 0.0, 0.0, 0.0])

        lmp.command("run 5 post no")
        self.assertAlmostEqual(lmp.get_thermo("temp"),4.0/525.0,14)
        self.assertAlmostEqual(lmp.get_thermo("pe"),1.0/16.0,14)
        self.assertAlmostEqual(lmp.get_thermo("press"),0.06916666666666667,14)
        if NUMPY_INSTALLED:
            npforce = lmp.numpy.fix_external_get_force("ext")
            self.assertEqual(len(npforce),8)
            self.assertEqual(len(npforce[0]),3)
            self.assertEqual(npforce[1][1],0.0)

        force = lmp.fix_external_get_force("ext");
        nlocal = lmp.extract_setting("nlocal");
        for i in range(nlocal):
            force[i][0] = 6.0
            force[i][1] = 6.0
            force[i][2] = 6.0
        lmp.fix_external_set_energy_global("ext", 1.0)
        lmp.fix_external_set_virial_global("ext",[1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        lmp.command("run 5 post no")
        self.assertAlmostEqual(lmp.get_thermo("temp"),1.0/30.0,14)
        self.assertAlmostEqual(lmp.get_thermo("pe"),1.0/8.0,14)
        self.assertAlmostEqual(lmp.get_thermo("press"),0.15416666666666667,14)
        if NUMPY_INSTALLED:
            npforce = lmp.numpy.fix_external_get_force("ext")
            self.assertEqual(npforce[0][0],6.0)
            self.assertEqual(npforce[3][1],6.0)
            self.assertEqual(npforce[7][2],6.0)

##############################
if __name__ == "__main__":
    unittest.main()
