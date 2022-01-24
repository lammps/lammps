from __future__ import print_function
from lammps import lammps
import numpy as np

class LAMMPSFix(object):
    def __init__(self, ptr, group_name="all"):
        self.lmp = lammps(ptr=ptr)
        self.group_name = group_name

class LAMMPSFixMove(LAMMPSFix):
    def __init__(self, ptr, group_name="all"):
        super(LAMMPSFixMove, self).__init__(ptr, group_name)

    def init(self):
        pass

    def initial_integrate(self, vflag):
        pass

    def final_integrate(self):
        pass

    def initial_integrate_respa(self, vflag, ilevel, iloop):
        pass

    def final_integrate_respa(self, ilevel, iloop):
        pass

    def reset_dt(self):
        pass


class NVE(LAMMPSFixMove):
    """ Python implementation of fix/nve """
    def __init__(self, ptr, group_name="all"):
        super(NVE, self).__init__(ptr)
        assert(self.group_name == "all")

    def init(self):
        dt = self.lmp.extract_global("dt")
        ftm2v = self.lmp.extract_global("ftm2v")
        self.ntypes = self.lmp.extract_global("ntypes")
        self.dtv = dt
        self.dtf = 0.5 * dt * ftm2v

    def initial_integrate(self, vflag):
        mass = self.lmp.numpy.extract_atom("mass")
        atype = self.lmp.numpy.extract_atom("type")
        x = self.lmp.numpy.extract_atom("x")
        v = self.lmp.numpy.extract_atom("v")
        f = self.lmp.numpy.extract_atom("f")

        for i in range(x.shape[0]):
            dtfm = self.dtf / mass[int(atype[i])]
            v[i,:]+= dtfm * f[i,:]
            x[i,:] += self.dtv * v[i,:]

    def final_integrate(self):
        mass = self.lmp.numpy.extract_atom("mass")
        atype = self.lmp.numpy.extract_atom("type")
        v = self.lmp.numpy.extract_atom("v")
        f = self.lmp.numpy.extract_atom("f")

        for i in range(v.shape[0]):
            dtfm = self.dtf / mass[int(atype[i])]
            v[i,:] += dtfm * f[i,:]


class NVE_Opt(LAMMPSFixMove):
    """ Performance-optimized Python implementation of fix/nve """
    def __init__(self, ptr, group_name="all"):
        super(NVE_Opt, self).__init__(ptr)
        assert(self.group_name == "all")

    def init(self):
        dt = self.lmp.extract_global("dt")
        ftm2v = self.lmp.extract_global("ftm2v")
        self.ntypes = self.lmp.extract_global("ntypes")
        self.dtv = dt
        self.dtf = 0.5 * dt * ftm2v

    def initial_integrate(self, vflag):        
        mass = self.lmp.numpy.extract_atom("mass")
        atype = self.lmp.numpy.extract_atom("type")
        x = self.lmp.numpy.extract_atom("x")
        v = self.lmp.numpy.extract_atom("v")
        f = self.lmp.numpy.extract_atom("f")
        dtf = self.dtf
        dtv = self.dtv

        dtfm = dtf / np.take(mass, atype)

        for d in range(x.shape[1]):
            v[:,d] += dtfm * f[:,d]
            x[:,d] += dtv * v[:,d]

    def final_integrate(self):
        mass = self.lmp.numpy.extract_atom("mass")
        atype = self.lmp.numpy.extract_atom("type")
        v = self.lmp.numpy.extract_atom("v")
        f = self.lmp.numpy.extract_atom("f")
        dtf = self.dtf

        dtfm = dtf / np.take(mass, atype)

        for d in range(v.shape[1]):
            v[:,d] += dtfm * f[:,d]
