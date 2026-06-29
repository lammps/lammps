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
        super(NVE, self).__init__(ptr, group_name)
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
        nlocal = self.lmp.extract_setting("nlocal")

        for i in range(nlocal):
            dtfm = self.dtf / mass[int(atype[i])]
            v[i,:]+= dtfm * f[i,:]
            x[i,:] += self.dtv * v[i,:]

    def final_integrate(self):
        mass = self.lmp.numpy.extract_atom("mass")
        atype = self.lmp.numpy.extract_atom("type")
        v = self.lmp.numpy.extract_atom("v")
        f = self.lmp.numpy.extract_atom("f")
        nlocal = self.lmp.extract_setting("nlocal")

        for i in range(nlocal):
            dtfm = self.dtf / mass[int(atype[i])]
            v[i,:] += dtfm * f[i,:]


class NVE_Opt(LAMMPSFixMove):
    """ Performance-optimized Python implementation of fix/nve """
    def __init__(self, ptr, group_name="all"):
        super(NVE_Opt, self).__init__(ptr, group_name)
        assert(self.group_name == "all")

    def init(self):
        dt = self.lmp.extract_global("dt")
        ftm2v = self.lmp.extract_global("ftm2v")
        self.ntypes = self.lmp.extract_global("ntypes")
        self.dtv = dt
        self.dtf = 0.5 * dt * ftm2v

    def initial_integrate(self, vflag):
        nlocal = self.lmp.extract_setting("nlocal")
        mass = self.lmp.numpy.extract_atom("mass")
        atype = self.lmp.numpy.extract_atom("type")
        x = self.lmp.numpy.extract_atom("x")[:nlocal,:]
        v = self.lmp.numpy.extract_atom("v")[:nlocal,:]
        f = self.lmp.numpy.extract_atom("f")[:nlocal,:]
        dtf = self.dtf
        dtv = self.dtv

        dtfm = dtf / np.take(mass, atype[:nlocal])

        for d in range(x.shape[1]):
            v[:,d] += dtfm * f[:,d]
            x[:,d] += dtv * v[:,d]

    def final_integrate(self):
        nlocal = self.lmp.extract_setting("nlocal")
        mass = self.lmp.numpy.extract_atom("mass")
        atype = self.lmp.numpy.extract_atom("type")
        v = self.lmp.numpy.extract_atom("v")[:nlocal,:]
        f = self.lmp.numpy.extract_atom("f")[:nlocal,:]

        dtf = self.dtf
        dtfm = dtf / np.take(mass, atype[:nlocal])

        for d in range(v.shape[1]):
            v[:,d] += dtfm * f[:,d]

class NVE_Group(LAMMPSFixMove):
    """ Python implementation of fix/nve with group"""
    def __init__(self, ptr, group_name="half"):
        super(NVE_Group, self).__init__(ptr, group_name)
        assert(self.group_name == "half")

    def init(self):
        dt = self.lmp.extract_global("dt")
        ftm2v = self.lmp.extract_global("ftm2v")
        self.ntypes = self.lmp.extract_global("ntypes")
        self.dtv = dt
        self.dtf = 0.5 * dt * ftm2v
        group_index = self.lmp.available_ids("group").index(self.group_name)
        self.group_mask = 1 << group_index

    def initial_integrate(self, vflag):
        mass = self.lmp.numpy.extract_atom("mass")
        atype = self.lmp.numpy.extract_atom("type")
        mask = self.lmp.numpy.extract_atom("mask")
        x = self.lmp.numpy.extract_atom("x")
        v = self.lmp.numpy.extract_atom("v")
        f = self.lmp.numpy.extract_atom("f")
        nlocal = self.lmp.extract_setting("nlocal")

        for i in range(nlocal):
            if mask[i] & self.group_mask:
                dtfm = self.dtf / mass[int(atype[i])]
                v[i,:]+= dtfm * f[i,:]
                x[i,:] += self.dtv * v[i,:]

    def final_integrate(self):
        mass = self.lmp.numpy.extract_atom("mass")
        mask = self.lmp.numpy.extract_atom("mask")
        atype = self.lmp.numpy.extract_atom("type")
        v = self.lmp.numpy.extract_atom("v")
        f = self.lmp.numpy.extract_atom("f")
        nlocal = self.lmp.extract_setting("nlocal")

        for i in range(nlocal):
            if mask[i] & self.group_mask:
                dtfm = self.dtf / mass[int(atype[i])]
                v[i,:] += dtfm * f[i,:]



