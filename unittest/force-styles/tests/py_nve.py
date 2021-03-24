from __future__ import print_function
from lammps import lammps

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
        nlocal = self.lmp.extract_global("nlocal")
        mass = self.lmp.extract_atom("mass")
        atype = self.lmp.extract_atom("type")
        x = self.lmp.extract_atom("x")
        v = self.lmp.extract_atom("v")
        f = self.lmp.extract_atom("f")

        for i in range(nlocal):
            dtfm = self.dtf / mass[int(atype[i])]
            v[i][0] += dtfm * f[i][0]
            v[i][1] += dtfm * f[i][1]
            v[i][2] += dtfm * f[i][2]
            x[i][0] += self.dtv * v[i][0]
            x[i][1] += self.dtv * v[i][1]
            x[i][2] += self.dtv * v[i][2]

    def final_integrate(self):
        nlocal = self.lmp.extract_global("nlocal")
        mass = self.lmp.extract_atom("mass")
        atype = self.lmp.extract_atom("type")
        v = self.lmp.extract_atom("v")
        f = self.lmp.extract_atom("f")

        for i in range(nlocal):
            dtfm = self.dtf / mass[int(atype[i])]
            v[i][0] += dtfm * f[i][0]
            v[i][1] += dtfm * f[i][1]
            v[i][2] += dtfm * f[i][2]


