from __future__ import print_function

class LAMMPSLJCutPotential(object):

    def __init__(self):
        self.pmap=dict()
        # set coeffs: eps, sig, 48*eps*sig**12, 24*eps*sig**6,
        #                        4*eps*sig**12,  4*eps*sig**6
        self.coeff = {'lj'  : {'lj'  : (1.0,1.0,48.0,24.0,4.0,4.0),
                               'NULL': (0.0,1.0, 0.0, 0.0,0.0,0.0)},
                      'NULL': {'lj'  : (0.0,1.0, 0.0, 0.0,0.0,0.0),
                               'NULL': (0.0,1.0, 0.0, 0.0,0.0,0.0)}}

    def map_coeff(self,name,type):
        if self.coeff.has_key(name):
           self.pmap[type] = name
        else:
           raise Exception("cannot match atom type %s" % name)

    def compute_force(self,rsq,itype,jtype):
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj1 = self.coeff[self.pmap[itype]][self.pmap[jtype]][2]
        lj2 = self.coeff[self.pmap[itype]][self.pmap[jtype]][3]
        return (r6inv * (lj1*r6inv - lj2))

    def compute_energy(self,rsq,itype,jtype):
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj3 = self.coeff[self.pmap[itype]][self.pmap[jtype]][4]
        lj4 = self.coeff[self.pmap[itype]][self.pmap[jtype]][5]
        return (r6inv * (lj3*r6inv - lj4))

lammps_pair_style = LAMMPSLJCutPotential()

