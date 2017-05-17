from __future__ import print_function

class LJCutMelt(object):

    def __init__(self):
        self.pmap=dict()
        # set coeffs: eps, sig, 48*eps*sig**12, 24*eps*sig**6,
        #                        4*eps*sig**12,  4*eps*sig**6
        self.coeff = {'lj'  : {'lj'  : (1.0,1.0,48.0,24.0,4.0,4.0),
                               'NULL': (0.0,1.0, 0.0, 0.0,0.0,0.0)},
                      'NULL': {'lj'  : (0.0,1.0, 0.0, 0.0,0.0,0.0),
                               'NULL': (0.0,1.0, 0.0, 0.0,0.0,0.0)}}

    def map_coeff(self,name,type):
        if name in self.coeff:
           self.pmap[type] = name
        else:
           raise Exception("cannot match atom type %s" % name)

    def compute_force(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj1 = coeff[2]
        lj2 = coeff[3]
        return (r6inv * (lj1*r6inv - lj2))

    def compute_energy(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj3 = coeff[4]
        lj4 = coeff[5]
        return (r6inv * (lj3*r6inv - lj4))

class LJCutSPCE(object):

    def __init__(self):
        self.pmap=dict()
        # SPCE oxygen in real units
        eps=0.15535
        sig=3.166

        # set coeffs: eps, sig, 48*eps*sig**12, 24*eps*sig**6,
        #                        4*eps*sig**12,  4*eps*sig**6
        self.coeff = {'OW'  : {'OW'  : (1.0,1.0,
                                48.0*eps*sig**12,24.0*eps*sig**6,
                                 4.0*eps*sig**12, 4.0*eps*sig**6),
                               'NULL': (0.0,1.0, 0.0, 0.0,0.0,0.0)},
                      'NULL': {'OW'  : (0.0,1.0, 0.0, 0.0,0.0,0.0),
                               'NULL': (0.0,1.0, 0.0, 0.0,0.0,0.0)}}

    def map_coeff(self,name,type):
        if name in self.coeff:
           self.pmap[type] = name
        else:
           raise Exception("cannot match atom type %s" % name)

    def compute_force(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj1 = coeff[2]
        lj2 = coeff[3]
        return (r6inv * (lj1*r6inv - lj2))

    def compute_energy(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj3 = coeff[4]
        lj4 = coeff[5]
        return (r6inv * (lj3*r6inv - lj4))
