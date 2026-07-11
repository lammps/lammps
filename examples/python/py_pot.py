from __future__ import print_function
import math

class LAMMPSPairPotential(object):
    def __init__(self):
        self.pmap=dict()
        self.units='lj'
    def map_coeff(self,name,ltype):
        self.pmap[ltype]=name
    def check_units(self,units):
        if (units != self.units):
           raise Exception("Conflicting units: %s vs. %s" % (self.units,units))

class Harmonic(LAMMPSPairPotential):
    def __init__(self):
        super(Harmonic,self).__init__()
        self.units = 'real'
        # set coeffs: K, r0
        self.coeff = {'A'  : {'A'  : (0.2,9.0),
                              'B'  : (math.sqrt(0.2*0.4),9.0)},
                      'B'  : {'A'  : (math.sqrt(0.2*0.4),9.0),
                              'B'  : (0.4,9.0)}}

    def compute_force(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r = math.sqrt(rsq)
        delta = coeff[1]-r
        if (r <= coeff[1]):
          return 2.0*delta*coeff[0]/r
        else:
          return 0.0

    def compute_energy(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r = math.sqrt(rsq)
        delta = coeff[1]-r
        if (r <= coeff[1]):
          return delta*delta*coeff[0]
        else:
          return 0.0

class LJCutMelt(LAMMPSPairPotential):
    def __init__(self):
        super(LJCutMelt,self).__init__()
        # set coeffs: 48*eps*sig**12, 24*eps*sig**6,
        #              4*eps*sig**12,  4*eps*sig**6
        self.units = 'lj'
        self.coeff = {'lj'  : {'lj'  : (48.0,24.0,4.0,4.0)}}

    def compute_force(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj1 = coeff[0]
        lj2 = coeff[1]
        return (r6inv * (lj1*r6inv - lj2))*r2inv

    def compute_energy(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj3 = coeff[2]
        lj4 = coeff[3]
        return (r6inv * (lj3*r6inv - lj4))


class LJCutSPCE(LAMMPSPairPotential):
    def __init__(self):
        super(LJCutSPCE,self).__init__()
        self.units='real'
        # SPCE oxygen LJ parameters in real units
        eps=0.15535
        sig=3.166
        self.coeff = {'OW'  : {'OW'  : (48.0*eps*sig**12,24.0*eps*sig**6,
                                         4.0*eps*sig**12, 4.0*eps*sig**6),
                               'HW'  : (0.0,0.0, 0.0,0.0)},
                      'HW'  : {'OW'  : (0.0,0.0, 0.0,0.0),
                               'HW'  : (0.0,0.0, 0.0,0.0)}}

    def compute_force(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj1 = coeff[0]
        lj2 = coeff[1]
        return (r6inv * (lj1*r6inv - lj2))*r2inv

    def compute_energy(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj3 = coeff[2]
        lj4 = coeff[3]
        return (r6inv * (lj3*r6inv - lj4))
