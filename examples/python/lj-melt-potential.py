from __future__ import print_function

class LAMMPSLJCutPotential(object):

    def __init__(self):
        self.pmap=dict()
        # coefficients: epsilon, sigma
        self.coeff = {'lj'  : {'lj'  : (1.0,1.0),
                               'NULL': (0.0,1.0)},
                      'NULL': {'lj'  : (0.0,1.0),
                               'NULL': (0.0,1.0)}}

    def map_coeff(self,type,name):
        if self.coeff.has_key(name):
           print("map type %d to name %s" % (type,name))
           self.pmap[type] = name
        else:
           print("cannot match atom type",name)

    def compute_force(self,r,itype,jtype,factor_lj):
        return 0.0

    def compute_energy(self,r,itype,jtype,factor_lj):
        return 0.0

lammps_pair_style = LAMMPSLJCutPotential()
print ("lj-melt potential file loaded")

