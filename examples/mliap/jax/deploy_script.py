import lammps
import lammps.mliap

#from lammps.mliap.mliap_unified_lj import MLIAPUnifiedLJ
from mliap_unified_jax import MLIAPUnifiedJAX

def create_pickle():
    unified = MLIAPUnifiedJAX(["Ar"])
    unified.pickle('mliap_unified_jax_Ar.pkl')

create_pickle()