import lammps
import lammps.mliap

from lammps.mliap.mliap_unified_lj import MLIAPUnifiedLJ


if __name__ == '__main__':
    unified = MLIAPUnifiedLJ()
    unified.pickle('mliap_unified_lj_Ar.pkl')
