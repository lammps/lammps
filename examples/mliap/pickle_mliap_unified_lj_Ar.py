import lammps
import lammps.mliap

from lammps.mliap.mliap_unified_lj import MLIAPUnifiedLJ


if __name__ == '__main__':
    unified = MLIAPUnifiedLJ()
    unified.element_types = ["Ar"]
    unified.ndescriptors = 1
    unified.nparams = 3
    # Mimicking the LJ pair-style:
    # pair_style lj/cut 2.5
    # pair_coeff * * 1 1
    unified.epsilon = 1.0
    unified.sigma = 1.0
    unified.rcutfac = 1.25

    unified.pickle('mliap_unified_lj_Ar.pkl')
