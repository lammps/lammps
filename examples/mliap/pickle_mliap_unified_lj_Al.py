import pickle
import lammps
import lammps.mliap

from lammps.mliap.mliap_unified_lj import MLIAPUnifiedLJ


if __name__ == '__main__':
    unified = MLIAPUnifiedLJ()
    unified.element_types = ["Al"]
    unified.ndescriptors = 1
    unified.nparams = 3
    # Mimicking the LJ pair-style:
    # pair_style lj/cut 2.5
    # pair_coeff * * 1 1
    unified.epsilon = 1.0
    unified.sigma = 2.5
    unified.rcutfac = 5.0

    with open('mliap_unified_lj_Al.pkl', 'wb') as fp:
        pickle.dump(unified, fp)
