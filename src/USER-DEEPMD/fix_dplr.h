#ifdef FIX_CLASS

FixStyle(dplr,FixDPLR)

#else

#ifndef LMP_FIX_DPLR_H
#define LMP_FIX_DPLR_H

#include <stdio.h>
#include "fix.h"
#include "pair_deepmd.h"
#include "deepmd/DeepTensor.h"
#include "deepmd/DataModifier.h"

#ifdef HIGH_PREC
#define FLOAT_PREC double
#else
#define FLOAT_PREC float
#endif

namespace LAMMPS_NS {
  class FixDPLR : public Fix {
public:
    FixDPLR(class LAMMPS *, int, char **);
    virtual ~FixDPLR() {};
    int setmask();
    void init();
    void setup(int);
    void post_integrate();
    void pre_force(int);
    void post_force(int);
    int pack_reverse_comm(int, int, double *);
    void unpack_reverse_comm(int, int *, double *);
    double compute_scalar(void);
    double compute_vector(int);
private:
    PairDeepMD * pair_deepmd;
    deepmd::DeepTensor dpt;
    deepmd::DipoleChargeModifier dtm;
    std::string model;
    int ntypes;
    std::vector<int > sel_type;
    std::vector<int > dpl_type;
    std::vector<int > bond_type;
    std::map<int,int > type_asso;
    std::map<int,int > bk_type_asso;
    std::vector<FLOAT_PREC> dipole_recd;
    std::vector<double> dfcorr_buff;
    std::vector<double> efield;
    std::vector<double> efield_fsum, efield_fsum_all;
    int efield_force_flag;
    void get_valid_pairs(std::vector<std::pair<int,int> >& pairs);
  };
}

#endif // LMP_FIX_DPLR_H
#endif // FIX_CLASS
