
#ifndef SHIPs_RADIAL_FUNCTIONS_H
#define SHIPs_RADIAL_FUNCTIONS_H

#include "ace_arraynd.h"
#include "ace_types.h"
#include "ace_radial.h"


class SHIPsRadialFunctions : public AbstractRadialBasis {
public:
    // transform parameters
    int p = 0;
    DOUBLE_TYPE r0 = 0.0;

    // cutoff parameters
    DOUBLE_TYPE rcut = 0.0;
    DOUBLE_TYPE xl = 0.0;
    DOUBLE_TYPE xr = 0.0;
    int pl = 0;
    int pr = 0;

    // basis size
    size_t maxn = 0;

    // recursion parameters
    Array1D<DOUBLE_TYPE> A = Array1D<DOUBLE_TYPE>("SHIPs radial basis: A");
    Array1D<DOUBLE_TYPE> B = Array1D<DOUBLE_TYPE>("SHIPs radial basis: B");
    Array1D<DOUBLE_TYPE> C = Array1D<DOUBLE_TYPE>("SHIPs radial basis: C");

    // temporary storage for evaluating the basis
    Array1D<DOUBLE_TYPE> P = Array1D<DOUBLE_TYPE>("SHIPs radial basis: P");
    Array1D<DOUBLE_TYPE> dP_dr = Array1D<DOUBLE_TYPE>("SHIPs radial basis: dP");
//////////////////////////////////

    SHIPsRadialFunctions() = default;

    ~SHIPsRadialFunctions() override = default;

// distance transform
    void transform(const DOUBLE_TYPE r, DOUBLE_TYPE &x_out, DOUBLE_TYPE &dx_out) const;

    // cutoff function
    void fcut(const DOUBLE_TYPE x, DOUBLE_TYPE &f_out, DOUBLE_TYPE &df_out) const;


    void fread(FILE *fptr);

    void load(string fname);

    void _init(DOUBLE_TYPE r0, int p, DOUBLE_TYPE rcut,
               DOUBLE_TYPE xl, DOUBLE_TYPE xr,
               int pl, int pr, size_t maxn);

    void calcP(DOUBLE_TYPE r, size_t maxn, SPECIES_TYPE z1, SPECIES_TYPE z2);


    size_t get_maxn();

    //TODO
    void init(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, DOUBLE_TYPE deltaSplineBins, SPECIES_TYPE nelements,
              DOUBLE_TYPE cutoff,
              string radbasename) override;

    //TODO
    void
    evaluate(DOUBLE_TYPE r, NS_TYPE nradbase_c, NS_TYPE nradial_c, SPECIES_TYPE mu_i, SPECIES_TYPE mu_j,
             bool calc_second_derivatives = false) override;

    //TODO
    void setuplookupRadspline() override;


    SHIPsRadialFunctions *clone() const override {
        return new SHIPsRadialFunctions(*this);
    };

    /**
     * Helper method, that populate `fr` and `dfr` 2D-arrays (n,l) with P(n), dP_dr  for given coordinate r
     * @param r
     * @param maxn
     * @param z1
     * @param z2
     */
    void fill_Rnl(DOUBLE_TYPE r, NS_TYPE maxn, SPECIES_TYPE z1, SPECIES_TYPE z2);

    void fill_gk(DOUBLE_TYPE r, NS_TYPE maxn, SPECIES_TYPE z1, SPECIES_TYPE z2);
};


#endif
