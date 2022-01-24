#ifndef STRESS_H
#define STRESS_H

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"
#include "ATC_TypeDefs.h"
#include "NonLinearSolver.h"

namespace ATC {
  enum ElasticityTensorType {FIRST_ELASTICITY_TENSOR=0, SECOND_ELASTICITY_TENSOR};
  /**
   * @class Stress
   * @brief Base class that defines interface for a constitutive law
   * @brief that computes stress given all field and gradient information.
   */
  class Stress
  {
    public:
      Stress()  {};
      virtual ~Stress() {};
      virtual void initialize(void){};
      //* Returns parameter values, (Nothing uses this).
      virtual void parameters(std::map<std::string,double> & /* parameters */) {}
      //* Computes stress given a displacement gradient.
      //* Units: mvv/L^3 (i.e. for units Real: g/(mol ps^2 A^2) )
      virtual void stress(const FIELD_MATS &fields,
                          const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT_VEC &stress)=0;
      //* Computes free (T>0)/potential(T=0) energy density
      //* Units: mvv/L^3 (i.e. for units Real: g/(mol ps^2 A^2) )
      virtual void elastic_energy(const FIELD_MATS &fields,
                                  const GRAD_FIELD_MATS &gradFields,
                                  DENS_MAT &energy) const;
      //* Returns the material tangent at a given deformation gradient.
      virtual void tangent(const MATRIX & /* F */, MATRIX & /* C */) const
        {throw ATC_Error("Stress::tangent: unimplemented function");}
  };


  /**
   *  @class  StressCubicElastic
   *  @brief  Class for computing stress for a cubic elastic material
   */

  class StressCubicElastic : public Stress
  {
    public:
      StressCubicElastic():c11_(0),c12_(0),c44_(0){};
      StressCubicElastic(std::fstream &matfile);
      StressCubicElastic(double c11, double c12, double c44)
        : c11_(c11), c12_(c12), c44_(c44) { }
      void stress(const FIELD_MATS &fields,
                  const GRAD_FIELD_MATS &gradFields,
                  DENS_MAT_VEC &flux);
      virtual void elastic_energy(const FIELD_MATS &fields,
                          const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT &energy) const;
      virtual void tangent(const MATRIX & /* F */, MATRIX &C) const {C=C_;}
    protected:
      double c11_, c12_, c44_;
      DENS_MAT C_;
      void set_tangent();
  };

  /**
   *  @class  StressCubicElasticDamped
   *  @brief  Class for computing stress for a cubic elastic material w/ damping
   */

  class StressCubicElasticDamped : public StressCubicElastic
  {
    public:
      StressCubicElasticDamped(std::fstream &matfile);
      StressCubicElasticDamped(double c11, double c12, double c44, double gamma)
        : StressCubicElastic(c11,c12,c44), gamma_(gamma) { }
      void stress(const FIELD_MATS &fields,
                  const GRAD_FIELD_MATS &gradFields,
                  DENS_MAT_VEC &flux);
    protected:
      double gamma_;
  };

  /**
   *  @class  StressLinearElastic
   *  @brief  Class for computing stress for a linear elastic material
   */

  class StressLinearElastic : public StressCubicElastic
  {
    public:
      StressLinearElastic(std::fstream &matfile);
      void stress(const FIELD_MATS &fields,
                  const GRAD_FIELD_MATS &gradFields,
                  DENS_MAT_VEC &flux);
    protected:
      double E_, nu_;
      double mu_, lambda_;
  };

  // forward declarations needed by StressCauchyBorn
  class CbPotential;
  class CBLattice;

  /**
   * Structure of lattice properties needed by StressCauchyBorn.
   */
  struct CbData {
    double e2mvv;           //*> Energy conversion factor (1/mvv2e).
    double boltzmann;       //*> Boltzmann constant (in LAMMPS units)
    double hbar;            //*> Planck's constant (in LAMMPS units)
    double inv_atom_volume; //*> Volume of atom.
    double atom_mass;       //*> Mass of an atom.
    DENS_MAT cell_vectors;  //*> Unit vectors for lattice cells.
    DENS_MAT basis_vectors; //*> Positions of atoms within a lattice cell.
  };

  /**
   * @class StressCauchyBorn
   * @brief Class for computing the stress and elastic constants for a
   * @brief Cauchy-Born material.
   */


  class StressCauchyBorn : public Stress
  {
    public:
      StressCauchyBorn(std::fstream &matfile, CbData &cb);
      virtual ~StressCauchyBorn();
      virtual void initialize(void);
      //* Returns the stress computed from a 0K Cauchy-Born approxmation.
      virtual void stress(const FIELD_MATS &fields, const GRAD_FIELD_MATS &gradFields,
                  DENS_MAT_VEC &flux);
      //* Computes free (T>0)/potential(T=0) energy density
      //* Units: mvv/L^3 (i.e. for units Real: g/(mol ps^2 A^2) )
      virtual void elastic_energy(const FIELD_MATS &fields,
                          const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT &energy) const;
      //* Computes entropic energy density
      void entropic_energy(const FIELD_MATS &fields, const GRAD_FIELD_MATS &gradFields,
                          DENS_MAT &energy) const;
      //* Returns the material tangent at a given deformation gradient.
      virtual void tangent(const MATRIX &F, MATRIX &C) const;
      double stiffness() const;
      //* Creates a linearization for a deformation gradient.
      DENS_VEC elasticity_tensor(const VECTOR &Fv, MATRIX &C, const ElasticityTensorType type=FIRST_ELASTICITY_TENSOR) const;
      DENS_VEC elasticity_tensor(const MATRIX &F, MATRIX &C, const ElasticityTensorType type=FIRST_ELASTICITY_TENSOR) const;
    protected:
      void linearize(MATRIX *F=nullptr);
      CBLattice   *cblattice_;       //*> CbLattice -> makes atom clusters.
      CbPotential *potential_;       //*> CbPotential -> interatomic forces.
      bool makeLinear_;
      StressCubicElastic *cubicMat_; //*> Stores optional linear elastic law.
      bool initialized_;
      double fixed_temperature_;     //*> Specifies a uniform temperature.
      CbData cbdata_;                //*> Lattice & atom volume/mass.
  };


    // adaptor to NonLinearSolver
    class CBElasticTangentOperator : public TangentOperator {
      public:
        CBElasticTangentOperator (StressCauchyBorn * cauchyBornStress,
                                  DENS_VEC & targetP) :
          TangentOperator(),
          cauchyBornStress_(cauchyBornStress),
          targetP_(targetP) {};
        void function(const VECTOR & F, DENS_VEC & R)
        {
          DENS_MAT B;
          tangent(F,R,B);
        }
        void tangent(const VECTOR & F, DENS_VEC & R, MATRIX & B)
        {
          cbP_ = cauchyBornStress_->elasticity_tensor(F, B);

          R = cbP_ - targetP_;
        }
      private:
        StressCauchyBorn * cauchyBornStress_;
        DENS_VEC targetP_, cbP_;
    };

    // adaptor to NonLinearSolver
    class CB2ndElasticTangentOperator : public TangentOperator {
      public:
        CB2ndElasticTangentOperator (StressCauchyBorn * cauchyBornStress,
                                  DENS_VEC & targetS) :
          TangentOperator(),
          cauchyBornStress_(cauchyBornStress),
          targetS_(targetS) {};
        void function(const VECTOR & U, DENS_VEC & r)
        {
          DENS_MAT C;
          tangent(U,r,C);
        }
        void tangent(const VECTOR & U, DENS_VEC & r, MATRIX & C)
        {
          cbS_ = cauchyBornStress_->elasticity_tensor(U, C, SECOND_ELASTICITY_TENSOR);

          r = cbS_ - targetS_;
        }
      private:
        StressCauchyBorn * cauchyBornStress_;
        DENS_VEC targetS_, cbS_;
    };

}
#endif
