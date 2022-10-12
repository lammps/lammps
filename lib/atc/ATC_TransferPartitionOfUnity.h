#ifndef ATC_TRANSFER_PARTITION_OF_UNITY_H
#define ATC_TRANSFER_PARTITION_OF_UNITY_H

// ATC headers
#include "ATC_Transfer.h"

namespace ATC {

class ATC_TransferPartitionOfUnity : public ATC_Transfer {

 public:

  // constructor
  ATC_TransferPartitionOfUnity(std::string groupName,
                               double **& perAtomArray,
                               LAMMPS_NS::Fix * thisFix,
                               std::string matParamFile = "none");

  // destructor
  virtual ~ATC_TransferPartitionOfUnity();

 protected:
  /** routine to calculate matrix of kernel functions */
  virtual void compute_kernel_matrix_molecule() {};

  /** calculate projection on the fly*/
  virtual void compute_projection(const DENS_MAT & atomData,
                                  DENS_MAT & nodeData);

  /** routine to calculate matrix of bond functions */
  virtual void compute_bond_matrix();

  /** routine to calculate kinetic part of stress */
  virtual void compute_kinetic_stress(DENS_MAT& stress);

  /** routine to calculate stress on-the-fly */
  virtual void compute_potential_stress(DENS_MAT& stress);

  /** routine to calculate kinetic part of heat flux */
  virtual void compute_kinetic_heatflux(DENS_MAT& flux);

  /** routine to calculate force part of the heat flux on-the-fly */
  virtual void compute_potential_heatflux(DENS_MAT& flux);

  /** compute velocity & difference betw atomic and (local) mean velocity */
  void compute_variation_velocity();

  /** calculate dislocation density tensor from DXA output */
  virtual void compute_dislocation_density(DENS_MAT & dislocation_density);

 private:

  DENS_MAT variationVelocity_;
  DENS_MAT vbar_;
};

}

#endif
