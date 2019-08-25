/** ATC_Transfer :  smoothing  */
#ifndef ATC_TRANSFER_KERNEL_H
#define ATC_TRANSFER_KERNEL_H

// ATC headers
#include "ATC_Transfer.h"

namespace ATC {

class KernelFunction;

class ATC_TransferKernel : public ATC_Transfer {

 public:
  
  // constructor
  ATC_TransferKernel(std::string groupName, 
                     double **& perAtomArray,
                     LAMMPS_NS::Fix * thisFix,
                     std::string matParamFile = "none");

  // destructor
  virtual ~ATC_TransferKernel();

  /** parser/modifier */
  virtual bool modify(int narg, char **arg);

 protected:
  /** routine to calculate matrix of kernel functions */
  
  virtual void compute_kernel_matrix_molecule();

  /** calculate projection on the fly*/
  virtual void compute_projection(const DENS_MAT & atomData,
                                  DENS_MAT & nodeData);

  /** routine to calculate matrix of bond functions */
  virtual void compute_bond_matrix();

  /** routine to calculate stress on-the-fly */
  virtual void compute_potential_stress(DENS_MAT& stress);

  /** routine to calculate force part of the heat flux on-the-fly */
  virtual void compute_potential_heatflux(DENS_MAT& flux);

  /** calculate dislocation density tensor from DXA output */
  virtual void compute_dislocation_density(DENS_MAT & dislocation_density);

};

};

#endif
