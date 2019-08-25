#ifndef CBLATTICE_H
#define CBLATTICE_H

#include <set>
#include <vector>
#include <stack>
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"

namespace ATC
{
  class CbPotential; // forward definition for potentials.
  class StressArgs;

  /**
   *  @class  AtomCluster
   *  @brief  Base class for a container for a cluster of atoms around an atom at the origin 
   *          (which is not included in the lists).
   *          Provides routines for outputting the atoms in the cluster to a vtk file,
   *          and checking for any overlapping atoms (which should not happen).
   */

  class AtomCluster
  {
    friend class CBLattice;
    friend DENS_MAT_VEC compute_dynamical_derivative(StressArgs &args); 
    friend int test_FCB(const StressArgs &args);
  public:
    //* Returns the number of atoms in the cluster.
    INDEX size() const { return cur_bond_len_.size(); }
    //* removes all atoms from the cluster
    void clear() { cur_bond_len_.clear(), ref_coords_.clear(); }
    //* Removes any overlapping atoms (if they exist).
    INDEX remove_overlap();
    //* Writes the atom positions of the cluster to a dat file.
    void write_to_dat(std::string path, bool current_config=true);
    //* Writes the atom positions of the cluster to a vtk file.
    void write_to_vtk(std::string path, bool current_config=true);

    //* Returns an atom coordinate in the reference configuration.
    const DENS_VEC &R(INDEX i) const { return ref_coords_[i]; }
    //* Returns an atom coordinate in the current configuration.
    DENS_VEC r(INDEX i) const { return F_*R(i); }

    //* Returns the ith atoms bond length to the central atom.
    double bond_length(INDEX i) const { return cur_bond_len_[i]; }
    //* Returns a reference to the deformation gradient tensor.
    const DENS_MAT& deformation_gradient() const { return F_; }
  
    //* Computes forces on central atom, with atom I displaced by u. 
    DENS_VEC perturbed_force(const CbPotential *p, int I, DENS_VEC *u) const;
    //* Computes the force constant matrix between atoms I and 0.
    DENS_MAT force_constants(INDEX I, const CbPotential *p) const;

  private:
    
    std::vector<double>   cur_bond_len_; //*> Bond lengths (current)
    std::vector<DENS_VEC> ref_coords_;   //*> Atom coordinates (ref)
    DENS_MAT F_;                         //*> Deformation gradient
  };

  /**
   *  @class  CBLattice 
   *  @brief  Base class that generates a virtual atom clusters given a lattice and 
   *          a deformation gradient.
   */
 
  class CBLattice
  {
  protected:
    double RC2_;             //*> cutoff radius^2
    std::stack<int> queue0_; //*> cell nghbrs of representative cell
    DENS_MAT n_, b_;         //*> cols are def base and basis vectors
    const MATRIX &N_, &B_;   //*> cols are ref base and basis vectors

  public:
    //* Operator that outputs the lattice and basis to a stream.
    friend std::ostream& operator<<(std::ostream& o, const CBLattice& lattice);
    CBLattice(const MATRIX &N, const MATRIX &B); 
    //* generates the virtual atom cluster
    void atom_cluster(const MATRIX &F, double cutoff, AtomCluster &v);


  protected:
    void _FindAtomsInCutoff(AtomCluster &v);
    bool _CheckUnitCell(char a, char b, char c, AtomCluster &v);
  };
  // hash functions: a, b, c must in range [-128, 127]
  inline int hash(int a,int b,int c) {return(a+128)|((b+128)<<8)|((c+128)<<16);}
  inline void unhash(int r, int &a, int &b, int &c)   
  {
     a = (r&0xFF) - 128;
     b = ((r>>8)&0xFF) - 128;
     c = ((r>>16)&0xFF) - 128;
  }
}
#endif
