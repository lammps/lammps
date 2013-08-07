#include "CBLattice.h"
#include "CbPotential.h"

namespace ATC {
  // removes any overlapping atoms (avoid calling because it scales n^2.)
  INDEX AtomCluster::remove_overlap()
  {
    INDEX removed_count = 0;
    std::vector<double>::iterator r(cur_bond_len_.begin());
    std::vector<DENS_VEC>::iterator R(ref_coords_.begin()), Rp;
    const double TOL = 1.0e-6 * R->dot(*R);

    for (; R!=ref_coords_.end(); R++, r++) { 
      for (Rp=R+1; Rp!=ref_coords_.end(); Rp++) {
        if (sum_difference_squared(*Rp, *R) < TOL) {
          ref_coords_.erase(R--);
          cur_bond_len_.erase(r--);
          ++removed_count;
          break;
        }
      }
    }
    return removed_count;
  }
  //=========================================================================
  // writes cluster to 3 column data format
  //=========================================================================
  void AtomCluster::write_to_dat(std::string path, bool current_config)
  {
    const int npts = int(ref_coords_.size());
    if (path.substr(path.size()-5,4) != ".dat") path += ".dat";
    fstream fid(path.c_str(), std::ios::out);

    for (int i=0; i<npts; i++) {
      DENS_VEC x (current_config ? r(i):R(i));
      for (INDEX j=0; j<x.size(); j++) fid << x(j) << " ";
      fid << " " << x.norm() << "\n";
    }
  }
  //=========================================================================
  // writes cluster to vtk format, (in either reference or current config)
  //=========================================================================
  void AtomCluster::write_to_vtk(std::string path, bool current_config)
  {
    const int npts = int(ref_coords_.size());
    if (path.substr(path.size()-5,4) != ".vtk") path += ".vtk";
    fstream fid(path.c_str(), std::ios::out);

    // write header
    fid << "# vtk DataFile Version 2.0\nWritten from FE-LAMMPS\n";
    fid << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
    fid << "POINTS " << npts << " float\n";

    for (int i=0; i<npts; i++) {
      DENS_VEC x (current_config ? r(i):R(i));
      for (INDEX j=0; j<x.size(); j++) fid << x(j) << " ";
      fid << ((i+1)%3==0 ? "\n" : " ");
    }
    fid << "\nCELLS "<<npts<<" "<<2*npts<<"\n";
    for (int i=0; i<npts; i++) fid << "1" << " " << i << "\n"; 

    fid << "CELL_TYPES " << npts << "\n";
    for (int i=0; i<npts; i++) fid << "1" << "\n";
  }
  //===========================================================================
  // constructor 
  // @param N 3x3 DenseMatrix with each column as a base vector
  // @param B 3xn DenseMatrix with each column being a basis vector
  // @param R vector of 3D bond Vectors to representative atom in ref config
  // @param r vector of bond lengths to representative atom in current config
  // @param RC cutoff radius of bond potential
  //===========================================================================
  CBLattice::CBLattice(const MATRIX &N, const MATRIX &B)
    : n_(N), b_(B), N_(N), B_(B)
  {
    //  builds the default queue
    for (int a=-1; a<2; a++)
      for (int b=-1; b<2; b++)
        for (int c=-1; c<2; c++)
           if ( a!=0 || b!=0 || c!=0) queue0_.push(hash(a,b,c));  
  }
  //=============================================================================
  // writes out default lattice parameters
  //=============================================================================
  std::ostream & operator<<(std::ostream& o, const CBLattice& lattice)
  {
    o<<"cutoff radius = "<<sqrt(lattice.RC2_)<<"\n";
    lattice.N_.print("Reference base vectors");
    lattice.B_.print("Reference basis vectors");
    return o;
  }
  //===========================================================================
  // Constructs the virtual atom cluster of neighbors within cutoff
  // @param F the deformation gradient tensor
  //===========================================================================
  void CBLattice::atom_cluster(const MATRIX &F, double cutoff, AtomCluster &v)
  {
    RC2_ = cutoff*cutoff;
    // compute new base and basis vectors
    v.clear();
    v.F_ = F;
    n_ = F*N_;
    b_ = F*B_;
    // add basis from the center cell (not including representative atom)
    for (int i=1; i<B_.nCols(); i++) {
     v.cur_bond_len_.push_back(b_.col_norm(i));
     v.ref_coords_.push_back(column(B_,i));
    }
    _FindAtomsInCutoff(v);
  }
  //=============================================================================
  // Computes forces on central atom, with atom I displaced by u. 
  //=============================================================================
  DENS_VEC AtomCluster::perturbed_force(const CbPotential *p, int I, DENS_VEC *u) const
  {
    DENS_VEC f(3);
    for (INDEX i=0; i<size(); i++) {
      DENS_VEC ri = r(i);
      if (u && i+1==I) ri += *u;
      if (u && I==0)   ri -= *u;
      const double d = ri.norm();
      f.add_scaled(ri, -p->phi_r(d)/d);
    }
    return f;
  }
  //=============================================================================
  // Computes the force constant matrix between atoms I and 0.
  //=============================================================================
  DENS_MAT AtomCluster::force_constants(INDEX I, const CbPotential *p) const
  {
    DENS_MAT D(3,3);
    for (INDEX i=0; i<3; i++) {
      DENS_VEC du(3);
      du(i) = 1.0e-6;  // take central difference
      row(D,i)  = perturbed_force(p, I, &du);
      du(i)     = -du(i);
      row(D,i) -= perturbed_force(p, I, &du);
    }
    D *= 0.5e6;
    return D;
  }
  //=============================================================================
  // performs an iterative search for all neighbors within cutoff
  //=============================================================================
  void CBLattice::_FindAtomsInCutoff(AtomCluster &v)
  {
    static const int dir[2] = {-1, 1};
    int a, b, c, abc;
    std::stack<int> queue(queue0_);
    std::set<int> done;
    // search in each direction
    while (!queue.empty()) {
       abc = queue.top();
       queue.pop();
       if (done.find(abc) == done.end()) { // value not in set
         unhash(abc,a,b,c);             // convert abc to a,b,c
         if (_CheckUnitCell(a, b, c, v)) {
           done.insert(abc);
           // add direct 'outward' neighbors to queue
           queue.push(hash(a+dir[a>0], b, c));
           queue.push(hash(a, b+dir[b>0], c));
           queue.push(hash(a, b, c+dir[c>0]));
         }
       }
    }
  }
  //===========================================================================
  // Computes \f$r^2 = \Vert a n_1 + b n_2 +c n_3 + b_d \Vert^2 \f$ 
  // and adds atom (a,b,c,d) if \f$r^2 < r_{cutoff}^2 \f$
  // @param a cell x-index
  // @param b cell y-index
  // @param c cell z-index
  //===========================================================================
  bool CBLattice::_CheckUnitCell(char a, char b, char c, AtomCluster &v) 
  {
    const int nsd = n_.nRows();  // number of spatial dimensions
    const double A=double(a), B=double(b), C=double(c);
    bool found=false;
    DENS_VEC r0(nsd,false), R0(nsd,false), Rd(nsd,false);  // don't initialize
  
    for (int i=0; i<nsd; i++) { // precompute locations of cell
      R0(i) = A*N_(0,i) + B*N_(1,i) + C*N_(2,i); // reference
    }
    r0 = v.F_*R0; // deformed
    for (int d=0; d<b_.nCols(); d++) {
      double ri = r0(0) + b_(0,d);
      double r2 = ri*ri;
      for (int i=1; i<nsd; i++) {
        ri = r0(i) + b_(i,d);
        r2 += ri*ri;
      }

      if (r2 <= RC2_) {
        v.ref_coords_.push_back(R0);
        v.ref_coords_.back() += column(B_,d);  // position is R0 + B_[d]
        v.cur_bond_len_.push_back(sqrt(r2));
        found = true;
      }
    }
    return found;
  }
}  // end ATC
