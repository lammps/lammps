#include "PerPairQuantity.h"
#include "PerAtomQuantity.h"
#include "KernelFunction.h"
#include "FE_Mesh.h"
#include "Utility.h"
#include "Quadrature.h"

using ATC::HeartBeat;
using std::pair;
using std::map;

namespace ATC {

//==========================================================
PairMap::PairMap(LammpsInterface * lammpsInterface, int groupbit ):
  lammpsInterface_(lammpsInterface),
  groupbit_(groupbit),
  nPairs_(0), nBonds_(0)
{
};
PairMap::~PairMap(void)
{
};
//==========================================================
PairMapNeighbor::PairMapNeighbor(LammpsInterface * lammpsInterface, int groupbit):
  PairMap(lammpsInterface,groupbit)
{
};

void PairMapNeighbor::reset(void) const
{
  int inum = lammpsInterface_->neighbor_list_inum();
  int *ilist = lammpsInterface_->neighbor_list_ilist();
  int *numneigh = lammpsInterface_->neighbor_list_numneigh();
  int **firstneigh = lammpsInterface_->neighbor_list_firstneigh();
  const int * mask = lammpsInterface_->atom_mask();

  pairMap_.clear();
  int pairIndex = nBonds_;
  std::pair< int,int > pair_ij;
  for (int i = 0; i < inum; i++) {
    int lammps_i = ilist[i];
    if (mask[lammps_i] & groupbit_) {
      for (int j = 0; j < numneigh[lammps_i]; j++) {
        int lammps_j = firstneigh[lammps_i][j];
        lammpsInterface_->neighbor_remap(lammps_j);
        pair_ij.first  = lammps_i; // alpha
        pair_ij.second = lammps_j; // beta
        pairMap_[pair_ij] = pairIndex;
        pairIndex++;
      }
    }
  }
  nPairs_ = pairIndex;
  needReset_ = false;
}

//==========================================================
PairMapBond::PairMapBond(LammpsInterface * lammpsInterface, int groupbit):
  PairMap(lammpsInterface,groupbit)
{
};

//==========================================================
PairMapBoth::PairMapBoth(LammpsInterface * lammpsInterface, int groupbit):
  PairMapNeighbor(lammpsInterface,groupbit)
{
};
//==========================================================
DensePerPairMatrix::DensePerPairMatrix(LammpsInterface * lammpsInterface,
  const PairMap & pairMap,
  int nCols):
  lammpsInterface_(lammpsInterface),
  pairMap_(pairMap),
  nCols_(nCols)
{
};

//==========================================================
PairVirial::PairVirial(LammpsInterface * lammpsInterface,
   const PairMap & pairMap, int nCols):
  DensePerPairMatrix(lammpsInterface,pairMap,nCols)
{
};
//==========================================================
PairVirialEulerian::PairVirialEulerian(LammpsInterface * lammpsInterface,
   const PairMap & pairMap):
  PairVirial(lammpsInterface,pairMap,6)
{
};


void PairVirialEulerian::reset(void) const
{
  int nPairs = pairMap_.size();
  quantity_.reset(nPairs,nCols_);

  double **xatom = lammpsInterface_->xatom();
  for (ATOM_PAIR apair = pairMap_.start();
       ! pairMap_.finished(); apair=pairMap_++){
    int lammps_a = (apair.first).first ;
    int lammps_b = (apair.first).second;
    int pairIndex = apair.second;
    double * xa = xatom[lammps_a];
    double * xb = xatom[lammps_b];
    double delx = xa[0] - xb[0];
    double dely = xa[1] - xb[1];
    double delz = xa[2] - xb[2];
    double rsq = delx*delx + dely*dely + delz*delz;
    double fforce = 0;
    lammpsInterface_->pair_force(apair,rsq,fforce);
    quantity_(pairIndex,0)=-delx*delx*fforce;
    quantity_(pairIndex,1)=-dely*dely*fforce;
    quantity_(pairIndex,2)=-delz*delz*fforce;
    quantity_(pairIndex,3)=-delx*dely*fforce;
    quantity_(pairIndex,4)=-delx*delz*fforce;
    quantity_(pairIndex,5)=-dely*delz*fforce;
  }
}
//==========================================================
PairVirialLagrangian::PairVirialLagrangian(LammpsInterface * lammpsInterface,
   const PairMap & pairMap,
   double ** xRef):
// const PerAtomQuantity<double> * xRef):
  PairVirial(lammpsInterface,pairMap,9),
  xRef_(xRef)
{

};


void PairVirialLagrangian::reset(void) const
{
  int nPairs = pairMap_.size();
  quantity_.reset(nPairs,nCols_);

  double **xatom = lammpsInterface_->xatom();

  double ** xref = xRef_;

  for (ATOM_PAIR apair = pairMap_.start();
       ! pairMap_.finished(); apair=pairMap_++){
    int lammps_a = (apair.first).first ;
    int lammps_b = (apair.first).second;
    int pairIndex = apair.second;
    double * xa = xatom[lammps_a];
    double * xb = xatom[lammps_b];
    double delx = xa[0] - xb[0];
    double dely = xa[1] - xb[1];
    double delz = xa[2] - xb[2];
    double * Xa = xref[lammps_a];
    double * Xb = xref[lammps_b];
    double delX = Xa[0] - Xb[0];
    double delY = Xa[1] - Xb[1];
    double delZ = Xa[2] - Xb[2];
    double rsq = delx*delx + dely*dely + delz*delz;
    double fforce = 0;
    lammpsInterface_->pair_force(apair,rsq,fforce);
    quantity_(pairIndex,0)=-delx*fforce*delX;
    quantity_(pairIndex,1)=-delx*fforce*delY;
    quantity_(pairIndex,2)=-delx*fforce*delZ;
    quantity_(pairIndex,3)=-dely*fforce*delX;
    quantity_(pairIndex,4)=-dely*fforce*delY;
    quantity_(pairIndex,5)=-dely*fforce*delZ;
    quantity_(pairIndex,6)=-delz*fforce*delX;
    quantity_(pairIndex,7)=-delz*fforce*delY;
    quantity_(pairIndex,8)=-delz*fforce*delZ;
  }
}
//==========================================================
PairPotentialHeatFlux::PairPotentialHeatFlux(LammpsInterface * lammpsInterface,
   const PairMap & pairMap):
  DensePerPairMatrix(lammpsInterface,pairMap,3)
{
};
//==========================================================
PairPotentialHeatFluxEulerian::PairPotentialHeatFluxEulerian(LammpsInterface * lammpsInterface,
   const PairMap & pairMap):
  PairPotentialHeatFlux(lammpsInterface,pairMap)
{

};

void PairPotentialHeatFluxEulerian::reset(void) const
{
  int nPairs = pairMap_.size();
  quantity_.reset(nPairs,nCols_);

  double **xatom = lammpsInterface_->xatom();
  double **vatom = lammpsInterface_->vatom();
  for (ATOM_PAIR apair = pairMap_.start();
       ! pairMap_.finished(); apair=pairMap_++){
    int lammps_a = (apair.first).first ;
    int lammps_b = (apair.first).second;
    int pairIndex = apair.second;
    double * xa = xatom[lammps_a];
    double * xb = xatom[lammps_b];
    double delx = xa[0] - xb[0];
    double dely = xa[1] - xb[1];
    double delz = xa[2] - xb[2];
    double rsq = delx*delx + dely*dely + delz*delz;
    double fforce = 0;
    lammpsInterface_->pair_force(apair,rsq,fforce);
    double* v = vatom[lammps_a];
    fforce *=delx*v[0] + dely*v[1] + delz*v[2];
    quantity_(pairIndex,0)=fforce*delx;
    quantity_(pairIndex,1)=fforce*dely;
    quantity_(pairIndex,2)=fforce*delz;
  }
}
//==========================================================
PairPotentialHeatFluxLagrangian::PairPotentialHeatFluxLagrangian(LammpsInterface * lammpsInterface,
   const PairMap & pairMap, double ** xRef):
  PairPotentialHeatFlux(lammpsInterface,pairMap),
  xRef_(xRef)
{

};

void PairPotentialHeatFluxLagrangian::reset(void) const
{
  int nPairs = pairMap_.size();
  quantity_.reset(nPairs,nCols_);

  double **xatom = lammpsInterface_->xatom();
  double **vatom = lammpsInterface_->vatom();

  for (ATOM_PAIR apair = pairMap_.start();
       ! pairMap_.finished(); apair=pairMap_++){
    int lammps_a = (apair.first).first ;
    int lammps_b = (apair.first).second;
    int pairIndex = apair.second;
    double * xa = xatom[lammps_a];
    double * xb = xatom[lammps_b];
    double delx = xa[0] - xb[0];
    double dely = xa[1] - xb[1];
    double delz = xa[2] - xb[2];
    double * Xa = xRef_[lammps_a];
    double * Xb = xRef_[lammps_b];
    double delX = Xa[0] - Xb[0];
    double delY = Xa[1] - Xb[1];
    double delZ = Xa[2] - Xb[2];
    double rsq = delx*delx + dely*dely + delz*delz;
    double fforce = 0;
    lammpsInterface_->pair_force(apair,rsq,fforce);
    double* v = vatom[lammps_a];
    fforce *=delx*v[0] + dely*v[1] + delz*v[2];
    quantity_(pairIndex,0)=fforce*delX;
    quantity_(pairIndex,1)=fforce*delY;
    quantity_(pairIndex,2)=fforce*delZ;
  }
}
//==========================================================
SparsePerPairMatrix::SparsePerPairMatrix(LammpsInterface * lammpsInterface,
  const PairMap & pairMap):
  lammpsInterface_(lammpsInterface),
  pairMap_(pairMap)
{
};
//==========================================================
BondMatrix::BondMatrix(LammpsInterface * lammpsInterface,
  const PairMap & pairMap, double ** x, const FE_Mesh * feMesh):
  SparsePerPairMatrix(lammpsInterface,pairMap), x_(x), feMesh_(feMesh)
{
};
//==========================================================
BondMatrixKernel::BondMatrixKernel(LammpsInterface * lammpsInterface,
  const PairMap & pairMap,
  double ** x,
  const FE_Mesh * feMesh,
  const KernelFunction * kernelFunction):
  BondMatrix(lammpsInterface,pairMap,x,feMesh),
  kernelFunction_(kernelFunction)
{
  if (kernelFunction_ == nullptr)
    throw ATC_Error("No AtC kernel function initialized");
};
void BondMatrixKernel::reset(void) const
{
  int nPairs = pairMap_.size(); // needs to come after quantity for reset
  int nNodes = feMesh_->num_nodes_unique();
  quantity_.reset(nNodes,nPairs);
  double lam1,lam2;
  int heartbeatFreq = (nNodes <= 10 ? 1 : (int) nNodes / 10);
  HeartBeat beat("computing bond matrix ",heartbeatFreq);
  beat.start();
  DENS_VEC xa(3),xI(3),xaI(3),xb(3),xbI(3),xba(3);
  double invVol = kernelFunction_->inv_vol();
  for (int I = 0; I < nNodes; I++) {
    beat.next();
    xI = feMesh_->nodal_coordinates(I);
    if (!kernelFunction_->node_contributes(xI)) { continue; }
    for (ATOM_PAIR apair = pairMap_.start();
       ! pairMap_.finished(); apair=pairMap_++){
      int lammps_a = (apair.first).first ;
      int lammps_b = (apair.first).second;
      xa.copy(x_[lammps_a],3);
      xaI = xa - xI;
      lammpsInterface_->periodicity_correction(xaI.ptr());
      xb.copy(x_[lammps_b],3);
      xba = xb - xa;
      xbI = xba + xaI;
      kernelFunction_->bond_intercepts(xaI,xbI,lam1,lam2);
      if (lam1 < lam2) {
        double bondValue = invVol*(kernelFunction_->bond(xaI,xbI,lam1,lam2));
        int pairIndex = apair.second;
        quantity_.add(I,pairIndex,bondValue);
      } // if lam1 < lam2
    } // pair map
  } // end nodes loop
  quantity_.compress();
  beat.finish();
}
//==========================================================
BondMatrixPartitionOfUnity::BondMatrixPartitionOfUnity(LammpsInterface * lammpsInterface,
  const PairMap & pairMap, double ** x,  const FE_Mesh * feMesh,
  const DIAG_MAN * invVols):
  BondMatrix(lammpsInterface,pairMap,x,feMesh),
  invVols_(invVols)
{
  ATC::Quadrature::instance()->set_line_quadrature(lineNgauss_,lineXg_,lineWg_);
  double lam1 = 0.0, lam2 = 1.0;
  double del_lambda = 0.5*(lam2 - lam1);
  double avg_lambda = 0.5*(lam2 + lam1);
  for (int i = 0; i < lineNgauss_; i++) {
    double lambda = del_lambda*lineXg_[i] +avg_lambda;
    lineXg_[i] = lambda;
    lineWg_[i] *= 0.5;
  }
};
void BondMatrixPartitionOfUnity::reset(void) const
{
  int nNodes = feMesh_->num_nodes_unique();
  int nPairs = pairMap_.size();
  quantity_.reset(nNodes,nPairs);
  int nodes_per_element = feMesh_->num_nodes_per_element();
  Array<int> node_list(nodes_per_element);
  DENS_VEC shp(nodes_per_element);
  int heartbeatFreq = (int) nPairs / 10;
  HeartBeat beat("computing bond matrix ",heartbeatFreq);
  beat.start();
  DENS_VEC xa(3),xb(3),xab(3),xlambda(3);
  for (ATOM_PAIR apair = pairMap_.start();
     ! pairMap_.finished(); apair=pairMap_++){
    beat.next();
    int lammps_a = (apair.first).first ;
    int lammps_b = (apair.first).second;
    int pairIndex = apair.second;
    xa.copy(x_[lammps_a],3);
    xb.copy(x_[lammps_b],3);
    xab = xa - xb;
    for (int i = 0; i < lineNgauss_; i++) {
      double lambda = lineXg_[i];
      xlambda = lambda*xab + xb;
      lammpsInterface_->periodicity_correction(xlambda.ptr());
      feMesh_->shape_functions(xlambda,shp,node_list);
      // accumulate to nodes whose support overlaps the integration point
      for (int I = 0; I < nodes_per_element; I++) {
        int Inode = node_list(I);
        double inv_vol = (invVols_->quantity())(Inode,Inode);
        double val = inv_vol*shp(I)*lineWg_[i];
        quantity_.add(Inode,pairIndex,val);
      }
    }
  }
  quantity_.compress();
  beat.finish();
}
//==========================================================
}
