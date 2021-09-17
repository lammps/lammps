// ATC headers
#include "ATC_TransferKernel.h"
#include "ATC_Error.h"
#include "FE_Engine.h"
#include "KernelFunction.h"
#include "LammpsInterface.h"
#include "Quadrature.h"
#include "PerPairQuantity.h"
#include "Stress.h"
#ifdef HAS_DXA
#include "DislocationExtractor.h"
#endif

// Other Headers
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <fstream>
#include <exception>

using namespace std;

namespace ATC {

using ATC_Utility::to_string;

  ATC_TransferKernel::ATC_TransferKernel(string groupName,
                                         double **& perAtomArray,
                                         LAMMPS_NS::Fix * thisFix,
                                         string matParamFile)
    : ATC_Transfer(groupName,perAtomArray,thisFix,matParamFile)
  {
    kernelBased_ = true;
  }

  //-------------------------------------------------------------------
  ATC_TransferKernel::~ATC_TransferKernel()
  {
  }

  //-------------------------------------------------------------------
  bool ATC_TransferKernel::modify(int narg, char **arg)
  {
    bool match = false;

    /*! \page man_hardy_kernel fix_modify AtC kernel
      \section syntax
      fix_modify AtC kernel <type> <parameters>
      - type (keyword) = step, cell, cubic_bar, cubic_cylinder, cubic_sphere,
                         quartic_bar, quartic_cylinder, quartic_sphere \n
      - parameters :\n
      step = radius (double) \n
      cell = hx, hy, hz (double) or h (double) \n
      cubic_bar = half-width (double) \n
      cubic_cylinder = radius (double) \n
      cubic_sphere = radius (double) \n
      quartic_bar = half-width (double) \n
      quartic_cylinder = radius (double) \n
      quartic_sphere = radius (double) \n
      \section examples
      <TT> fix_modify AtC kernel cell 1.0 1.0 1.0 </TT> \n
      <TT> fix_modify AtC kernel quartic_sphere 10.0 </TT>
      \section description

      \section restrictions
      Must be used with the hardy AtC fix \n
      For bar kernel types, half-width oriented along x-direction \n
      For cylinder kernel types, cylindrical axis is assumed to be in z-direction \n
      ( see \ref man_fix_atc )
      \section related
      \section default
      No default
    */

    // no match, call base class parser
    if (!match) {
      match = ATC_Transfer::modify(narg, arg);
    }

    return match;

  }

  //-------------------------------------------------------------------
  void ATC_TransferKernel::compute_kernel_matrix_molecule(void) // KKM add
  {
    int nLocalMol = smallMoleculeSet_->local_molecule_count();
    if (nLocal_>0) {
      SPAR_MAT & N(kernelAccumulantMol_.set_quantity());
      N.reset(nLocalMol,nNodes_);
      SPAR_MAT & dN(kernelAccumulantMolGrad_.set_quantity());
      dN.reset(nLocalMol,nNodes_);
      DENS_VEC derivKer(nsd_);
      DENS_VEC xI(nsd_),xm(nsd_),xmI(nsd_);
      const DENS_MAT & centroidMolMatrix(moleculeCentroid_->quantity());
      ATC::LammpsInterface::instance()->stream_msg_once("computing kernel matrix molecule ",true,false);
      int heartbeatFreq = (nNodes_ <= 10 ? 1 : (int) nNodes_ / 10);
      for (int i = 0; i < nNodes_; i++) {
        if (i % heartbeatFreq == 0 ) {
          ATC::LammpsInterface::instance()->stream_msg_once(".",false,false);
        }
        xI = (feEngine_->fe_mesh())->nodal_coordinates(i);
        for (int j = 0; j < nLocalMol; j++) {
          for (int k = 0; k < nsd_; k++) {
            xm(k) = centroidMolMatrix(j,k);
          }
          xmI = xm - xI;
          lammpsInterface_->periodicity_correction(xmI.ptr());
          double val = kernelFunction_->value(xmI);
          if (val > 0) N.add(j,i,val);
          kernelFunction_->derivative(xmI,derivKer);
          double val_grad = derivKer(2);
          if (val_grad!= 0) dN.add(j,i,val_grad);
        }
      }
      // reset kernelShpFunctions with the weights of molecules on processors
      DENS_VEC fractions(N.nRows());
      DENS_VEC fractions_deriv(dN.nRows());
      for (int i = 0; i < nLocalMol; i++) {
        fractions(i) = smallMoleculeSet_->local_fraction(i);
      }
        N.row_scale(fractions);
        N.compress();
        dN.row_scale(fractions);
        dN.compress();
      if (lammpsInterface_->rank_zero()) {
        ATC::LammpsInterface::instance()->stream_msg_once("done",false,true);
      }
    }
  }

  //-------------------------------------------------------------------
  void ATC_TransferKernel::compute_projection(const DENS_MAT & atomData,
                                              DENS_MAT & nodeData)
  {
    DENS_MAT workNodeArray(nNodes_, atomData.nCols());
    workNodeArray.zero();
    nodeData.reset(workNodeArray.nRows(),workNodeArray.nCols());
    nodeData.zero();

    if (nLocal_>0) {
      set_xPointer();
      DENS_VEC xI(nsd_),xa(nsd_),xaI(nsd_);
      double val;
      for (int i = 0; i < nNodes_; i++) {
        xI = (feEngine_->fe_mesh())->nodal_coordinates(i);
        for (int j = 0; j < nLocal_; j++) {
          int lammps_j = internalToAtom_(j);
          xa.copy(xPointer_[lammps_j],3);
          xaI = xa - xI;
          lammpsInterface_->periodicity_correction(xaI.ptr());
          val = kernelFunction_->value(xaI);
          if (val > 0) {
            for (int k=0; k < atomData.nCols(); k++) {
              workNodeArray(i,k) += val*atomData(j,k);
            }
          }
        }
      }
    }

    // accumulate across processors
    int count = workNodeArray.nRows()*workNodeArray.nCols();
    lammpsInterface_->allsum(workNodeArray.ptr(),nodeData.ptr(),count);
  }


  //-------------------------------------------------------------------
  void ATC_TransferKernel::compute_bond_matrix()
  {
    atomicBondMatrix_=bondMatrix_->quantity();
  }

  //-------------------------------------------------------------------
  // on-the-fly calculation of stress
  void ATC_TransferKernel::compute_potential_stress(DENS_MAT& stress)
  {
    set_xPointer();
    stress.zero();
    // neighbor lists
    int *numneigh = lammpsInterface_->neighbor_list_numneigh();
    int **firstneigh = lammpsInterface_->neighbor_list_firstneigh();
    double ** xatom    = lammpsInterface_->xatom();
    double lam1,lam2;
    double bond_value;
    // process differently for mesh vs translation-invariant kernels
    ATC::LammpsInterface::instance()->stream_msg_once("computing potential stress: ",true,false);
    int heartbeatFreq = (nNodes_ <= 10 ? 1 : (int) nNodes_ / 10);
    // "normal" kernel functions
    DENS_VEC xa(nsd_),xI(nsd_),xaI(nsd_),xb(nsd_),xbI(nsd_),xba(nsd_);
    double kernel_inv_vol = kernelFunction_->inv_vol();
    for (int i = 0; i < nNodes_; i++) {
      if (i % heartbeatFreq == 0 ) {
        ATC::LammpsInterface::instance()->stream_msg_once(".",false,false);
        }
      //  point
      xI = (feEngine_->fe_mesh())->nodal_coordinates(i);
      if (!kernelFunction_->node_contributes(xI)) {
        continue;
      }
      int inode = i;
      for (int j = 0; j < nLocal_; j++) {
        // second (neighbor) atom location
        int lammps_j = internalToAtom_(j);
        xa.copy(xPointer_[lammps_j],3);
        // difference vector
        xaI = xa - xI;
        lammpsInterface_->periodicity_correction(xaI.ptr());
        for (int k = 0; k < numneigh[lammps_j]; ++k) {
          int lammps_k = firstneigh[lammps_j][k];
          // first atom location
          xb.copy(xPointer_[lammps_k],3);
          // difference vector
          xba = xb - xa;
          xbI = xba + xaI;
          kernelFunction_->bond_intercepts(xaI,xbI,lam1,lam2);
          // compute virial
          if (lam1 < lam2) {
            bond_value
              = kernel_inv_vol*(kernelFunction_->bond(xaI,xbI,lam1,lam2));
            double delx = xatom[lammps_j][0] - xatom[lammps_k][0];
            double dely = xatom[lammps_j][1] - xatom[lammps_k][1];
            double delz = xatom[lammps_j][2] - xatom[lammps_k][2];
            double rsq = delx*delx + dely*dely + delz*delz;
            double fforce = 0;
            lammpsInterface_->pair_force(lammps_j,lammps_k,rsq,fforce);
            fforce *= 0.5; // dbl count
            if (atomToElementMapType_ == LAGRANGIAN) {
              double delX = xref_[lammps_j][0] - xref_[lammps_k][0];
              double delY = xref_[lammps_j][1] - xref_[lammps_k][1];
              double delZ = xref_[lammps_j][2] - xref_[lammps_k][2];
              stress(inode,0) +=-delx*fforce*delX*bond_value;
              stress(inode,1) +=-delx*fforce*delY*bond_value;
              stress(inode,2) +=-delx*fforce*delZ*bond_value;
              stress(inode,3) +=-dely*fforce*delX*bond_value;
              stress(inode,4) +=-dely*fforce*delY*bond_value;
              stress(inode,5) +=-dely*fforce*delZ*bond_value;
              stress(inode,6) +=-delz*fforce*delX*bond_value;
              stress(inode,7) +=-delz*fforce*delY*bond_value;
              stress(inode,8) +=-delz*fforce*delZ*bond_value;
            }
            else { //EULERIAN
              stress(inode,0) +=-delx*delx*fforce*bond_value;
              stress(inode,1) +=-dely*dely*fforce*bond_value;
              stress(inode,2) +=-delz*delz*fforce*bond_value;
              stress(inode,3) +=-delx*dely*fforce*bond_value;
              stress(inode,4) +=-delx*delz*fforce*bond_value;
              stress(inode,5) +=-dely*delz*fforce*bond_value;
            }
          }
        }
      }
    }
    ATC::LammpsInterface::instance()->stream_msg_once("done",false,true);
  }

  //-------------------------------------------------------------------
  // on-the-fly calculation of the heat flux
  void ATC_TransferKernel::compute_potential_heatflux(DENS_MAT& flux)
  {
    set_xPointer();
    flux.zero();
    // neighbor lists
    int *numneigh = lammpsInterface_->neighbor_list_numneigh();
    int **firstneigh = lammpsInterface_->neighbor_list_firstneigh();
    double ** xatom    = lammpsInterface_->xatom();
    double ** vatom    = lammpsInterface_->vatom();

    double lam1,lam2;
    double bond_value;
    // process differently for mesh vs translation-invariant kernels
    // "normal" kernel functions
    DENS_VEC xa(nsd_),xI(nsd_),xaI(nsd_),xb(nsd_),xbI(nsd_),xba(nsd_);
    double kernel_inv_vol = kernelFunction_->inv_vol();
    for (int i = 0; i < nNodes_; i++) {
      int inode = i;
      //  point
      xI = (feEngine_->fe_mesh())->nodal_coordinates(i);
      if (!kernelFunction_->node_contributes(xI)) {
        continue;
      }
      for (int j = 0; j < nLocal_; j++) {
        int lammps_j = internalToAtom_(j);
        xa.copy(xPointer_[lammps_j],3);
        // difference vector
        xaI = xa - xI;
        lammpsInterface_->periodicity_correction(xaI.ptr());
        for (int k = 0; k < numneigh[lammps_j]; ++k) {
          int lammps_k = firstneigh[lammps_j][k];
          // first atom location
          xb.copy(xPointer_[lammps_k],3);
          // difference vector
          xba = xb - xa;
          xbI = xba + xaI;
          kernelFunction_->bond_intercepts(xaI,xbI,lam1,lam2);
          // compute virial
          if (lam1 < lam2) {
            bond_value
              = kernel_inv_vol*(kernelFunction_->bond(xaI,xbI,lam1,lam2));
            double delx = xatom[lammps_j][0] - xatom[lammps_k][0];
            double dely = xatom[lammps_j][1] - xatom[lammps_k][1];
            double delz = xatom[lammps_j][2] - xatom[lammps_k][2];
            double rsq = delx*delx + dely*dely + delz*delz;
            double fforce = 0;
            lammpsInterface_->pair_force(lammps_j,lammps_k,rsq,fforce);
            fforce *= 0.5; // dbl count
            double * v = vatom[lammps_j];
            fforce *= (delx*v[0] + dely*v[1] + delz*v[2]);
            if (atomToElementMapType_ == LAGRANGIAN) {
              double delX = xref_[lammps_j][0] - xref_[lammps_k][0];
              double delY = xref_[lammps_j][1] - xref_[lammps_k][1];
              double delZ = xref_[lammps_j][2] - xref_[lammps_k][2];
              flux(inode,0) +=fforce*delX*bond_value;
              flux(inode,1) +=fforce*delY*bond_value;
              flux(inode,2) +=fforce*delZ*bond_value;
            }
            else { // EULERIAN
              flux(inode,0) +=fforce*delx*bond_value;
              flux(inode,1) +=fforce*dely*bond_value;
              flux(inode,2) +=fforce*delz*bond_value;
            }
          }
        }
      }
    }
  }

  //-------------------------------------------------------------------
  // calculation of the dislocation density tensor
  void ATC_TransferKernel::compute_dislocation_density(DENS_MAT & A)
  {
   A.reset(nNodes_,9);
#ifdef HAS_DXA
   double cnaCutoff = lammpsInterface_->near_neighbor_cutoff();
   // Extract dislocation lines within the processor's domain.
   DXADislocationExtractor extractor(lammpsInterface_->lammps_pointer(),dxaExactMode_);
   extractor.extractDislocations(lammpsInterface_->neighbor_list(), cnaCutoff);

   // Calculate scalar dislocation density and density tensor.
   double dislocationDensity = 0.0;
   double dislocationDensityTensor[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


   const std::vector<DislocationSegment*>& segments = extractor.getSegments();
   int localNumberLines = (int) segments.size();
   int totalNumberLines;
   lammpsInterface_->int_allsum(&localNumberLines,&totalNumberLines,1);
   if (totalNumberLines == 0) {
     ATC::LammpsInterface::instance()->print_msg_once("no dislocation lines found");
     return;
   }

   // for output
   int nPt = 0, nSeg = 0;
   for(unsigned segmentIndex = 0; segmentIndex < segments.size(); segmentIndex++) {
     DislocationSegment* segment = segments[segmentIndex];
     const std::deque<Point3>& line = segment->line;
     nPt  += line.size();
     nSeg += line.size()-1;
   }

   DENS_MAT segCoor(3,nPt);
   Array2D<int> segConn(2,nSeg);
   DENS_MAT segBurg(nPt,3);


   DENS_MAT local_A(nNodes_,9);
   local_A.zero();
   DENS_VEC xa(nsd_),xI(nsd_),xaI(nsd_),xb(nsd_),xbI(nsd_),xba(nsd_);
   double kernel_inv_vol = kernelFunction_->inv_vol();
   int iPt = 0, iSeg= 0;
   for(unsigned segmentIndex = 0; segmentIndex < segments.size(); segmentIndex++) {

     DislocationSegment* segment = segments[segmentIndex];
     const std::deque<Point3>& line = segment->line;

     Vector3 burgers = segment->burgersVectorWorld;
     Point3 x1, x2;
     for(std::deque<Point3>::const_iterator p1 = line.begin(), p2 = line.begin() + 1; p2 < line.end(); ++p1, ++p2) {
       x1 = (*p1);
       x2 = (*p2);
       Vector3 delta = x2 - x1;
       // totals
       dislocationDensity += Length(delta);
       for(int i = 0; i < 3; i++) {
         for(int j = 0; j < 3; j++) {
           dislocationDensityTensor[3*j+i] += delta[i] * burgers[j];
         }
       }
       // nodal partition
       for(int k = 0; k < 3; k++) {
         xa(k) = x1[k];
         xb(k) = x2[k];
         xba(k) = delta[k];
       }
       for (int I = 0; I < nNodes_; I++) {
         xI = (feEngine_->fe_mesh())->nodal_coordinates(I);
         if (!kernelFunction_->node_contributes(xI)) {
           continue;
         }
         xaI = xa - xI;
         lammpsInterface_->periodicity_correction(xaI.ptr());
         xbI = xba + xaI;
         double lam1=0,lam2=0;
         kernelFunction_->bond_intercepts(xaI,xbI,lam1,lam2);
         if (lam1 < lam2) {
           double bond_value
             = kernel_inv_vol*(kernelFunction_->bond(xaI,xbI,lam1,lam2));
           local_A(I,0) += xba(0)*burgers[0]*bond_value;
           local_A(I,1) += xba(0)*burgers[1]*bond_value;
           local_A(I,2) += xba(0)*burgers[2]*bond_value;
           local_A(I,3) += xba(1)*burgers[0]*bond_value;
           local_A(I,4) += xba(1)*burgers[1]*bond_value;
           local_A(I,5) += xba(1)*burgers[2]*bond_value;
           local_A(I,6) += xba(2)*burgers[0]*bond_value;
           local_A(I,7) += xba(2)*burgers[1]*bond_value;
           local_A(I,8) += xba(2)*burgers[2]*bond_value;
         }
       }
       segCoor(0,iPt) = x1[0];
       segCoor(1,iPt) = x1[1];
       segCoor(2,iPt) = x1[2];
       segBurg(iPt,0) = burgers[0];
       segBurg(iPt,1) = burgers[1];
       segBurg(iPt,2) = burgers[2];
       segConn(0,iSeg) = iPt;
       segConn(1,iSeg) = iPt+1;
       iPt++;
       iSeg++;
     }
     segCoor(0,iPt) = x2[0];
     segCoor(1,iPt) = x2[1];
     segCoor(2,iPt) = x2[2];
     segBurg(iPt,0) = burgers[0];
     segBurg(iPt,1) = burgers[1];
     segBurg(iPt,2) = burgers[2];
     iPt++;
   }

   int count = nNodes_*9;
   lammpsInterface_->allsum(local_A.ptr(),A.ptr(),count);

   double totalDislocationDensity;
   lammpsInterface_->allsum(&dislocationDensity,&totalDislocationDensity,1);
   double totalDislocationDensityTensor[9];
   lammpsInterface_->allsum(dislocationDensityTensor,totalDislocationDensityTensor,9);
   int totalNumberSegments;
   lammpsInterface_->int_allsum(&nSeg,&totalNumberSegments,1);

   // output
   double volume = lammpsInterface_->domain_volume();
   stringstream ss;
   ss << "total dislocation line length = " << totalDislocationDensity;
   ss << " lines = " << totalNumberLines << " segments = " << totalNumberSegments;
   ss << "\n      ";
   ss << "total dislocation density tensor = \n";
   for(int i = 0; i < 3; i++) {
     ss << "   ";
     for(int j = 0; j < 3; j++) {
       totalDislocationDensityTensor[3*j+i] /= volume;
       ss << totalDislocationDensityTensor[3*j+i] << " ";
     }
     ss << "\n";
   }
   ATC::LammpsInterface::instance()->print_msg_once(ss.str());
   if (nSeg > 0) {
     set<int> otypes;
     otypes.insert(VTK);
     otypes.insert(FULL_GNUPLOT);
     string name = "dislocation_segments_step=" ;
     name += to_string(output_index());
     OutputManager segOutput(name,otypes);
     segOutput.write_geometry(&segCoor,&segConn);
     OUTPUT_LIST segOut;
     segOut["burgers_vector"] = &segBurg;
     segOutput.write_data(0,&segOut);
   }
#else
  throw ATC_Error("unimplemented function compute_dislocation_density (DXA support not included");
#endif
  }

} // end namespace ATC

