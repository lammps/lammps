// ATC headers 
#include "ATC_TransferPartitionOfUnity.h"
#include "ATC_Error.h"
#include "FE_Engine.h"
#include "LammpsInterface.h"
#include "Quadrature.h"
#include "PerPairQuantity.h"
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


static const int line_ngauss = 10;
static double line_xg[line_ngauss], line_wg[line_ngauss];

namespace ATC {

  ATC_TransferPartitionOfUnity::ATC_TransferPartitionOfUnity(
    string groupName, 
    double ** & perAtomArray,
    LAMMPS_NS::Fix * thisFix,
    string matParamFile)
    : ATC_Transfer(groupName,perAtomArray,thisFix,matParamFile)
  {
    ATC::Quadrature::instance()->set_line_quadrature(line_ngauss,line_xg,line_wg);
    // transform gauss points from [-1,1] to [0,1]
    double lam1 = 0.0, lam2 = 1.0;
    double del_lambda = 0.5*(lam2 - lam1);
    double avg_lambda = 0.5*(lam2 + lam1);
    for (int i = 0; i < line_ngauss; i++) {
      double lambda = del_lambda*line_xg[i] +avg_lambda;
      line_xg[i] = lambda;
      line_wg[i] *= 0.5;
    }
  }

  //-------------------------------------------------------------------
  ATC_TransferPartitionOfUnity::~ATC_TransferPartitionOfUnity()
  {
    // clear out all managed memory to avoid conflicts with dependencies on class member data
    interscaleManager_.clear();
  }

  //-------------------------------------------------------------------
  void ATC_TransferPartitionOfUnity::compute_projection(
    const DENS_MAT & /* atomData */, DENS_MAT & /* nodeData */)
  {
    throw ATC_Error("unimplemented function");
  }

  //-------------------------------------------------------------------
  void ATC_TransferPartitionOfUnity::compute_bond_matrix()
  {
    atomicBondMatrix_ = bondMatrix_->quantity();
  }

  //-------------------------------------------------------------------
  // kinetic energy portion of stress
  
 /**
 *  @class  KineticTensor
 *  @brief  Class for computing the quantity - m v' (x) v'
 */

  void ATC_TransferPartitionOfUnity::compute_kinetic_stress(
                                      DENS_MAT& stress)
  {
    compute_variation_velocity();
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double mvv2e = lammpsInterface_->mvv2e(); // [MV^2]-->[Energy]

    DENS_MAT & v = variationVelocity_;
 
    atomicTensor_.reset(nLocal_,6);
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      double ma =  mass ? mass[type[atomIdx]]: rmass[atomIdx];
      ma *= mvv2e; // convert mass to appropriate units
      atomicTensor_(i,0) -= ma*v(i,0)*v(i,0);
      atomicTensor_(i,1) -= ma*v(i,1)*v(i,1);
      atomicTensor_(i,2) -= ma*v(i,2)*v(i,2);
      atomicTensor_(i,3) -= ma*v(i,0)*v(i,1);
      atomicTensor_(i,4) -= ma*v(i,0)*v(i,2);
      atomicTensor_(i,5) -= ma*v(i,1)*v(i,2);
    }
    project_volume_normalized(atomicTensor_, stress);
  }

  //-------------------------------------------------------------------
  // on-the-fly calculation of stress
  void ATC_TransferPartitionOfUnity::compute_potential_stress(DENS_MAT& stress)
  {
    int nCols;
    if (atomToElementMapType_ == LAGRANGIAN)
      nCols = 9;
    else // EULERIAN
      nCols = 6;
    stress.reset(nNodes_,nCols);

    // neighbor lists
    int *numneigh = lammpsInterface_->neighbor_list_numneigh();
    int **firstneigh = lammpsInterface_->neighbor_list_firstneigh();
    double ** xatom    = lammpsInterface_->xatom();
    Array<bool> latticePeriodicity(3);
    latticePeriodicity(0) = (bool) periodicity[0];
    latticePeriodicity(1) = (bool) periodicity[1];
    latticePeriodicity(2) = (bool) periodicity[2];
    // process differently for mesh vs translation-invariant kernels
    ATC::LammpsInterface::instance()->stream_msg_once("computing potential stress: ",true,false);
    int heartbeatFreq = (nLocal_ <= 10 ? 1 : (int) nLocal_ / 10);
    // mesh-based kernel functions
    int nodesPerElement = feEngine_->fe_mesh()->num_nodes_per_element();
    Array<int> node_list(nodesPerElement);
    DENS_VEC shp(nodesPerElement);
    DENS_VEC xa(nsd_),xb(nsd_),xab(nsd_),xlambda(nsd_);
    DENS_VEC virial(nCols);
    for (int j = 0; j < nLocal_; j++) {
      if (j % heartbeatFreq == 0 ) {
        ATC::LammpsInterface::instance()->stream_msg_once(".",false,false);
      }
      // first atom location
      int lammps_j = internalToAtom_(j); 
      xa.copy(xPointer_[lammps_j],3);
      for (int k = 0; k < numneigh[lammps_j]; ++k) {
        int lammps_k = firstneigh[lammps_j][k];
        //if (lammps_k < lammps_j) continue; // full neighbor list
        // second (neighbor) atom location
        xb.copy(xPointer_[lammps_k],3);
        double delx = xatom[lammps_j][0] - xatom[lammps_k][0];
        double dely = xatom[lammps_j][1] - xatom[lammps_k][1];
        double delz = xatom[lammps_j][2] - xatom[lammps_k][2];
        double rsq = delx*delx + dely*dely + delz*delz;
        double fforce = 0;
        lammpsInterface_->pair_force(lammps_j,lammps_k,rsq,fforce);
        fforce *= 0.5; // 1/2 sum_ab = sum_(ab)
        if (atomToElementMapType_ == LAGRANGIAN) {
          double delX = xref_[lammps_j][0] - xref_[lammps_k][0]; 
          double delY = xref_[lammps_j][1] - xref_[lammps_k][1]; 
          double delZ = xref_[lammps_j][2] - xref_[lammps_k][2]; 
          virial[0] =-delx*fforce*delX;
          virial[1] =-delx*fforce*delY;
          virial[2] =-delx*fforce*delZ;
          virial[3] =-dely*fforce*delX;
          virial[4] =-dely*fforce*delY;
          virial[5] =-dely*fforce*delZ;
          virial[6] =-delz*fforce*delX;
          virial[7] =-delz*fforce*delY;
          virial[8] =-delz*fforce*delZ;
        }
        else {// EULERIAN
          virial[0] =-delx*delx*fforce;
          virial[1] =-dely*dely*fforce;
          virial[2] =-delz*delz*fforce;
          virial[3] =-delx*dely*fforce;
          virial[4] =-delx*delz*fforce;
          virial[5] =-dely*delz*fforce;
        }
        xab = xa - xb;
        for (int i = 0; i < line_ngauss; i++) {
          double lambda = line_xg[i];
          xlambda = lambda*xab + xb;
          
          lammpsInterface_->periodicity_correction(xlambda.ptr());
          feEngine_->shape_functions(xlambda,shp,node_list);
          // accumulate to nodes whose support overlaps the integration point
          for (int I = 0; I < nodesPerElement; I++) {
            int inode = node_list(I);
            double inv_vol = (accumulantInverseVolumes_->quantity())(inode,inode);
            double bond_value  = inv_vol*shp(I)*line_wg[i];
            for (int j = 0; j < nCols; j++)
              stress(inode,j) += virial[j]*bond_value;
          }
        }
      }
    }
    if (lammpsInterface_->comm_rank() == 0) {
      ATC::LammpsInterface::instance()->stream_msg_once("done",false,true);
    }
  }

  //-------------------------------------------------------------------
  // compute kinetic part of heat flux
  void ATC_TransferPartitionOfUnity::compute_kinetic_heatflux(
                                      DENS_MAT& flux)
  {
    compute_variation_velocity();
    int * type     = lammpsInterface_->atom_type();
    double * mass  = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    double mvv2e = lammpsInterface_->mvv2e();
    double * atomPE = lammpsInterface_->compute_pe_peratom();
    double atomKE, atomEnergy;
    atomicVector_.reset(nLocal_,3);
    for (int i = 0; i < nLocal_; i++) {
      int atomIdx = internalToAtom_(i);
      double ma = mass ? mass[type[atomIdx]]: rmass[atomIdx];
      ma *= mvv2e; // convert mass to appropriate units
      atomKE = 0.0;
      for (int j = 0; j < nsd_; j++) {
        atomKE += 0.5*ma*(variationVelocity_(i,j)*variationVelocity_(i,j));
      }
      atomEnergy = atomKE + atomPE[atomIdx];
      for (int j = 0; j < nsd_; j++) {
        atomicVector_(i,j) += atomEnergy*variationVelocity_(i,j);
      }
    }
    project_volume_normalized(atomicVector_,flux);
  }

  //-------------------------------------------------------------------
  // on-the-fly calculation of the heat flux
  void ATC_TransferPartitionOfUnity::compute_potential_heatflux(DENS_MAT& flux)
  {
    compute_variation_velocity();
    flux.zero();
    // neighbor lists
    int *numneigh = lammpsInterface_->neighbor_list_numneigh();
    int **firstneigh = lammpsInterface_->neighbor_list_firstneigh();
    double ** xatom    = lammpsInterface_->xatom();
    Array<bool> latticePeriodicity(3);
    latticePeriodicity(0) = (bool) periodicity[0];
    latticePeriodicity(1) = (bool) periodicity[1];
    latticePeriodicity(2) = (bool) periodicity[2];
    // process differently for mesh vs translation-invariant kernels
    // mesh-based kernel functions
    int nodesPerElement = feEngine_->fe_mesh()->num_nodes_per_element();
    Array<int> node_list(nodesPerElement);
    DENS_VEC shp(nodesPerElement);
    DENS_VEC xa(nsd_),xb(nsd_),xab(nsd_),xlambda(nsd_);
    for (int j = 0; j < nLocal_; j++) {
      // first atom location
      int lammps_j = internalToAtom_(j); 
      xa.copy(xPointer_[lammps_j],3);
      for (int k = 0; k < numneigh[lammps_j]; ++k) {
        int lammps_k = firstneigh[lammps_j][k];
        // second (neighbor) atom location
        xb.copy(xPointer_[lammps_k],3);
        double delx = xatom[lammps_j][0] - xatom[lammps_k][0];
        double dely = xatom[lammps_j][1] - xatom[lammps_k][1];
        double delz = xatom[lammps_j][2] - xatom[lammps_k][2];
        double rsq = delx*delx + dely*dely + delz*delz;
        double fforce = 0;
        lammpsInterface_->pair_force(lammps_j,lammps_k,rsq,fforce);
        fforce *= 0.5; // 1/2 sum_ab = sum_(ab)
        fforce *= (delx*variationVelocity_(j,0) +
                   dely*variationVelocity_(j,1) +
                   delz*variationVelocity_(j,2));
        double flux_vec[3];
        if (atomToElementMapType_ == LAGRANGIAN) {
          double delX = xref_[lammps_j][0] - xref_[lammps_k][0]; 
          double delY = xref_[lammps_j][1] - xref_[lammps_k][1]; 
          double delZ = xref_[lammps_j][2] - xref_[lammps_k][2]; 
          flux_vec[0] =fforce*delX;
          flux_vec[1] =fforce*delY;
          flux_vec[2] =fforce*delZ;
        }
        else {// EULERIAN
          flux_vec[0] =fforce*delx;
          flux_vec[1] =fforce*dely;
          flux_vec[2] =fforce*delz;
        }
        xab = xa - xb;
        for (int i = 0; i < line_ngauss; i++) {
          double lambda = line_xg[i];
          xlambda = lambda*xab + xb;
          
          lammpsInterface_->periodicity_correction(xlambda.ptr());
          feEngine_->shape_functions(xlambda,shp,node_list);
          // accumulate to nodes whose support overlaps the integration point
          for (int I = 0; I < nodesPerElement; I++) {
            int inode = node_list(I);
            double inv_vol = (accumulantInverseVolumes_->quantity())(inode,inode);
            double bond_value = inv_vol*shp(I)*line_wg[i];
            flux(inode,0) += flux_vec[0]*bond_value;
            flux(inode,1) += flux_vec[1]*bond_value;
            flux(inode,2) += flux_vec[2]*bond_value;
          }
        }
      }
    }
  }

  //-------------------------------------------------------------------
  void ATC_TransferPartitionOfUnity::compute_variation_velocity()
  {
    // now compute v'_a = v_a - N_Ia v_I
    variationVelocity_.reset(nLocal_,nsd_);
    if (nLocal_>0) {
      // interpolate nodal velocities to the atoms
      vbar_.reset(nLocal_,nsd_);
      double ** v    = lammpsInterface_->vatom();
      PerAtomQuantity<double> * vbar = interscaleManager_.per_atom_quantity(field_to_prolongation_name(VELOCITY));
      if (!vbar) {
        DENS_MAN * nodeVelocity = interscaleManager_.dense_matrix(field_to_string(VELOCITY));
        if (this->kernel_on_the_fly()) {
          vbar = new OnTheFlyShapeFunctionProlongation(this,
                   nodeVelocity,this->atom_coarsegraining_positions());
        } else {
          vbar = new FtaShapeFunctionProlongation(this,
                   nodeVelocity,this->interpolant());
        }
        interscaleManager_.add_per_atom_quantity(vbar,
                      field_to_prolongation_name(VELOCITY));
      }
      // use of prolong assumes atom system contained within mesh
      vbar_ = vbar->quantity(); 
      // compute and store variation velocities of atoms
      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        for (int j = 0; j < nsd_; j++) {
          variationVelocity_(i,j) = v[atomIdx][j] - vbar_(i,j);
        }
      }
    }
  }

  //-------------------------------------------------------------------
  // calculation of the dislocation density tensor 
  void ATC_TransferPartitionOfUnity::compute_dislocation_density(DENS_MAT & A)
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
   Array<bool> latticePeriodicity(3);
   latticePeriodicity(0) = (bool) periodicity[0];
   latticePeriodicity(1) = (bool) periodicity[1];
   latticePeriodicity(2) = (bool) periodicity[2];
   // mesh-based kernel functions
   int nodesPerElement = feEngine_->fe_mesh()->num_nodes_per_element();
   Array<int> node_list(nodesPerElement);
   DENS_VEC shp(nodesPerElement);
   DENS_VEC xa(nsd_),xb(nsd_),xba(nsd_),xlambda(nsd_);
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
       for (int i = 0; i < line_ngauss; i++) {
         double lambda = line_xg[i];
         xlambda = lambda*xba + xa;
         
         lammpsInterface_->periodicity_correction(xlambda.ptr());
         feEngine_->shape_functions(xlambda,shp,node_list);
         // accumulate to nodes whose support overlaps the integration point
         for (int I = 0; I < nodesPerElement; I++) {
           int inode = node_list(I);
           double inv_vol = (accumulantInverseVolumes_->quantity())(inode,inode);
           double bond_value  = inv_vol*shp(I)*line_wg[i];
           local_A(inode,0) += xba(0)*burgers[0]*bond_value;
           local_A(inode,1) += xba(0)*burgers[1]*bond_value;
           local_A(inode,2) += xba(0)*burgers[2]*bond_value;
           local_A(inode,3) += xba(1)*burgers[0]*bond_value;
           local_A(inode,4) += xba(1)*burgers[1]*bond_value;
           local_A(inode,5) += xba(1)*burgers[2]*bond_value;
           local_A(inode,6) += xba(2)*burgers[0]*bond_value;
           local_A(inode,7) += xba(2)*burgers[1]*bond_value;
           local_A(inode,8) += xba(2)*burgers[2]*bond_value;
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
   ss.str("");
   DENS_VEC A_avg(9); 
   for (int i = 0; i < nNodes_; i++) {
     for (int j = 0; j < 9; j++) {
       A_avg(j) += A(i,j);
     }
   }
   A_avg /= nNodes_;
   ss << "average nodal dislocation density tensor = \n";
   ss << A_avg(0) << "  " << A_avg(1) << " " << A_avg(2) << "\n"; 
   ss << A_avg(3) << "  " << A_avg(4) << " " << A_avg(5) << "\n"; 
   ss << A_avg(6) << "  " << A_avg(7) << " " << A_avg(8) << "\n"; 
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
    throw ATC_Error("TransferParititionOfUnity::compute_dislocaton_density - unimplemented function");
#endif
  }

} // end namespace ATC

