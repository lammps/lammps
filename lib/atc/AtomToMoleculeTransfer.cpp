// ATC headers 
#include "AtomToMoleculeTransfer.h"
#include "ATC_Method.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class SmallMoleculeCentroid
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  SmallMoleculeCentroid::SmallMoleculeCentroid(ATC_Method * atc, PerAtomQuantity<double> * source, SmallMoleculeSet * smallMoleculeSet, PerAtomQuantity<double> * atomPositions) : AtomToSmallMoleculeTransfer<double>(atc, source, smallMoleculeSet), atomPositions_(atomPositions)
  {
    atomPositions_->register_dependence(this);
  }
 
  //--------------------------------------------------------
  //  Destructor 
  //--------------------------------------------------------
  SmallMoleculeCentroid::~SmallMoleculeCentroid()
  {
    atomPositions_->remove_dependence(this);
  }
  
  //--------------------------------------------------------
  //  Quantity
  //--------------------------------------------------------
  void SmallMoleculeCentroid::reset_quantity() const
  {
    const LammpsInterface * lammps(atc_->lammps_interface());
    const DENS_MAT & sourceMatrix(source_->quantity()); // source is always a scalar quantity here \sum m_i \x_i
    double xi[3], xj[3], xjImage[3];
    const DENS_MAT & atomPosMatrix(atomPositions_->quantity());
    int nLocalMol = smallMoleculeSet_->local_molecule_count();
    quantity_.reset(nLocalMol,atc_->nsd());
    for (int i = 0; i < nLocalMol; i++) {
      const set<int> & atomsLocalMolArray =  smallMoleculeSet_->atoms_by_local_molecule(i);
      set<int>::const_iterator atomsLocalMolID;
      double totalSourceMol = 0.0; // for total source 
      for (atomsLocalMolID = atomsLocalMolArray.begin(); atomsLocalMolID != atomsLocalMolArray.end(); atomsLocalMolID++) {
        totalSourceMol += sourceMatrix(*atomsLocalMolID,0);
      } // compute total source
      atomsLocalMolID = atomsLocalMolArray.begin();
      for (int j = 0; j < atc_->nsd(); j++) {
         xi[j] = atomPosMatrix(*atomsLocalMolID,j);
      }
      for (atomsLocalMolID = atomsLocalMolArray.begin(); atomsLocalMolID != atomsLocalMolArray.end(); atomsLocalMolID++) {

        for (int j = 0; j < atc_->nsd(); j++){
            xj[j] = atomPosMatrix(*atomsLocalMolID,j);
        }
        lammps->closest_image(xi,xj,xjImage);
        for (int j = 0; j < atc_->nsd() ; j++) {
        quantity_(i,j) += sourceMatrix(*atomsLocalMolID,0) * xjImage[j]/totalSourceMol;
        }
      }
    }
  }
      
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class SmallMoleculeDipoleMoment
  //--------------------------------------------------------  
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------

  SmallMoleculeDipoleMoment::SmallMoleculeDipoleMoment(ATC_Method * atc, PerAtomQuantity<double> * source, SmallMoleculeSet * smallMoleculeSet, PerAtomQuantity<double> * atomPositions, SmallMoleculeCentroid * centroid) : SmallMoleculeCentroid(atc, source, smallMoleculeSet, atomPositions), centroid_(centroid) //check here
  {
    centroid_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------

  SmallMoleculeDipoleMoment::~SmallMoleculeDipoleMoment()
  {
    centroid_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  Quantity
  //--------------------------------------------------------
  void SmallMoleculeDipoleMoment::reset_quantity() const
  {
    const LammpsInterface * lammps(atc_->lammps_interface());
    const DENS_MAT & sourceMatrix(source_->quantity()); // source is always a scalar quantity here \sum m_i
    const DENS_MAT & atomPosMatrix(atomPositions_->quantity());
    int nLocalMol = smallMoleculeSet_->local_molecule_count();
    int nsd = atc_->nsd();
    quantity_.reset(nLocalMol,nsd);
    double dx[3];

    //call the SmallMoleculeCentroid here to find Centroid .... 
    const DENS_MAT & centroidMolMatrix(centroid_->quantity()); 
    for (int i = 0; i < nLocalMol; i++) {
      const set<int> & atomsLocalMolArray =  smallMoleculeSet_->atoms_by_local_molecule(i);
      set<int>::const_iterator atomsLocalMolID;;
      for (atomsLocalMolID = atomsLocalMolArray.begin(); atomsLocalMolID != atomsLocalMolArray.end();atomsLocalMolID++) {
        for (int j = 0; j < nsd; j++) {
          dx[j] = atomPosMatrix(*atomsLocalMolID,j) - centroidMolMatrix(i,j);
        }
        lammps->minimum_image(dx[0], dx[1], dx[2]);
        for(int j = 0; j < nsd; j++) {
          quantity_(i,j) += sourceMatrix(*atomsLocalMolID,0) * dx[j]; 
        } 
      }
    }
  }

  //--------------------------------------------------------
  //  Class SmallMoleculeQuadrupoleMoment
  //--------------------------------------------------------  
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------

  SmallMoleculeQuadrupoleMoment::SmallMoleculeQuadrupoleMoment(ATC_Method * atc, PerAtomQuantity<double> * source, SmallMoleculeSet * smallMoleculeSet, PerAtomQuantity<double> * atomPositions, SmallMoleculeCentroid * centroid) : SmallMoleculeCentroid(atc, source, smallMoleculeSet, atomPositions), centroid_(centroid)
  {
    centroid_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------

  SmallMoleculeQuadrupoleMoment::~SmallMoleculeQuadrupoleMoment()
  {
    centroid_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  Quantity
  //--------------------------------------------------------

  void SmallMoleculeQuadrupoleMoment::reset_quantity() const
  {
    const LammpsInterface * lammps(atc_->lammps_interface());
    const DENS_MAT & sourceMatrix(source_->quantity()); // source is always a scalar quantity here \sum m_i
    const DENS_MAT & atomPosMatrix(atomPositions_->quantity());
    int nLocalMol = smallMoleculeSet_->local_molecule_count();
    int nsd = atc_->nsd();
    quantity_.reset(nLocalMol,nsd);
    double dx[3];
    //call the SmallMoleculeCentroid here to find Centroid ....
    const DENS_MAT & centroidMolMatrix(centroid_->quantity());
    for (int i = 0; i < nLocalMol; i++) {
      const set<int> & atomsLocalMolArray =  smallMoleculeSet_->atoms_by_local_molecule(i);
      set<int>::const_iterator atomsLocalMolID;;
      for (atomsLocalMolID = atomsLocalMolArray.begin(); atomsLocalMolID != atomsLocalMolArray.end();atomsLocalMolID++) {
        for (int j = 0; j < nsd; j++) {
          dx[j] = atomPosMatrix(*atomsLocalMolID,j) - centroidMolMatrix(i,j);
        }
        lammps->minimum_image(dx[0], dx[1], dx[2]);
        for(int j = 0; j < nsd; j++) {
          quantity_(i,j) += 0.5*sourceMatrix(*atomsLocalMolID,0) * dx[j] * dx[2];
          /* quantity_(i,3*j)   += 0.5*sourceMatrix(*atomsLocalMolID,0) * dx[j]*dx[0];
          quantity_(i,3*j+1) += 0.5*sourceMatrix(*atomsLocalMolID,0) * dx[j]*dx[1];
          quantity_(i,3*j+2) += 0.5*sourceMatrix(*atomsLocalMolID,0) * dx[j]*dx[2]; */
        }
      }
    }
  }

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  MotfShapeFunctionRestriction::MotfShapeFunctionRestriction(PerMoleculeQuantity<double> * source,
                                                             SPAR_MAN * shapeFunction) :
    MatToMatTransfer<double>(source),
    shapeFunction_(shapeFunction)
  {
    shapeFunction_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  MotfShapeFunctionRestriction::~MotfShapeFunctionRestriction()
  {
    shapeFunction_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void MotfShapeFunctionRestriction::reset_quantity() const
  {
    // computes nodeData = N*atomData where N are the shape functions
    const DENS_MAT & sourceMatrix(source_->quantity());
    // reallocate memory only if sizing has changed
    const SPAR_MAT & shapeFunctionMatrix(shapeFunction_->quantity());
    quantity_.resize(shapeFunctionMatrix.nCols(),sourceMatrix.nCols());
    
    local_restriction(sourceMatrix,shapeFunctionMatrix);
    
    // communicate for total restriction
    int count = quantity_.nRows()*quantity_.nCols();
    lammpsInterface_->allsum(_workspace_.ptr(),quantity_.ptr(),count);
  }

  //--------------------------------------------------------
  //  local_restriction
  //--------------------------------------------------------
  void MotfShapeFunctionRestriction::local_restriction(const DENS_MAT & sourceMatrix,
                                                      const SPAR_MAT & shapeFunctionMatrix) const
  {
    if (sourceMatrix.nRows() > 0)
      _workspace_ = shapeFunctionMatrix.transMat(sourceMatrix);
    else
      _workspace_.reset(quantity_.nRows(),quantity_.nCols());
  }

} // end namespace 
