// A class for defining transfer operations molecular centers of mass (centroid), dipole moments, quadrupole moments 

#ifndef ATOM_TO_MOLECULE_TRANSFER_H
#define ATOM_TO_MOLECULE_TRANSFER_H

// ATC_Method headers
#include "TransferOperator.h"
#include "MoleculeSet.h"
#include "PaqAtcUtility.h"

using namespace std;

namespace ATC {

  // forward declarations
  class ATC_Method;

  /** 
  * @class PerMoleculeQuantity
  *
  */
  
  template <typename T>
  class PerMoleculeQuantity : public DenseMatrixTransfer<T> {
    
  public:

  PerMoleculeQuantity(ATC_Method * atc):DenseMatrixTransfer<T>(), atc_(atc) {};
    
  virtual ~PerMoleculeQuantity() {};

  protected:

  /** utility object for atc information */
  ATC_Method * atc_;

  private: 
  
  //do not define
  PerMoleculeQuantity();

  };


  /**
  *  @class  AtomtoSmallMoleculeTransfer
  *  @brief  Class for defining objects to transfer total mass or charge of a molecule
  */
  template <typename T>
  class AtomToSmallMoleculeTransfer : public PerMoleculeQuantity<T> {
  
  public:

    //constructor
    // had to define all functions in header, not sure why (JAT, 12/14/11)
    AtomToSmallMoleculeTransfer(ATC_Method * atc, PerAtomQuantity<T> * source, SmallMoleculeSet * smallMoleculeSet)  :
      PerMoleculeQuantity<T>(atc),
      source_(source),
    smallMoleculeSet_(smallMoleculeSet)
    {
      source_->register_dependence(this);
      smallMoleculeSet_->register_dependence(this);
    };
  
  //destructor 
  virtual ~AtomToSmallMoleculeTransfer()
    {
      source_->remove_dependence(this);
      smallMoleculeSet_->remove_dependence(this);
    };

  //  apply transfer operator 
  void reset_quantity() const
  {
    const DenseMatrix<T> & sourceMatrix(source_->quantity());
    int nLocalMol = smallMoleculeSet_->local_molecule_count();
    (this->quantity_).reset(nLocalMol,sourceMatrix.nCols()); 
    for (int i = 0;  i < nLocalMol ; i++) {
      const set<int> & atomsLocalMolArray =  smallMoleculeSet_->atoms_by_local_molecule(i);
      set<int>::const_iterator atomsLocalMolID;
      for (atomsLocalMolID = atomsLocalMolArray.begin(); atomsLocalMolID!= atomsLocalMolArray.end();atomsLocalMolID++) {
        for (int j = 0; j < sourceMatrix.nCols(); j++){
          (this->quantity_)(i,j) += sourceMatrix(*atomsLocalMolID,j);
        }
      }
    }
  };

  protected: 

  // pointer to source atomic quantity data 
  PerAtomQuantity<T> * source_;
  // pointer to molecule data 
  SmallMoleculeSet  * smallMoleculeSet_;

  private:

  // do not define
  AtomToSmallMoleculeTransfer();
  
  };

  /**
   *  @class SmallMoleculeCentroid 
   *  @brief  Class for defining objects to transfer molecular centroid (center of mass)
   */


  class SmallMoleculeCentroid : public AtomToSmallMoleculeTransfer<double> {

  public: 

  //constructor
  SmallMoleculeCentroid(ATC_Method * atc, PerAtomQuantity<double> * source, SmallMoleculeSet * smallMoleculeSet, PerAtomQuantity<double> * atomPositions);

  //destructor
  virtual ~SmallMoleculeCentroid();

  // apply transfer operator  
  virtual void reset_quantity() const;
  
  protected: 
  
  // pointer to source atomic quantity date : positions of atoms in a molecule
  PerAtomQuantity<double> * atomPositions_;

  private: 

  //do not define
  SmallMoleculeCentroid();

  };

  /**
   *  @class SmallMoleculeDipoleMoment
   *  @brief  Class for defining objects to transfer molecular dipole moments 
   */

  class SmallMoleculeDipoleMoment : public SmallMoleculeCentroid {

  public: 
  
  //constructor
  SmallMoleculeDipoleMoment(ATC_Method * atc, PerAtomQuantity<double> * source, SmallMoleculeSet * smallMoleculeSet, PerAtomQuantity<double> * atomPositions, SmallMoleculeCentroid * centroid);
  
  //destructor
  virtual ~SmallMoleculeDipoleMoment();

  // apply transfer operator  
  virtual void reset_quantity() const;

  protected: 

  //pointer to the centroid data 
  SmallMoleculeCentroid * centroid_;

  private: 

  //do not define
  SmallMoleculeDipoleMoment();

};

  /**
   *  @class  AtomToFeTransfer 
   *  @brief  Class for defining objects to transfer molecular quadrupole moments
   */


  class SmallMoleculeQuadrupoleMoment : public SmallMoleculeCentroid {

  public: 

  //constructor
  SmallMoleculeQuadrupoleMoment(ATC_Method * atc, PerAtomQuantity<double> * source, SmallMoleculeSet * smallMoleculeSet, PerAtomQuantity<double> * atomPositions, SmallMoleculeCentroid * centroid);

  //destructor
  virtual ~SmallMoleculeQuadrupoleMoment();

  //apply transfer operator
  virtual void reset_quantity() const;
  
  protected:

  //pointer to the centroid data 
  SmallMoleculeCentroid * centroid_;

  private:

  //do not define
  SmallMoleculeQuadrupoleMoment();

  }; 

    /**
   *  @class  MotfShapeFunctionRestriction
   *  @brief  Class for defining objects that transfer atomistic quantities to FE using shape functions
   *          (implements restrict_volumetric_quantity)
   */

  class MotfShapeFunctionRestriction : public MatToMatTransfer<double> {

  public:
    
    // constructor
    MotfShapeFunctionRestriction(PerMoleculeQuantity<double> * source,
                                 SPAR_MAN * shapeFunction);
    
    // destructor
    virtual ~MotfShapeFunctionRestriction();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** reference to shape function matrix */
    SPAR_MAN * shapeFunction_;

    /** persistant workspace */
    
    
    mutable DENS_MAT _workspace_;

    /** applies restriction operation on this processor */
    virtual void local_restriction(const DENS_MAT & sourceMatrix,
                                   const SPAR_MAT & shapeFunctionMatrix) const;

  private:

    // do not define
    MotfShapeFunctionRestriction();

  };

}
#endif

  
