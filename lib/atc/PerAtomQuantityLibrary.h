// A library for defining various atomic quantities

#ifndef PER_ATOM_QUANTITY_LIBRARY_H
#define PER_ATOM_QUANTITY_LIBRARY_H

#include "PerAtomQuantity.h"
#include "FundamentalAtomicQuantity.h"
#include <set>
#include <map>
#include <vector>
#include <string>

namespace ATC {

  // forward declarations
  class LammpsInterface;
  class FE_Mesh;
  class FE_Engine;
  template <class T> class DenseMatrixTransfer;

  // need to add capability to take in group bit (JAT, 04/02/11)
  /**
   *  @class  AtomNumber
   *  @brief  Class for identifying atoms based on a specified group
   */

  class AtomNumber : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomNumber(ATC_Method * atc, AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomNumber() {};
    
    /** reset the quantity */
    virtual void reset() const;

  protected:
//  int groupBit_;
    ATC_Method * atc_;

  private:

    // do not define
    AtomNumber();

  };

  /**
   *  @class  AtomTypeVector
   *  @brief  Class for identifying atoms based on a specified group
   */

  class AtomTypeVector : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomTypeVector(ATC_Method * atc, std::vector<int> typeList,
      AtomType atomType = INTERNAL);
    AtomTypeVector(ATC_Method * atc, std::vector<int> typeList, std::vector<int> grpList,
      AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomTypeVector() {};
    
    /** reset the quantity */
    virtual void reset() const;

  protected:
    ATC_Method * atc_;
    int ntypes_;
    std::vector<int> typeList_,index_; // lammps->atc & atc->lammps
    std::vector<int> groupList_; 

  private:
    AtomTypeVector(); // do not define

  };

  
  //      inherited classes are used for this task because
  //      lammps changes pointer location so it can only be
  //      accessed by functions
  /**
   *  @class  XrefWrapper
   *  @brief  Class for wrapping the xref_ array
   */
  
  class XrefWrapper : public ProtectedClonedAtomQuantity<double> {
    
  public:
    
    // constructor
    XrefWrapper(ATC_Method * atc, AtomType atomType=INTERNAL);
    
    // destructor
    virtual ~XrefWrapper() {};
    
  protected:
    
    /** pointer to atc to access raw pointer */
    ATC_Method * atc_;

    /** gets appropriate pointer for lammps data */
    double ** lammps_vector() const;
    
  private:

    // do not define
    XrefWrapper();
    
  };
  
  /**
   *  @class  AtomToElementMap
   *  @brief  Class for identifying the element associated with an atom
   */

  class AtomToElementMap : public ProtectedAtomQuantity<int> {

  public:

    // constructor
    AtomToElementMap(ATC_Method * atc,
                     PerAtomQuantity<double> * atomPositions = NULL,
                     AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomToElementMap();

  protected:

    /** resets the data if necessary */
    virtual void reset() const;

    /** atomic positions */
    PerAtomQuantity<double> * atomPositions_;

  private:

    // do not define
    AtomToElementMap();

  };

  /**
   *  @class  AtomToElementMap
   *  @brief  Class list of atoms in an element set
   */

  class AtomInElementSet : public DependencyManager {

  public:

    // constructor
    AtomInElementSet(ATC_Method * atc,
                     PerAtomQuantity<int> * map,
                     ESET eset, int type); 

    // destructor
    virtual ~AtomInElementSet();

    // accessors
    virtual const ID_LIST & quantity(); 
    virtual ID_LIST & set_quantity() {return list_;}
    int size() {if (needReset_) reset(); return list_.size(); }
    ID_PAIR item(int i) {if (needReset_) reset(); return list_[i]; }
  protected:

    /** resets the data if necessary */
    virtual void reset();

    PaqAtcUtility atc_;

    /** atom to element map */
    PerAtomQuantity<int> * map_;
    ESET eset_;
    int type_;
    const Array<int> & quantityToLammps_;
    ID_LIST list_;

  private:

    // do not define
    AtomInElementSet();

  };

  /**
   *  @class  AtomVolumeUser
   *  @brief  Class for defining the volume per atom based on a user specification
   */

  class AtomVolumeUser : public ProtectedAtomDiagonalMatrix<double> {

  public:

    // constructor
    AtomVolumeUser(ATC_Method * atc,
                   std::map<int,double> & atomGroupVolume,
                   AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomVolumeUser() {};

  protected:

    /** resets the data if necessary */
    virtual void reset() const;

    /** reference to the map of atom group ids to atom volumes */
    std::map<int,double> & atomGroupVolume_;

    /** pointer to lammps interface */
    const LammpsInterface * lammpsInterface_;

    /** map from atc indices to lammps indices */
    const Array<int> & atcToLammps_;

  private:

    // do not define
    AtomVolumeUser();

  };

  /**
   *  @class  AtomVolumeGroup
   *  @brief  Class for defining the volume per atom based on the atom count in a group and its volume
   */

  class AtomVolumeGroup : public AtomVolumeUser {

  public:

    // constructor
    AtomVolumeGroup(ATC_Method * atc,
                    std::map<int,double> & atomGroupVolume,
                    AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomVolumeGroup() {};

  protected:

    /** map from group to group atom volume */
    std::map<int,double> atomGroupVolume_;

    /** pointer to lammps interface */
    const LammpsInterface * lammpsInterface_;
    
    /** reference to array mapping atc indices to lammps indices */
    const Array<int> & atcToLammps_;

  private:

    // do not define
    AtomVolumeGroup();

  };

  /**
   *  @class  AtomVolumeLattice
   *  @brief  Class for defining the volume per atom based on the lattice type and size
   */

  class AtomVolumeLattice : public ProtectedAtomDiagonalMatrix<double> {

  public:

    // constructor
    AtomVolumeLattice(ATC_Method * atc,
                      AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomVolumeLattice() {};

  protected:

    /** resets the data if necessary */
    virtual void reset() const;

    /** pointer to lammps interface */
    const LammpsInterface * lammpsInterface_;

  private:

    // do not define
    AtomVolumeLattice();

  };

  /**
   *  @class  AtomVolumeElement
   *  @brief  Class for defining the volume per atom based on the atom count per element and elemental volume
   */

  class AtomVolumeElement : public ProtectedAtomDiagonalMatrix<double> {

  public:

    // constructor
    AtomVolumeElement(ATC_Method * atc,
                      PerAtomQuantity<int> * atomElement = NULL,
                      AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomVolumeElement();

  protected:

    /** resets the data if necessary */
    virtual void reset() const;

    /** pointer to the atom to element map */
    PerAtomQuantity<int> * atomElement_;

    /** pointer to lammps interface */
    const LammpsInterface * lammpsInterface_;

    /** pointer to mesh object */
    const FE_Mesh * feMesh_;

    // workspace variables
    mutable INT_VECTOR _elementAtomCountLocal_;
    mutable INT_VECTOR _elementAtomCount_;
    mutable DENS_VEC _elementAtomVolume_;
    mutable DENS_MAT _nodalCoordinates_;

  private:

    // do not define
    AtomVolumeElement();

  };

  /**
   *  @class  AtomVolumeRegion
   *  @brief  Class for defining the volume per atom based on the atom count in the MD regions and their volumes.
   *          It will only be meaningful if atoms completely fill all the regions.
   */

  class AtomVolumeRegion : public ProtectedAtomDiagonalMatrix<double> {

  public:

    // constructor
    AtomVolumeRegion(ATC_Method * atc,
                     DENS_MAN * atomCoarseGrainingPositions = NULL,
                     AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomVolumeRegion();

  protected:

    /** resets the data if necessary */
    virtual void reset() const;

    /** pointer to atomic coordinates data */
    DENS_MAN * atomCoarseGrainingPositions_;

    /** pointer to lammps interface */
    const LammpsInterface * lammpsInterface_;

    /** vector from region index to volume */
    DENS_VEC regionalAtomVolume_;

  private:

    // do not define
    AtomVolumeRegion();

  };

  /**
   *  @class  AtomVolumeFile
   *  @brief  Class for defining the volume per atom based on data read in from a file
   */

  class AtomVolumeFile : public ProtectedAtomDiagonalMatrix<double> {

  public:

    // constructor
    AtomVolumeFile(ATC_Method * atc,
                   const std::string & atomVolumeFile,
                   AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomVolumeFile() {};

  protected:

    /** resets the data if necessary */
    virtual void reset() const;

    /** file name containing the atomic information */
    const std::string & atomVolumeFile_;

    /** pointer to lammps interface */
    const LammpsInterface * lammpsInterface_;

  private:

    // do not define
    AtomVolumeFile();

  };

  /**
   *  @class  AtomicMassWeightedDisplacement
   *  @brief  Class for computing the precursor atomic quantity m*(x - x_ref)
   */
  
  class AtomicMassWeightedDisplacement : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomicMassWeightedDisplacement(ATC_Method * atc,
                          PerAtomQuantity<double> * atomPositions = NULL,
                          PerAtomQuantity<double> * atomMasses = NULL,
                          PerAtomQuantity<double> * atomReferencePositions = NULL,
                          AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomicMassWeightedDisplacement();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic positions */
    PerAtomQuantity<double> * atomPositions_;
    
    /** atomic masses */
    PerAtomQuantity<double> * atomMasses_;

    /** atomic reference positions */
    PerAtomQuantity<double> * atomReferencePositions_;

  private:
    
    // do not define
    AtomicMassWeightedDisplacement();

  };

  /**
   *  @class  FluctuatingVelocity 
   *  @brief  Class for computing the atomic quantity v - bar{v}
   */
  
  class FluctuatingVelocity : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    FluctuatingVelocity(ATC_Method * atc,
                       PerAtomQuantity<double> * atomVelocities = NULL,
                       PerAtomQuantity<double> * atomMeanVelocities = NULL,
                       AtomType atomType = INTERNAL);

    // destructor
    virtual ~FluctuatingVelocity();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic velocities */
    PerAtomQuantity<double> * atomVelocities_;

    /** atomic mean velocities */
    PerAtomQuantity<double> * atomMeanVelocities_;
  private:
    
    // do not define
    FluctuatingVelocity();

  };

  /**
   *  @class  ChargeVelcity
   *  @brief  Class for computing the atomic quantity q v'
   */
  
  class ChargeVelocity : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    ChargeVelocity(ATC_Method * atc,
                       PerAtomQuantity<double> * fluctuatingVelocities = NULL,
                       FundamentalAtomQuantity * atomCharges = NULL,
                       AtomType atomType = INTERNAL);

    // destructor
    virtual ~ChargeVelocity();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic velocities */
    PerAtomQuantity<double> * fluctuatingVelocities_;

    /** atomic mean velocities */
    FundamentalAtomQuantity * atomCharge_;
  private:
    
    // do not define
    ChargeVelocity();

  };

  /**
   *  @class  SpeciesVelcity
   *  @brief  Class for computing the atomic quantity m^(a) v'
   */
  
  class SpeciesVelocity : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    SpeciesVelocity(ATC_Method * atc,
                       PerAtomQuantity<double> * fluctuatingVelocities = NULL,
                       PerAtomQuantity<double> * atomTypeVector = NULL,
                       AtomType atomType = INTERNAL);

    // destructor
    virtual ~SpeciesVelocity();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic velocities */
    PerAtomQuantity<double> * fluctuatingVelocities_;

    /** atomic mean velocities */
    PerAtomQuantity<double> * atomTypeVector_;
  private:
    
    // do not define
    SpeciesVelocity();

  };

  /**
   *  @class  AtomicMomentum
   *  @brief  Class for computing the precursor atomic quantity m*v
   */
  
  class AtomicMomentum : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomicMomentum(ATC_Method * atc,
                   PerAtomQuantity<double> * atomVelocities = NULL,
                   PerAtomQuantity<double> * atomMasses = NULL,
                   AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomicMomentum();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic velocities */
    PerAtomQuantity<double> * atomVelocities_;

    /** atomic masses */
    PerAtomQuantity<double> * atomMasses_;

  private:
    
    // do not define
    AtomicMomentum();

  };

  /**
   *  @class  AtomicEnergyForTemperature 
   *  @brief  Base class for accessing quantities needed for computing temperature 
   */
  
  class AtomicEnergyForTemperature : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomicEnergyForTemperature(ATC_Method * atc,
                               AtomType atomType = INTERNAL)
      : ProtectedAtomQuantity<double>(atc, 1, atomType) {};

    // destructor
    virtual ~AtomicEnergyForTemperature() {};

    // returns coefficient which multiplies kinetic energy in temperature definition
    virtual double kinetic_energy_multiplier() const = 0;

  private:
    
    // do not define
    AtomicEnergyForTemperature();

  };

  /**
   *  @class  TwiceKineticEnergy 
   *  @brief  Class for computing the precursor atomic quantity m*v*v 
   *          (used when the kinetic definition of temperature is required)
   */
  
  class TwiceKineticEnergy : public AtomicEnergyForTemperature {

  public:

    // constructor
    TwiceKineticEnergy(ATC_Method * atc,
                       PerAtomQuantity<double> * atomVelocities = NULL,
                       PerAtomQuantity<double> * atomMasses = NULL,
                       AtomType atomType = INTERNAL);

    // destructor
    virtual ~TwiceKineticEnergy();

    // returns coefficient which multiplies kinetic energy in temperature definition
    virtual double kinetic_energy_multiplier() const {return 2.;};

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic velocities */
    PerAtomQuantity<double> * atomVelocities_;

    /** atomic masses */
    PerAtomQuantity<double> * atomMasses_;

  private:
    
    // do not define
    TwiceKineticEnergy();

  };

  /**
   *  @class  KineticTensor 
   *  @brief  Class for computing the atomic quantity m v (x) v 
   */
  
  class KineticTensor : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    KineticTensor(ATC_Method * atc,
                       PerAtomQuantity<double> * atomVelocities = NULL,
                       PerAtomQuantity<double> * atomMasses = NULL,
                       AtomType atomType = INTERNAL);

    // destructor
    virtual ~KineticTensor();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic velocities */
    PerAtomQuantity<double> * atomVelocities_;

    /** atomic masses */
    PerAtomQuantity<double> * atomMasses_;

  private:
    
    // do not define
    KineticTensor();

  };


  /**
   *  @class  FluctuatingKineticTensor 
   *  @brief  Class for computing the atomic quantity m v (x) v 
   */
  
  class FluctuatingKineticTensor : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    FluctuatingKineticTensor(ATC_Method * atc,
                       PerAtomQuantity<double> * atomVelocities = NULL,
                       PerAtomQuantity<double> * atomMasses = NULL,
                       PerAtomQuantity<double> * atomMeanVelocities = NULL,
                       AtomType atomType = INTERNAL);

    // destructor
    virtual ~FluctuatingKineticTensor();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic velocities */
    PerAtomQuantity<double> * atomVelocities_;

    /** atomic masses */
    PerAtomQuantity<double> * atomMasses_;

    /** atomic mean velocities */
    PerAtomQuantity<double> * atomMeanVelocities_;
  private:
    
    // do not define
    FluctuatingKineticTensor();

  };

  /**
   *  @class  TwiceFluctuatingKineticEnergy 
   *  @brief  Class for computing the precursor atomic quantity m*(v-vr)*(v-vr) 
   *          (used when the kinetic definition of temperature is required)
   */
  
  class TwiceFluctuatingKineticEnergy : public AtomicEnergyForTemperature {

  public:

    // constructor
    TwiceFluctuatingKineticEnergy(ATC_Method * atc,
                                  PerAtomQuantity<double> * atomVelocities = NULL,
                                  PerAtomQuantity<double> * atomMasses = NULL,
                                  PerAtomQuantity<double> * atomMeanVelocities = NULL,
                                  AtomType atomType = INTERNAL);

    // destructor
    virtual ~TwiceFluctuatingKineticEnergy();

    // returns coefficient which multiplies kinetic energy in temperature definition
    virtual double kinetic_energy_multiplier() const {return 2.;};

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic velocities */
    PerAtomQuantity<double> * atomVelocities_;

    /** atomic masses */
    PerAtomQuantity<double> * atomMasses_;

    /** atomic mean velocities */
    PerAtomQuantity<double> * atomMeanVelocities_;

  private:
    
    // do not define
    TwiceFluctuatingKineticEnergy();

  };

  /**
   *  @class  MixedKePeEnergy
   *  @brief  Class for computing the precursor atomic quantity for 
   *          a mixed temperature definition involving both KE and PE
   */
  
  class MixedKePeEnergy : public AtomicEnergyForTemperature {

  public:

    // constructor
    MixedKePeEnergy(ATC_Method * atc,
                    double keMultiplier,
                    double peMultiplier,
                    PerAtomQuantity<double> * twiceKineticEnergy = NULL,
                    PerAtomQuantity<double> * potentialEnergy = NULL,
                    AtomType atomType = INTERNAL);

    // destructor
    virtual ~MixedKePeEnergy();

    // returns coefficient which multiplies kinetic energy in temperature definition
    virtual double kinetic_energy_multiplier() const {return keMultiplier_;};

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** factor multiplying kinetic energy */
    double keMultiplier_;

    /** factor multiplying potential energy */
    double peMultiplier_;

    /** twice the kinetic energy of each atom */
    PerAtomQuantity<double> * twiceKineticEnergy_;

    /** potential energy of each atom */
    PerAtomQuantity<double> * potentialEnergy_;

  private:
    
    // do not define
    MixedKePeEnergy();

  };

  /**
   *  @class  TotalEnergy
   *  @brief  Class for the atomic total energy
   */
  
  class TotalEnergy : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    TotalEnergy(ATC_Method * atc,
                PerAtomQuantity<double> * twiceKineticEnergy = NULL,
                PerAtomQuantity<double> * potentialEnergy = NULL,
                AtomType atomType = INTERNAL);

    // destructor
    virtual ~TotalEnergy();

  protected:
    /** handles resetting of data */
    virtual void reset() const;

    /** twice the kinetic energy of each atom */
    PerAtomQuantity<double> * twiceKineticEnergy_;

    /** potential energy of each atom */
    PerAtomQuantity<double> * potentialEnergy_;

  private:
    TotalEnergy(); // do not define
  };

  /**
   *  @class  FluctuatingPotentialEnergy
   *  @brief  Class for computing the precursor atomic quantity for 
   *          a configurational (PE-based) temperature
   */
  
  class FluctuatingPotentialEnergy : public AtomicEnergyForTemperature {

  public:

    // constructor
    FluctuatingPotentialEnergy(ATC_Method * atc,
                               PerAtomQuantity<double> * potentialEnergy = NULL,
                               PerAtomQuantity<double> * referencePotential = NULL,
                               AtomType atomType = INTERNAL);

    // destructor
    virtual ~FluctuatingPotentialEnergy();

    // returns coefficient which multiplies kinetic energy in temperature definition
    virtual double kinetic_energy_multiplier() const {return 0.;;};

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** potential energy of each atom */
    PerAtomQuantity<double> * potentialEnergy_;

    /** twice the kinetic energy of each atom */
    PerAtomQuantity<double> * referencePotential_;

  private:
    
    // do not define
    FluctuatingPotentialEnergy();

  };

  /**
   *  @class  DotTwiceKineticEnergy 
   *  @brief  Class for computing the precursor atomic power 2*v*f 
   *          (used when the kinetic definition of temperature is required)
   */
 
  class DotTwiceKineticEnergy : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    DotTwiceKineticEnergy(ATC_Method * atc,
                          PerAtomQuantity<double> * atomForces = NULL,
                          PerAtomQuantity<double> * atomVelocities = NULL,
                          AtomType atomType = INTERNAL);

    // destructor
    virtual ~DotTwiceKineticEnergy();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic forces */
    PerAtomQuantity<double> * atomForces_;

    /** atomic velocities */
    PerAtomQuantity<double> * atomVelocities_;

  private:
    
    // do not define
    DotTwiceKineticEnergy();

  };

  /**
   *  @class  VelocitySquared
   *  @brief  Class for computing the quantity |v|^2
   *          (used for weights in the thermostat)
   */
 
  class VelocitySquared : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    VelocitySquared(ATC_Method *atc,
                    PerAtomQuantity<double> * atomVelocities = NULL,
                    AtomType atomType = INTERNAL);

    // destructor
    virtual ~VelocitySquared();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic velocities */
    PerAtomQuantity<double> * atomVelocities_;

  private:

    // do not define
    VelocitySquared();

  };

  /**
   *  @class  LambdaSquared
   *  @brief  Class for computing the 2nd order RHS fractional step
   *          contribution to the equation for lambda, with appropriate weights
   */
 
  class LambdaSquared : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    LambdaSquared(ATC_Method *atc,
                  PerAtomQuantity<double> * atomMasses = NULL,
                  PerAtomQuantity<double> * atomVelocitiesSquared = NULL,
                  PerAtomQuantity<double> * atomLambdas = NULL,
                  AtomType atomType = INTERNAL);

    // destructor
    virtual ~LambdaSquared();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic masses */
    PerAtomQuantity<double> * atomMasses_;

    /** atomic velocities squared */
    PerAtomQuantity<double> * atomVelocitiesSquared_;

    /** atomic lambdas */
    PerAtomQuantity<double> * atomLambdas_;

  private:

    // do not define
    LambdaSquared();

  };

  /**
   *  @class  LargeToSmallAtomMap
   *  @brief  mapping from a larger set of atoms to a smaller set
   *          this implementation maximizes storage but reduces execution times,
   *          including taking advantage of MPI communcation
   */

  class LargeToSmallAtomMap : public ProtectedAtomQuantity<int> {

  public:
    
    // constructor
    LargeToSmallAtomMap(ATC_Method * atc,
                        AtomType atomType = INTERNAL)
    : ProtectedAtomQuantity<int>(atc,1,atomType), size_(0) {};
    
    // destructor
    virtual ~LargeToSmallAtomMap() {};

    /** change map when atoms change */
    virtual void reset_nlocal() {this->set_reset();};

    /** get the number of elements associated with the map */
    virtual int size() const {this->quantity(); return size_;};

    /** sets quantity to lammps data, if needed, should be called in pre_exchange */
    virtual void prepare_exchange() {};

    /** sets quantity to lammps data, if needed */
    virtual void post_exchange() {this->set_reset();};

    /** returns how much lammps memory is used in this function */
    virtual int memory_usage() const {return 0;};

    /** packs up data for parallel transfer when atoms change processors */
    virtual int pack_exchange(int i, double *buffer) {return 0;};

    /** unpacks data after parallel transfer when atoms change processors */
    virtual int unpack_exchange(int i, double *buffer) {return 0;};

    /** packs up data for parallel transfer to ghost atoms on other processors */
    virtual int pack_comm(int index, double *buf, 
                          int pbc_flag, int *pbc) {return 0;};

    /** unpacks data after parallel transfer to ghost atoms on other processors */
    virtual int unpack_comm(int index, double *buf) {return 0;};

    /** returns size of per-atom communication */
    virtual int size_comm() const {return 0;};

    /** changes size of temperary lammps storage data if transfer is being used */
    virtual void grow_lammps_array(int nmax, const std::string & tag) {};

    /** rearrange memory of temporary lammps storage data, called from copy_array */
    virtual void copy_lammps_array(int i, int j) {};

  protected:

    /** number of nodes in the map */
    mutable int size_;

  };

  /**
   *  @class  AtomToType
   *  @brief  mapping from all atoms to the subset of atoms of a specified type
   */

  class AtomToType : public LargeToSmallAtomMap {

  public:
    
    // constructor
    AtomToType(ATC_Method * atc,
               int type,
               AtomType atomType = INTERNAL);
    
    // destructor
    virtual ~AtomToType() {};

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** tag for type */
    int type_;

  private:

    // do not define
    AtomToType();

  };

  /**
   *  @class  AtomToGroup
   *  @brief  mapping from all atoms to the subset of atoms of a specified group
   */

  class AtomToGroup : public LargeToSmallAtomMap {

  public:
    
    // constructor
    AtomToGroup(ATC_Method * atc,
                int group,
                AtomType atomType = INTERNAL);
    
    // destructor
    virtual ~AtomToGroup() {};

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** tag for group */
    int group_;

  private:

    // do not define
    AtomToGroup();

  };

  /**
   *  @class  AtomToNodeset
   *  @brief  mapping from all atoms to a subset of nodes
   */

  class AtomToNodeset : public LargeToSmallAtomMap {

  public:
    
    // constructor
    AtomToNodeset(ATC_Method * atc,
                  SetDependencyManager<int> * subsetNodes,
                  PerAtomQuantity<int> * atomElement = NULL,
                  AtomType atomType = INTERNAL);
    
    // destructor
    virtual ~AtomToNodeset() {
      atomElement_->remove_dependence(this);
      subsetNodes_->remove_dependence(this);
    };

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** set of nodes which are being regulated */
    SetDependencyManager<int> * subsetNodes_;

    /** map from atom to element in which it resides */
    PerAtomQuantity<int> * atomElement_;

    /** pointer to the finite element engine */
    const FE_Mesh * feMesh_;

    // workspace
    mutable Array<int> _nodes_; // nodes associated with an element

  private:

    // do not define
    AtomToNodeset();

  };

  /**
   *  @class  AtomToElementset
   *  @brief  mapping from all atoms to a subset of elements
   */

  class AtomToElementset : public LargeToSmallAtomMap {

  public:
    
    // constructor
    AtomToElementset(ATC_Method * atc,
                     MatrixDependencyManager<DenseMatrix, bool> * elementMask,
                     PerAtomQuantity<int> * atomElement = NULL,
                     AtomType atomType = INTERNAL);
    
    // destructor
    virtual ~AtomToElementset();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** set of nodes which are being regulated */
    MatrixDependencyManager<DenseMatrix, bool> * elementMask_;

    /** map from atom to element in which it resides */
    PerAtomQuantity<int> * atomElement_;

    /** pointer to the finite element engine */
    const FE_Mesh * feMesh_;

  private:

    // do not define
    AtomToElementset();

  };

  /**
   *  @class  MappedAtomQuantity
   *  @brief  generic reduced mapping
   */
  
  class MappedAtomQuantity : public ProtectedMappedAtomQuantity<double> {

  public:
    
    // constructor
    MappedAtomQuantity(ATC_Method * atc,
                       PerAtomQuantity<double> * source,
                       LargeToSmallAtomMap * map,
                       AtomType atomType = INTERNAL);
    
    // destructor
    virtual ~MappedAtomQuantity() {
      source_->remove_dependence(this);
      map_->remove_dependence(this);
    };

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** original quantity */
    PerAtomQuantity<double> * source_;

    /** mapping transfer */
    LargeToSmallAtomMap * map_;

  private:

    // do not define
    MappedAtomQuantity();

  };

  /**
   *  @class  VelocitySquaredMapped
   *  @brief  Class for computing the quantity |v|^2 on a subset of atoms
   *          (used for weights in the thermostat)
   */
 
  class VelocitySquaredMapped : public ProtectedMappedAtomQuantity<double> {

  public:

    // constructor
    VelocitySquaredMapped(ATC_Method *atc,
                          MatrixDependencyManager<DenseMatrix, int> * atomMap,
                          PerAtomQuantity<double> * atomVelocities = NULL,
                          AtomType atomType = INTERNAL);

    // destructor
    virtual ~VelocitySquaredMapped();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic velocities */
    PerAtomQuantity<double> * atomVelocities_;

  private:

    // do not define
    VelocitySquaredMapped();

  };

    /**
   *  @class  LambdaSquaredMapped
   *  @brief  Class for computing the 2nd order RHS fractional step
   *          contribution to the equation for lambda, with appropriate weights
   */
 
  class LambdaSquaredMapped : public ProtectedMappedAtomQuantity<double> {

  public:

    // constructor
    LambdaSquaredMapped(ATC_Method *atc,
                        MatrixDependencyManager<DenseMatrix, int> * atomMap,
                        PerAtomQuantity<double> * atomMasses = NULL,
                        PerAtomQuantity<double> * atomVelocitiesSquared = NULL,
                        PerAtomQuantity<double> * atomLambdas = NULL,
                        AtomType atomType = INTERNAL);

    // destructor
    virtual ~LambdaSquaredMapped();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic masses */
    PerAtomQuantity<double> * atomMasses_;

    /** atomic velocities squared */
    PerAtomQuantity<double> * atomVelocitiesSquared_;

    /** atomic lambdas */
    PerAtomQuantity<double> * atomLambdas_;

  private:

    // do not define
    LambdaSquaredMapped();

  };

  /**
   *  @class  HeatCapacity 
   *  @brief  Class for the classical atomic heat capacity 
   */
  
  class HeatCapacity : public ConstantQuantity<double> {

  public:

    // constructor
    HeatCapacity(ATC_Method * atc, AtomType atomType = INTERNAL);

    // destructor
    virtual ~HeatCapacity() {};

  protected:

  private:
    
    // do not define
    HeatCapacity();

  };

  /**
   *  @class  AtomicVelocityRescaleFactor 
   *  @brief  Class for computing the atomic rescaling induced by the rescaling thermostat 
   */
 
  class AtomicVelocityRescaleFactor : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomicVelocityRescaleFactor(ATC_Method * atc,
                                PerAtomQuantity<double> * atomLambdas = NULL,
                                AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomicVelocityRescaleFactor();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic lambdas */
    PerAtomQuantity<double> * atomLambdas_;

  private:
    
    // do not define
    AtomicVelocityRescaleFactor();

  };

  /**
   *  @class  AtomicFluctuatingVelocityRescaled 
   *  @brief  Class for computing the atomic rescaling of the velocity fluctuations by the rescaling thermostat 
   */
 
  class AtomicFluctuatingVelocityRescaled : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomicFluctuatingVelocityRescaled(ATC_Method * atc,
                                      PerAtomQuantity<double> * atomRescaleFactor = NULL,
                                      PerAtomQuantity<double> * atomFluctuatingVelocity = NULL,
                                      AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomicFluctuatingVelocityRescaled();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic rescaling factor */
    PerAtomQuantity<double> * atomRescaleFactor_;

    /** atomic fluctuating velocity */
    PerAtomQuantity<double> * atomFluctuatingVelocity_;

  private:
    
    // do not define
    AtomicFluctuatingVelocityRescaled();

  };

  /**
   *  @class  AtomicCombinedRescaleThermostatError 
   *  @brief  Class for computing the atomic error in the rescaling thermostat when used in combination with a specified streaming velocity 
   */
 
  class AtomicCombinedRescaleThermostatError : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomicCombinedRescaleThermostatError(ATC_Method * atc,
                                         PerAtomQuantity<double> * atomFluctuatingMomentumRescaled = NULL,
                                         PerAtomQuantity<double> * atomMeanVelocity = NULL,
                                         PerAtomQuantity<double> * atomStreamingVelocity = NULL,
                                         PerAtomQuantity<double> * atomMass = NULL,
                                         AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomicCombinedRescaleThermostatError();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic rescaled fluctuating momentum */
    PerAtomQuantity<double> * atomFluctuatingMomentumRescaled_;

    /** atomic mean (prolonged FE) velocity */
    PerAtomQuantity<double> * atomMeanVelocity_;

    /** atomic streaming velocity, as computed by rescaling kinetothermostat */
    PerAtomQuantity<double> * atomStreamingVelocity_;

    /** atomic masses */
    PerAtomQuantity<double> * atomMass_;

  private:
    
    // do not define
    AtomicCombinedRescaleThermostatError();

  };

  /**
   *  @class  AtomicThermostatForce 
   *  @brief  Class for computing the atomic force induced by the GLC-based thermostats 
   */
  
  class AtomicThermostatForce : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomicThermostatForce(ATC_Method * atc,
                          PerAtomQuantity<double> * atomLambdas = NULL,
                          PerAtomQuantity<double> * atomVelocities = NULL,
                          AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomicThermostatForce();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic lambdas */
    PerAtomQuantity<double> * atomLambdas_;

    /** atomic velocities */
    PerAtomQuantity<double> * atomVelocities_;

  private:
    
    // do not define
    AtomicThermostatForce();

  };

  /**
   *  @class  AtomicKinetostatForceDisplacement 
   *  @brief  Class for computing the atomic force induced by the GLC-based kinetostats 
   */
  
  class AtomicKinetostatForceDisplacement : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomicKinetostatForceDisplacement(ATC_Method * atc,
                                      PerAtomQuantity<double> * atomLambda = NULL,
                                      PerAtomQuantity<double> * atomMass = NULL,
                                      AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomicKinetostatForceDisplacement();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** computes the multiplication factor assocaited with the controlled quantity being an integral of the degrees of freedom */
    virtual double time_step_factor(double dt) const {return 1./dt/dt;};

    /** atomic lambdas */
    PerAtomQuantity<double> * atomLambda_;

    /** atomic velocities */
    PerAtomQuantity<double> * atomMass_;

  private:
    
    // do not define
    AtomicKinetostatForceDisplacement();

  };

  /**
   *  @class  AtomicKinetostatForceVelocity
   *  @brief  Class for computing the atomic force induced by the GLC-based kinetostats 
   */
  
  class AtomicKinetostatForceVelocity : public AtomicKinetostatForceDisplacement {

  public:

    // constructor
    AtomicKinetostatForceVelocity(ATC_Method * atc,
                                  PerAtomQuantity<double> * atomLambda = NULL,
                                  PerAtomQuantity<double> * atomMass = NULL,
                                  AtomType atomType = INTERNAL) :
      AtomicKinetostatForceDisplacement(atc,atomLambda,atomMass,atomType) {};

    // destructor
    virtual ~AtomicKinetostatForceVelocity() {};

  protected:

    /** computes the multiplication factor assocaited with the controlled quantity being an integral of the degrees of freedom */
    virtual double time_step_factor(double dt) const {return 1./dt;};

  private:

    // do not define
    AtomicKinetostatForceVelocity();

  };

  /**
   *  @class  AtomicKinetostatForceStress
   *  @brief  Class for computing the atomic force induced by the stress-based kinetostats 
   */
  
  class AtomicKinetostatForceStress : public ProtectedAtomQuantity<double> {

  public:

    // constructor
    AtomicKinetostatForceStress(ATC_Method * atc,
                                PerAtomQuantity<double> * atomLambda = NULL,
                                AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomicKinetostatForceStress();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic lambdas */
    PerAtomQuantity<double> * atomLambda_;

  private:
    
    // do not define
    AtomicKinetostatForceStress();

  };

  /**
   *  @class  PerAtomKernelFunction
   *  @brief  Class for computing the kernel function at each atom location
   */
  
  class PerAtomKernelFunction : public ProtectedAtomSparseMatrix<double> {

  public:

    // constructor
    PerAtomKernelFunction(ATC_Method * atc,
                          PerAtomQuantity<double> * atomPositions = NULL,
                          AtomType atomType = INTERNAL);

    // destructor
    virtual ~PerAtomKernelFunction();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic coarse-graining positions */
    PerAtomQuantity<double> * atomPositions_;

    /** finite element engine */
    const FE_Engine * feEngine_;

  private:
    
    // do not define
    PerAtomKernelFunction();

  };

  /**
   *  @class  PerAtomShapeFunction
   *  @brief  Class for computing the shape function at each atom location
   */
  
  class PerAtomShapeFunction : public ProtectedAtomSparseMatrix<double> {

  public:

    // constructor
    PerAtomShapeFunction(ATC_Method * atc,
                         PerAtomQuantity<double> * atomPositions = NULL,
                         PerAtomQuantity<int> * atomElements = NULL,
                         AtomType atomType = INTERNAL);

    // destructor
    virtual ~PerAtomShapeFunction();

  protected:

    /** handles resetting of data */
    virtual void reset() const;

    /** atomic coarse-graining positions */
    PerAtomQuantity<double> * atomPositions_;

    /** atom to element map */
    PerAtomQuantity<int> * atomElements_;

    /** finite element engine */
    const FE_Engine * feEngine_;

  private:
    
    // do not define
    PerAtomShapeFunction();

  };

  /**
   *  @class  LambdaCouplingMatrix
   *  @brief  constructs the coupling matrix needed to solve for lambda, i.e. N in N^T w N L = b
   */
  
  class LambdaCouplingMatrix : public ProtectedMappedAtomSparseMatrix<double> {

  public:

    // constructor
    LambdaCouplingMatrix(ATC_Method * atc,
                         MatrixDependencyManager<DenseMatrix, int> * nodeToOverlapMap = NULL,
                         SPAR_MAN * shapeFunction = NULL);
    
    // destructor
    virtual ~LambdaCouplingMatrix() {
      shapeFunction_->remove_dependence(this);
      nodeToOverlapMap_->remove_dependence(this);
    };

  protected:

    /** does the actual computation of the quantity */
    virtual void reset() const;

    /** base shape function */
    SPAR_MAN * shapeFunction_;

    /** map from all nodes to regulated ones */
    MatrixDependencyManager<DenseMatrix, int> * nodeToOverlapMap_;

  private:

    // do not define
    LambdaCouplingMatrix();

  };

  /**
   *  @class  LocalLambdaCouplingMatrix
   *  @brief  constructs the coupling matrix needed to solve for lambda, i.e. N in N^T w N L = b
   *          when localization is being used for the constraint
   */
  
  class LocalLambdaCouplingMatrix : public LambdaCouplingMatrix {

  public:

    // constructor
    LocalLambdaCouplingMatrix(ATC_Method * atc,
                              MatrixDependencyManager<DenseMatrix, int> * lambdaAtomMap = NULL,
                              MatrixDependencyManager<DenseMatrix, int> * nodeToOverlapMap = NULL,
                              SPAR_MAN * shapeFunction = NULL);
    
    // destructor
    virtual ~LocalLambdaCouplingMatrix() {
      lambdaAtomMap_->remove_dependence(this);
    };

  protected:

    /** does the actual computation of the quantity */
    virtual void reset() const;

    /** map from all atoms to regulated ones */
    MatrixDependencyManager<DenseMatrix, int> * lambdaAtomMap_;

  private:

    // do not define
    LocalLambdaCouplingMatrix();

  };

  /**
   *  @class  GhostCouplingMatrix
   *  @brief  constructs the modified shape functions used to compute the total forces between ghost and internal atoms
   */
  
  class GhostCouplingMatrix : public LambdaCouplingMatrix {

  public:

    // constructor
    GhostCouplingMatrix(ATC_Method * atc,
                        SPAR_MAN * shapeFunction,
                        SetDependencyManager<int> * subsetNodes,
                        MatrixDependencyManager<DenseMatrix, int> * nodeToOverlapMap = NULL);
    
    // destructor
    virtual ~GhostCouplingMatrix() {
      subsetNodes_->remove_dependence(this);
    };

  protected:

    /** does the actual computation of the quantity */
    virtual void reset() const;

    /** set of nodes which are being regulated */
    SetDependencyManager<int> * subsetNodes_;

    // workspace
    mutable DENS_VEC _activeNodes_; // nodes which are being regulated are 1, otherwise 0
    mutable DENS_VEC _weights_; // required weighting for each shape function row to enforce partition of unity
    mutable DIAG_MAT _weightMatrix_; // diagonal with necessary scaling for partition of unity

  private:

    // do not define
    GhostCouplingMatrix();

  };


}

#endif
