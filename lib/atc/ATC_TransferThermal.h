/** ATC_TransferThermal : a class for atom-continuum transfers &
    control involving heat transport */

#ifndef ATC_TRANSFER_THERMAL_H
#define ATC_TRANSFER_THERMAL_H

// ATC_Transfer headers
#include "ATC_Transfer.h"
#include "Thermostat.h"
#include "TimeIntegrator.h"

// Other headers
#include <map>

namespace ATC {

  class ATC_TransferThermal : public ATC_Transfer {

  public:
  
    // constructor
    ATC_TransferThermal(std::string groupName, std::string matParamFile,
                        ExtrinsicModelType extrinsic = NO_MODEL);

    // destructor
    ~ATC_TransferThermal();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);
  
    /** pre time integration */
    virtual void initialize();

    /** post time integration */
    virtual void finish();

    /** first time substep routines */
    virtual void pre_init_integrate();
    virtual void mid_init_integrate(){};
    virtual void post_init_integrate();

    /** second time substep routine */
    virtual void pre_final_integrate();
    virtual void post_final_integrate();

    /** compute vector for output */
    virtual double compute_vector(int n);

  private:

    /** sets the position/velocity of the ghost atoms */
    virtual void set_ghost_atoms();
    
    /** resets any arrays associated with local atom count */
    virtual void reset_nlocal();

    /** compute the mass matrix components coming from MD integration */
    virtual void compute_md_mass_matrix(FieldName thisField,
                                        map<FieldName,DIAG_MAT> & massMats);

    /** time integration flag and access method */
    TimeIntegrator::TimeIntegrationType integrationType_;
    virtual TimeIntegrator::TimeIntegrationType get_integration_type(){return integrationType_;};
        
    /** fractional step auxilliary storage */
    DENS_MAT delTheta;
    DENS_MAT delThetaV;
    DENS_MAT dot_atomicTemp;
    DENS_MAT dot_atomicTempOld;
    DENS_MAT dot_dot_atomicTemp;
    DENS_MAT dot_dot_atomicTempOld;

    /** physics specific filter initialization */
    void init_filter();

    /** field mask for velocity integration */
    Array2D<bool> temperatureMask_;

    void output();

    /** thermostat manager */
    Thermostat thermostat_;

    // Data for poor man's fractional step
    bool pmfcOn_;
    DENS_MAT oldFieldTemp_;

    // Add in fields for restarting
    virtual void read_restart_data(string fileName_, OUTPUT_LIST & data);
    virtual void write_restart_data(string fileName_, OUTPUT_LIST & data);
    void pack_thermal_fields(OUTPUT_LIST & data);
    
  };

};
#endif
