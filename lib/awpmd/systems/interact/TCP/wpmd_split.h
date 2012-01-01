# ifndef WPMD_SPLIT_H
# define WPMD_SPLIT_H
  
/** @file wpmd.h 
    @brief Representation of electrons by multiple wave packets within WPMD */
  
/*s****************************************************************************
 * $Log: wpmd_split.h,v $
 * Revision 1.2  2011/06/11 16:53:55  valuev
 * sync with LAMMPS
 *
 * Revision 1.1  2011/06/10 17:15:07  morozov
 * First Windows project with the correct directory structure
 *
 * Revision 1.17  2011/06/09 22:55:08  valuev
 * norm matrices
 *
 * Revision 1.16  2011/06/07 19:58:42  valuev
 * corrected partitions
 *
 * Revision 1.15  2011/06/07 17:43:00  valuev
 * added Y derivatives
 *
 * Revision 1.14  2011/06/03 08:13:33  valuev
 * added partitions to account for ghost atoms
 *
 * Revision 1.13  2011/06/01 23:45:35  valuev
 * modified for LAMMPS compatibility
 *
 * Revision 1.12  2011/05/28 17:16:22  valuev
 * fixed template<>, some fixes to UHF
 *
 * Revision 1.11  2011/05/27 08:43:52  valuev
 * fixed split packet antisymmetrized version
 *
 * Revision 1.10  2011/05/25 05:23:43  valuev
 * fixed variable transformation for norm matrix
 *
 * Revision 1.9  2011/05/24 19:54:32  valuev
 * fixed sqmatrix::iterator
 *
 * Revision 1.8  2011/05/20 21:39:49  valuev
 * separated norm calculation
 *
 * Revision 1.7  2011/05/14 18:56:19  valuev
 * derivative for ee split interactions
 *
 * Revision 1.6  2011/05/05 08:56:02  valuev
 * working split wp version
 *
 * Revision 1.5  2011/05/04 16:48:52  valuev
 * fixed syntax
 *
 * Revision 1.4  2011/05/04 09:04:48  valuev
 * completed wp_split (except for ee forces)
 *
 * Revision 1.3  2011/04/22 09:54:24  valuev
 * working on split WP
 *
 * Revision 1.1  2011/04/20 08:43:09  valuev
 * started adding split packet WPMD
 *
 *******************************************************************************/

#include "wpmd.h"


class AWPMD_split: public AWPMD {
protected:
  int s_add, spl_add;
public:
  vector<Vector_2> split_c[2]; ///<\en split coefficients for electrons (c_re, c_im)  or (psi,phi) depending on the norm mode
  vector<int> nspl[2]; ///<\en number of wave packets for each electron (size is ne[i])
  
  vector<double> wf_norm[2]; ///<\en norms for each electron
  vector<double> wf_norm_der[2]; ///<\en norm derivative 
  vector<cdouble> ovl_der[2]; ///<\en overlap derivative: \<psi|psi'\> 
  vector<double> E_der[2]; ///<\en energy derivative with respect to {a,b} coordinates


  vector< cdouble > Lh[2]; ///<\en Substitute for L in Hartree case: block matrices 1x(10*nspl[i]) 
  vector< sqmatrix<double> > Normh[2];  ///<\en Substitute for Norm in Hartree case: block matrices

  ///\en resizes all internal arrays according to new electrons (wavepackets) added
  virtual void resize(int flag);


  

public:
  AWPMD_split():s_add(0),spl_add(0){}
  
 
  ///\en Prepares to setup a new system of particles using \ref add_ion(),
  ///    \ref add_electron() and \ref add_split().
  ///    There is no need to call this function when using 
  ///    \ref set_electrons() and \ref set_ions() to setup particles.
  virtual void reset(){
    for(int s=0;s<2;s++){
      split_c[s].clear();
      nspl[s].clear();   
    }
    s_add=0;
    spl_add=0;
    AWPMD::reset();
  }



  ///\en Setup electrons: forms internal wave packet representations.
  ///    If PBCs are used the coords must be within the range [0, cell)
  ///    the \a splits array defines the number of wavepackets required for each electron
  ///    the data for splits should be placed in the corresponding data arrays
  ///    \a c array contains the splits mixing  coefficints 
  ///    \a n is the number of electrons of a given spin component
  ///    Electron velocity v is multiplied by mass to obtain momentum.
  ///    Default mass (-1) means me.
  ///    Electronic charges q are -1 by default (when q=NULL), otherwise the charges are assigned for each split
  int set_electrons(int s, int nel, Vector_3P x, Vector_3P v, double* w, double* pw, Vector_2 *c, int *splits, double mass=-1, double *q=NULL);

  
  ///\en Starts adding new electron: continue with \ref add_split functions.
  int add_electron(int s){
    if(s < 0 || s > 1)
      return LOGERR(-1,fmt("AWPMD_split.add_electron: invaid spin setting (%d)!",s),LINFO);
    s_add=s;
    spl_add=0;
    return ne[s_add];
  }

  ///\en Adds a new split to current electron.
  ///    May change the arguments according to the constraints set.
  ///    \return global id of the wavepacket (starting from 0 for each spin s)
  ///    Electron velocity v is multiplied by mass to obtain momentum.
  ///    Default mass (-1) means me.
  ///    The tags must be nonzero, >0 for the local particle, <0 for ghost particle. 
  ///    Unique particle id  is abs(tag).
  ///    Default tag (0) means inserting the current particle id as local particle. 
  int add_split(Vector_3 &x, Vector_3 &v, double &w, double &pw, Vector_2 &c, double mass=-1, double q=-1., int tag=0);
  
  
  ///\en gets current electronic coordinates, and (optionally) number of wave packets for each electron
  int get_electrons(int spin, Vector_3P x, Vector_3P v, double* w, double* pw, cdouble *c, int *splits=NULL, double mass=-1);


  void eterm_deriv(int ic1,int s1, int c1,int k1,int ic2,int s2, int c2,int j2,cdouble pref,
                              const OverlapDeriv &o,cdouble v,cdouble dv_aj_conj,
                              cdouble dv_ak,cVector_3 dv_bj_conj, cVector_3 dv_bk);
  
  ///\en adds the derivatives of Y for the term v*Y[s](c2,c1)
  void y_deriv(cdouble v,int s, int c2, int c1);
  

  ///\en Calculates block norms an derivatives
  void calc_norms(int flag);
  
  ///\en Prepares force arrays according to \a flag setting for interaction()
  virtual void clear_forces(int flagi,Vector_3P fi, Vector_3P fe_x, 
                    Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c);

  ///\en Calcualtes the overlap between two electrons taking all split WPs into account.
  ///    Norms must be pre-calculated.
  cdouble overlap(int ic1, int s1, int c1,int ic2, int s2, int c2);

  //e same as interaction, but using Hartee factorization (no antisymmetrization)
  int interaction_hartree(int flag=0, Vector_3P fi=NULL, Vector_3P fe_x=NULL, 
                                      Vector_3P fe_p=NULL, double *fe_w=NULL, double *fe_pw=NULL, Vector_2P fe_c=NULL);

  ///\en Calculates interaction in the system of ni ions + electrons 
  /// the electonic subsystem must be previously setup by set_electrons, ionic by set_ions
  /// 0x1   -- give back ion forces \n
  /// 0x2   -- add ion forces to the existing set \n
  /// 0x4   -- calculate electronic forces  \n
  /// 0x8   -- add electronic forces to the existing arrays \n
  /// 0x10  -- calculate internal electronic derivatives only: \n
  ///          will not update electronic force arrays, which may be NULL, \n
  ///          the forces may be obtained then using \ref get_el_forces() for all WPs \n
  ///          or separately for each WP using \ref get_wp_force()
  /// if PBCs are used the coords must be within a range [0, cell)
  int interaction(int flag=0, Vector_3P fi=NULL, Vector_3P fe_x=NULL, 
                            Vector_3P fe_p=NULL, double *fe_w=NULL, double *fe_pw=NULL, Vector_2P fe_c=NULL);

  ///\en Get electronic forcess in the arrays provided, using calculated internal representation
  ///    Valid flag settings are:\n
  ///    0x4   -- overwrite existing forces  \n
  ///    0x8   -- add electronic forces to the existing arrays \n
  void get_el_forces(int flag, Vector_3P fe_x, 
                               Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c);

 
  void get_wp_force(int s, int ispl, Vector_3P fe_x, Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c){
    WavePacket wk=wp[s][ispl];  
    int indw1=8*ispl;
    int indn1=(nvar[s]/10)*8+2*ispl;
    wk.int2phys_der< eq_minus_second >(E_der[s].begin()+indw1,(double *)fe_x,(double *)fe_p,fe_w,fe_pw,1./one_h);
    *fe_c=-Vector_2(E_der[s][indn1],E_der[s][indn1+1]);
  }
};





# endif
