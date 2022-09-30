#ifdef PAIR_CLASS
// clang-format off
PairStyle(pod,CPairPOD);
// clang-format on
#else

#ifndef LMP_PAIR_POD_H
#define LMP_PAIR_POD_H

#include "pair.h"
#include "pod.h"

namespace LAMMPS_NS {

  class CPairPOD : public Pair {
    public:
      class CPOD *podptr;
      
      CPairPOD(class LAMMPS *);
      ~CPairPOD() override;
      void compute(int, int) override;
      
      void settings(int, char **) override;
      void coeff(int, char **) override;
      void init_style() override;
      double init_one(int, int) override;
      double memory_usage() override;
      void *extract(const char *, int &) override;
    
      void InitPairPOD(std::string pod_file, std::string coeff_file); 
      
      int backend=1;
      int dim = 3;
      int atommemory = 0;
      int podpairlist=0;
      int lammpspairlist=0;
      int analysis = 0;
      int runMD = 0;
      int savecalculation = 0;
      int savefrequency = 0;
      int blocksize=8192;      // avereage size of computation blocks
      int randomrotation = 0;
      
      double energy=0.0;    // potential energy    
      double *forces=NULL;   // atomic forces
      double *stress=NULL;  // stress tensor
      int *atomtype=NULL;   // atomic types for all atoms
      double *pos=NULL;     // positions of atoms
      double *vel=NULL;     // velocity of atoms
      
      double *gd=NULL;         // global linear descriptors 
      double *podcoeff=NULL;   // POD coefficients  
      double *newpodcoeff=NULL;// normalized POD coefficients  
      double *energycoeff=NULL; // energy coefficients
      double *forcecoeff=NULL; // force coefficients
      
      double *y=NULL;      // [dim * nmaxatom] positions of own and ghost atoms        
      int *atomID=NULL;  // [nmaxatom]  IDs of owned and ghost atoms      
      int *pairlist=NULL;  // [nmaxpairs] indices of neighbors for owned atoms    
      int *pairnum=NULL;   // [nmaxatom] number of pairs for all atoms i 
      int *pairnumsum=NULL;// [nmaxatom+1] cumulative sum of pairnum
      
      int numatoms=0;    // number of atom in the simulation box    
      int nlocalatom=0;  // number of owned atoms    
      int nghostatom=0;  // number of ghost atoms
      int ntotalatom=0;  // number of owned + ghost atoms
      int nlocalmax=0;   // maximum number of owned atoms
      int nmaxatom=0;    // maximum number of atoms (nmaxatom >= ntotalatom)   
      int natompairs=0;  // total number of atom pairs for owned atoms    
      int nmaxpairs=0;   // maximum number of atom pairs for owned atoms    
      
      void check_tempmemory(int start, int end);       
      void check_tempmemory(double **x, int **firstneigh, int *numneigh, int *ilist, int start, int end);
      void estimate_tempmemory();    
      void free_tempmemory();
      void allocate_tempmemory();
      
      void check_atommemory(int inum, int nall);
      void free_atommemory();    
      void allocate_atommemory();    
      
      void free_memory();            
      void allocate_memory();  
                 
      void lammpsNeighPairs(double **x, int **firstneigh, int *atomtype, int *numneigh, int i); 
      
      void printinfo()
      {        
          printf("print pairpod information: \n"); 
          printf("backend %d \n", backend); 
          printf("podpairlist %d \n", podpairlist);         
          printf("lammpspairlist %d \n", lammpspairlist);         
          printf("analysis %d \n", analysis);         
          printf("runMD %d \n", runMD); 
          printf("savecalculation %d \n", savecalculation);
          printf("savefrequency %d \n", savefrequency);       
          printf("numatoms %d \n", numatoms);
          printf("nlocalatom %d \n", nlocalatom);
          printf("nghostatom %d \n", nghostatom);
          printf("ntotalatom %d \n", ntotalatom);
          printf("nlocalmax %d \n", nlocalmax);
          printf("nmaxatom %d \n", nmaxatom);
          printf("natompairs %d \n", natompairs);        
          printf("nmaxpairs %d \n", nmaxpairs);                       
          printf("numblocks %d \n", numblocks);        
          printf("blocksize %d \n", blocksize);
          printf("nablock %d \n", nablock);
          printf("nablockmax %d \n", nablockmax);       
          printf("nij %d \n", nij);       
          printf("nijmax %d \n", nijmax);       
          printf("szd %d \n", szd);       
      }    
          
  protected:
      int atomblocks[1001]; // store sizes of computation blocks   
      int numblocks=1;      // number of computation blocks                        
      int nablock=0;        // number of atoms per computation block    
      int nablockmax=0;     // maximum number of atoms per computation block    
      int nij=0;    //  number of atom pairs per computation block    
      int nijmax=0;  // maximum number of atom pairs per computation block              
      int szd=0;    // size of tmpmem
      
      // temporary arrays for computation blocks  
      double *tmpmem=NULL;  // temporary memory
      int *typeai=NULL;     // types of atoms I only
      int *numneighsum=NULL;// cumulative sum for an array of numbers of neighbors    
      double *rij=NULL;     // (xj - xi) for all pairs (I, J)
      int *idxi=NULL;       // storing linear indices for all pairs (I, J)
      int *ai=NULL;         // IDs of atoms I for all pairs (I, J)
      int *aj=NULL;         // IDs of atoms J for all pairs (I, J)
      int *ti=NULL;         // types of atoms I for all pairs (I, J) 
      int *tj=NULL;         // types of atoms J  for all pairs (I, J)
          
      double **scale;         // for thermodynamic integration        
  };

}    // namespace LAMMPS_NS

#endif
#endif
