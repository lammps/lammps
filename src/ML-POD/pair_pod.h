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
  private:
      std::vector<std::string> globVector(const std::string& pattern, std::vector<std::string> & files);

      bool is_a_number(std::string line);

      int latticecoords(double *y, int *alist, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx);

      int podneighborcount(double *r, double rcutsq, int nx, int N, int dim);
      int podneighborlist(int *neighlist, int *numneigh, double *r, double rcutsq, int nx, int N, int dim);

      void get_data(std::vector<std::string> species);

      void read_data_files(std::string data_file, std::vector<std::string> species);

  public:
      struct datastruct {
          std::string file_format;
          std::string file_extension;
          std::string data_path;
          std::vector<std::string> data_files;
          std::vector<std::string> filenames;

          std::vector<int> num_atom;
          std::vector<int> num_atom_cumsum;
          std::vector<int> num_atom_each_file;
          std::vector<int> num_config;
          std::vector<int> num_config_cumsum;
          int num_atom_sum;
          int num_atom_min;
          int num_atom_max;
          int num_config_sum;

          double *lattice=NULL;
          double *energy=NULL;
          double *stress=NULL;
          double *position=NULL;
          double *velocity=NULL;
          double *force=NULL;
          int *atomtype=NULL;

          void copydatainfo(datastruct &data) {
              data.data_path = data_path;
              data.file_format = file_format;
              data.file_extension = file_extension;
              data.data_files = data_files;
              data.filenames = filenames;
          }

          void freememory(int backend)
          {
              TemplateFree(lattice, backend);
              TemplateFree(energy, backend);
              TemplateFree(stress, backend);
              TemplateFree(position, backend);
              TemplateFree(velocity, backend);
              TemplateFree(force, backend);
              TemplateFree(atomtype, backend);
          }
      };

      datastruct data;
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

      void check_pairmemory(double *x, double *a1, double *a2, double *a3, int natom);
      void free_pairmemory();
      void allocate_pairmemory();

      void check_atommemory(int inum, int nall);
      void free_atommemory();
      void allocate_atommemory();

      void free_memory();
      void allocate_memory();

      //void lammpsNeighPairs(double **x, int **firstneigh, int *atomtype, int *numneigh,
      //        int *ilist, int istart, int iend);

      void lammpsNeighPairs(double **x, int **firstneigh, int *atomtype, int *map, int *numneigh, int i);

  protected:
      int nablockmax=0;     // maximum number of atoms per computation block
      int nij=0;    //  number of atom pairs
      int nijmax=0;  // maximum number of atom pairs
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
