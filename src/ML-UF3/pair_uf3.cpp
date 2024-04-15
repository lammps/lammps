/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ajinkya Hire (Univ. of Florida),
                        Hendrik Kraß (Univ. of Constance),
                        Matthias Rupp (Luxembourg Institute of Science and Technology),
                        Richard Hennig (Univ of Florida)
---------------------------------------------------------------------- */

#include "pair_uf3.h"

#include "uf3_pair_bspline.h"
#include "uf3_triplet_bspline.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "text_file_reader.h"

#include <cmath>
#include <unordered_map>

namespace LAMMPS_NS{
  struct UF3Impl {

    std::vector<std::vector<std::vector<double>>> n2b_knot, n2b_coeff;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> n3b_knot_matrix;
    std::unordered_map<std::string, std::vector<std::vector<std::vector<double>>>> n3b_coeff_matrix;
    std::vector<std::vector<uf3_pair_bspline>> UFBS2b;
    std::vector<std::vector<std::vector<uf3_triplet_bspline>>> UFBS3b;

  };
}

using namespace LAMMPS_NS;
using MathConst::THIRD;

/* ---------------------------------------------------------------------- */

PairUF3::PairUF3(LAMMPS *lmp) :
    Pair(lmp), setflag_3b(nullptr), knot_spacing_type_2b(nullptr), knot_spacing_type_3b(nullptr),
    cut(nullptr), cut_3b(nullptr), cut_3b_list(nullptr), min_cut_3b(nullptr)
{
  uf3_impl = new UF3Impl;
  single_enable = 1;    // 1 if single() routine exists
  one_coeff = 1;        // 1 if allows only one coeff * * call 
  restartinfo = 0;      // 1 if pair style writes restart info
  maxshort = 10;
  neighshort = nullptr;
  centroidstressflag = CENTROID_AVAIL;
  manybody_flag = 1;
  bsplines_created = 0;

  n2b_knots_array = nullptr;
  n2b_coeff_array = nullptr;
  n2b_knots_array_size = nullptr;
  n2b_coeff_array_size = nullptr;

  map_3b = nullptr;
  n3b_knots_array = nullptr;
  n3b_coeff_array = nullptr;
  n3b_knots_array_size = nullptr;
  n3b_coeff_array_size = nullptr;
}

/* ---------------------------------------------------------------------- */

PairUF3::~PairUF3()
{
  if (copymode) return;
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(knot_spacing_type_2b);
    memory->destroy(n2b_knots_array_size);
    memory->destroy(n2b_coeff_array_size);
    memory->destroy(n2b_knots_array);
    memory->destroy(n2b_coeff_array);

    if (pot_3b) {
      memory->destroy(setflag_3b);
      memory->destroy(cut_3b);
      memory->destroy(cut_3b_list);
      memory->destroy(min_cut_3b);
      memory->destroy(neighshort);
      memory->destroy(knot_spacing_type_3b);
      memory->destroy(map_3b);
      memory->destroy(n3b_knots_array_size);
      memory->destroy(n3b_coeff_array_size);
      memory->destroy(n3b_knots_array);
      memory->destroy(n3b_coeff_array);
    }
  }
  delete uf3_impl;
}

/* ----------------------------------------------------------------------
 *     global settings
 * ---------------------------------------------------------------------- */

void PairUF3::settings(int narg, char **arg)
{

  if (narg != 1)
    error->all(FLERR,
               "Invalid number of arguments for pair_style uf3"
               "  Are you using a 2-body or 2 & 3-body UF potential?");
  nbody_flag = utils::numeric(FLERR, arg[0], true, lmp);
  const int num_of_elements = atom->ntypes;
  if (nbody_flag == 2) {
    pot_3b = false;
    manybody_flag = 0;
    n2body_pot_files = num_of_elements * (num_of_elements + 1) / 2;
    tot_pot_files = n2body_pot_files;
  } else if (nbody_flag == 3) {
    pot_3b = true;
    single_enable = 0;
    n2body_pot_files = num_of_elements * (num_of_elements + 1) / 2;
    n3body_pot_files = num_of_elements * (num_of_elements * (num_of_elements + 1) / 2);
    tot_pot_files = n2body_pot_files + n3body_pot_files;
  } else
    error->all(FLERR, "Pair style uf3 not (yet) implemented for {}-body terms", nbody_flag);

}

/* ----------------------------------------------------------------------
 *    set coeffs for one or more type pairs
 * ---------------------------------------------------------------------- */
void PairUF3::coeff(int narg, char **arg)
{
  if (narg != 3+atom->ntypes)
    error->all(FLERR, "Invalid number of arguments uf3 in pair coeffs.");

  if (!allocated) allocate();

  map_element2type(narg-3, arg+3, false);

  if (comm->me == 0)
    uf3_read_unified_pot_file(arg[2]);
  communicate();
  //if (narg != 3 && narg != 5) error->all(FLERR, "Invalid number of arguments uf3 in pair coeffs.");

  /*int ilo, ihi, jlo, jhi, klo, khi;
  if (narg == 3) {
    utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
    utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo, i); j <= jhi; j++) uf3_read_pot_file(i, j, arg[2]);
    }
  } else if (narg == 5) {
    utils::bounds(FLERR, arg[1], 1, atom->ntypes, ilo, ihi, error);
    utils::bounds(FLERR, arg[2], 1, atom->ntypes, jlo, jhi, error);
    utils::bounds(FLERR, arg[3], 1, atom->ntypes, klo, khi, error);
    if (!utils::strmatch(arg[0], "^3b$"))
      error->all(FLERR, "Pair style uf3 3-body terms require the first argument to be 3b");

    for (int i = ilo; i <= ihi; i++) {
      for (int j = jlo; j <= jhi; j++) {
        for (int k = MAX(klo, jlo); k <= khi; k++) uf3_read_pot_file(i, j, k, arg[4]);
      }
    }
  }*/
}

//Broadcast data read from potential file to all processors
void PairUF3::communicate()
{
  const int num_of_elements = atom->ntypes;
  MPI_Bcast(&cut[0][0], (num_of_elements + 1)*(num_of_elements + 1),
            MPI_DOUBLE, 0, world);

  MPI_Bcast(&n2b_knots_array_size[0][0],
            (num_of_elements + 1)*(num_of_elements + 1), MPI_INT, 0, world);
  MPI_Bcast(&n2b_coeff_array_size[0][0],
            (num_of_elements + 1)*(num_of_elements + 1), MPI_INT, 0, world);
  
  MPI_Bcast(&max_num_knots_2b, 1, MPI_INT, 0, world);
  MPI_Bcast(&max_num_coeff_2b, 1, MPI_INT, 0, world);

  if (pot_3b){
    MPI_Bcast(&cut_3b_list[0][0],
              (num_of_elements + 1)*(num_of_elements + 1), MPI_DOUBLE, 0, world);

    MPI_Bcast(&cut_3b[0][0][0],
              (num_of_elements + 1)*(num_of_elements + 1)*(num_of_elements + 1),
              MPI_DOUBLE, 0, world);

    MPI_Bcast(&n3b_knots_array_size[0][0], tot_interaction_count_3b*3,
              MPI_INT, 0, world);
    MPI_Bcast(&n3b_coeff_array_size[0][0], tot_interaction_count_3b*3,
              MPI_INT, 0, world);

    MPI_Bcast(&max_num_knots_3b, 1, MPI_INT, 0, world);
    MPI_Bcast(&max_num_coeff_3b, 1, MPI_INT, 0, world);
  }

  if (comm->me != 0) {
    memory->destroy(n2b_knots_array);
    memory->destroy(n2b_coeff_array);
    
    memory->create(n2b_knots_array, num_of_elements + 1, num_of_elements + 1,
                   max_num_knots_2b, "pair:n2b_knots_array");
    memory->create(n2b_coeff_array, num_of_elements + 1, num_of_elements + 1,
                   max_num_coeff_2b, "pair:n2b_coeff_array");
    if (pot_3b) {
        memory->destroy(n3b_knots_array);
        memory->destroy(n3b_coeff_array);

        memory->create(n3b_knots_array, tot_interaction_count_3b, 3,
                       max_num_knots_3b, "pair:n3b_knots_array");
    
        memory->create(n3b_coeff_array, tot_interaction_count_3b, max_num_coeff_3b,
                       max_num_coeff_3b, max_num_coeff_3b, "pair:n3b_coeff_array");
    }
  }
  
  MPI_Bcast(&knot_spacing_type_2b[0][0],
            (num_of_elements + 1)*(num_of_elements + 1), MPI_INT, 0, world);

  MPI_Bcast(&n2b_knots_array[0][0][0],
            (num_of_elements + 1)*(num_of_elements + 1)*max_num_knots_2b, MPI_DOUBLE, 0, world);
  MPI_Bcast(&n2b_coeff_array[0][0][0],
            (num_of_elements + 1)*(num_of_elements + 1)*max_num_coeff_2b, MPI_DOUBLE, 0, world);

  MPI_Bcast(&setflag[0][0],
            (num_of_elements + 1)*(num_of_elements + 1), MPI_INT, 0, world);

  if (pot_3b) {
    MPI_Bcast(&knot_spacing_type_3b[0][0][0],
              (num_of_elements + 1)*(num_of_elements + 1)*(num_of_elements + 1),
              MPI_INT, 0, world);
    MPI_Bcast(&n3b_knots_array[0][0][0],
              tot_interaction_count_3b*3*max_num_knots_3b, MPI_DOUBLE, 0, world);
    MPI_Bcast(&n3b_coeff_array[0][0][0][0],
              tot_interaction_count_3b*max_num_coeff_3b*max_num_coeff_3b*max_num_coeff_3b,
              MPI_DOUBLE, 0, world);
    MPI_Bcast(&setflag_3b[0][0][0],
              (num_of_elements + 1)*(num_of_elements + 1)*(num_of_elements + 1),
              MPI_INT, 0, world);
    MPI_Bcast(&min_cut_3b[0][0][0][0],
              (num_of_elements + 1)*(num_of_elements + 1)*(num_of_elements + 1)*3,
              MPI_DOUBLE, 0, world);
  }
}

void PairUF3::allocate()
{
  allocated = 1;
  const int num_of_elements = atom->ntypes;

  map = new int[num_of_elements+1]; //No need to delete map as ~Pair deletes map

  // Contains info about wether UF potential were found for type i and j
  memory->create(setflag, num_of_elements + 1, num_of_elements + 1, "pair:setflag");

  // Contains info about 2-body cutoff distance for type i and j
  // cutsq is the global variable
  // Even though we are making cutsq don't manually change the default values
  // Lammps take care of setting the value
  memory->create(cutsq, num_of_elements + 1, num_of_elements + 1, "pair:cutsq");
  // cut is specific to this pair style. We will set the values in cut
  memory->create(cut, num_of_elements + 1, num_of_elements + 1, "pair:cut");
  //Contains info about type of knot_spacing--> 0 = uniform knot spacing (default)
  //1 = non-uniform knot spacing
  memory->create(knot_spacing_type_2b, num_of_elements + 1, num_of_elements + 1,
                 "pair:knot_spacing_2b");

  //Contains size of 2b knots vectors and 2b coeff matrices
  memory->create(n2b_knots_array_size, num_of_elements + 1, num_of_elements + 1,
                "pair:n2b_knots_array_size");
  memory->create(n2b_coeff_array_size, num_of_elements + 1, num_of_elements + 1,
                "pair:n2b_coeff_array_size");

  // Contains knot_vect of 2-body potential for type i and j
  uf3_impl->n2b_knot.resize(num_of_elements + 1);
  uf3_impl->n2b_coeff.resize(num_of_elements + 1);
  uf3_impl->UFBS2b.resize(num_of_elements + 1);
  for (int i = 1; i < num_of_elements + 1; i++) {
    uf3_impl->n2b_knot[i].resize(num_of_elements + 1);
    uf3_impl->n2b_coeff[i].resize(num_of_elements + 1);
    uf3_impl->UFBS2b[i].resize(num_of_elements + 1);
  }
  if (pot_3b) {
    // Contains info about wether UF potential were found for type i, j and k
    memory->create(setflag_3b, num_of_elements + 1, num_of_elements + 1,
            num_of_elements + 1, "pair:setflag_3b");
    // Contains info about 3-body cutoff distance for type i, j and k
    memory->create(cut_3b, num_of_elements + 1, num_of_elements + 1,
            num_of_elements + 1, "pair:cut_3b");
    // Contains info about 3-body cutoff distance for type i, j and k
    // for constructing 3-body list
    memory->create(cut_3b_list, num_of_elements + 1, num_of_elements + 1,
            "pair:cut_3b_list");
    // Contains info about minimum 3-body cutoff distance for type i, j and k
    memory->create(min_cut_3b, num_of_elements + 1, num_of_elements + 1,
            num_of_elements + 1, 3, "pair:min_cut_3b");
    //Contains info about type of knot_spacing--> 0 = uniform knot spacing (default)
    //1 = non-uniform knot spacing
    memory->create(knot_spacing_type_3b, num_of_elements + 1, num_of_elements + 1,
                   num_of_elements + 1, "pair:knot_spacing_3b");

    tot_interaction_count_3b = 0;
    //conatins map of I-J-K interaction
    memory->create(map_3b, num_of_elements + 1, num_of_elements + 1,
                   num_of_elements + 1, "pair:map_3b");

    // setting cut_3b, setflag = 0 and map_3b
    for (int i = 1; i < num_of_elements + 1; i++) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        cut_3b_list[i][j] = 0;
        for (int k = 1; k < num_of_elements + 1; k++) {
          cut_3b[i][j][k] = 0;
          min_cut_3b[i][j][k][0] = 0;
          min_cut_3b[i][j][k][1] = 0;
          min_cut_3b[i][j][k][2] = 0;

          setflag_3b[i][j][k] = 0;

          map_3b[i][j][k] = tot_interaction_count_3b;
          tot_interaction_count_3b++;
        }
      }
    }

    //contains sizes of 3b knots vectors and 3b coeff matrices
    memory->create(n3b_knots_array_size, tot_interaction_count_3b, 3,
                   "pair:n3b_knots_array_size");
    memory->create(n3b_coeff_array_size, tot_interaction_count_3b, 3,
                   "pair:n3b_coeff_array_size");

    uf3_impl->n3b_knot_matrix.resize(num_of_elements + 1);
    uf3_impl->UFBS3b.resize(num_of_elements + 1);
    for (int i = 1; i < num_of_elements + 1; i++) {
      uf3_impl->n3b_knot_matrix[i].resize(num_of_elements + 1);
      uf3_impl->UFBS3b[i].resize(num_of_elements + 1);
      for (int j = 1; j < num_of_elements + 1; j++) {
        uf3_impl->n3b_knot_matrix[i][j].resize(num_of_elements + 1);
        uf3_impl->UFBS3b[i][j].resize(num_of_elements + 1);
      }
    }
    memory->create(neighshort, maxshort, "pair:neighshort");
    
  }
}

void PairUF3::uf3_read_unified_pot_file(char *potf_name)
{
  //Go through the entire file and get the sizes of knot vectors and
  //coeff vectors/matrices
  //
  //Create arrays
  //
  //Go through the file again and read the knots and coefficients
  //

  const int num_of_elements = atom->ntypes;

  //if (true) { 
  FILE *fp = utils::open_potential(potf_name, lmp, nullptr);
  if (!fp)
    error->all(FLERR,
               "Cannot open UF3 potential file {}: {}",
               potf_name, utils::getsyserror());

  TextFileReader txtfilereader(fp, "UF3:POTFP");
  txtfilereader.ignore_comments = false;

  //while loop over the entire file, find blocks starting with #UF3 POT
  //if block found read the very next line to determine 2B or 3B block
  //if 2B read the knot vector and coeff vector size
  //if 3B read the knot vectors and coeff matrix size
  int line_counter = 1;
  char *line;
  while((line = txtfilereader.next_line(1))){
    Tokenizer line_token(line);

    //Detect start of a block
    if (line_token.contains("#UF3 POT")) {
      //Block start detected
      if (line_token.contains("UNITS:") == 0)
        error->all(FLERR,
                   "UF3: {} file does not contain the 'UNITS:' metadata in "
                   "the header",
                   potf_name);

      //Read the 2nd line of the block
      std::string temp_line = txtfilereader.next_line(1);
      line_counter++;
      ValueTokenizer fp2nd_line(temp_line);

      std::string nbody_on_file = fp2nd_line.next_string();
      if (nbody_on_file == "2B") {
        //2B block
        if (fp2nd_line.count() != 6)
          error->all(FLERR, "UF3: Expected 6 words on line {} of {} file "
                            "but found {} word/s",
                            line_counter, potf_name, fp2nd_line.count());
      
        //get the elements
        std::string element1 = fp2nd_line.next_string();
        std::string element2 = fp2nd_line.next_string();
        int itype = 0;
        int jtype = 0;
        for (int i=1; i<num_of_elements+1; i++){
          if (std::string(elements[map[i]]) == element1) {
            itype = i;
            break;
          }
        }
        for (int i=1; i<num_of_elements+1; i++){
          if (std::string(elements[map[i]]) == element2) {
            jtype = i;
            break;
          }
        }
        
        if ((itype != 0) && (jtype != 0)) {
          //Trailing and leading trim check
          int leading_trim = fp2nd_line.next_int();
          int trailing_trim = fp2nd_line.next_int();
          if (leading_trim != 0)
            error->all(FLERR,
                       "UF3: Current implementation is throughly tested only "
                       "for leading_trim=0");
          if (trailing_trim != 3)
            error->all(FLERR,
                       "UF3: Current implementation is throughly tested only "
                       "for trailing_trim=3");

          //read next line, should contain cutoff and size of knot vector
          temp_line = txtfilereader.next_line(1);
          line_counter++;
          ValueTokenizer fp3rd_line(temp_line);
          if (fp3rd_line.count() != 2)
            error->all(FLERR,
                       "UF3: Expected only 2 words on 3rd line => "
                       "Rij_CUTOFF NUM_OF_KNOTS. Found {} word/s",
                       fp3rd_line.count());

          //cut is used in init_one which is called by pair.cpp at line 267
          //where the return of init_one is squared
          cut[itype][jtype] = fp3rd_line.next_double(); //MPI_Bcast
          cut[jtype][itype] = cut[itype][jtype];        

          int num_knots_2b = fp3rd_line.next_int();
          n2b_knots_array_size[itype][jtype] = num_knots_2b; //MPI_Bcast
          n2b_knots_array_size[jtype][itype] = num_knots_2b; //MPI_Bcast
          max_num_knots_2b = std::max(max_num_knots_2b, num_knots_2b); //MPI_Bcast

          //skip next line
          txtfilereader.skip_line();
          line_counter++;

          //read number of coeff
          temp_line = txtfilereader.next_line(1);
          line_counter++;
          ValueTokenizer fp5th_line(temp_line);
        
          int num_coeff_2b = fp5th_line.next_int();
          n2b_coeff_array_size[itype][jtype] = num_coeff_2b; //MPI_Bcast
          n2b_coeff_array_size[jtype][itype] = num_coeff_2b; //MPI_Bcast
          max_num_coeff_2b = std::max(max_num_coeff_2b, num_coeff_2b); //MPI_Bcast
        }
      }
      else if ((nbody_on_file == "3B") && (pot_3b)) {
        //3B block
        if (fp2nd_line.count() != 7)
          error->all(FLERR, "UF3: Expected 7 words on line {} of {} file"
                            "but found {} word/s",
                            line_counter, potf_name, fp2nd_line.count());
    
        if (nbody_on_file == "3B") {
          //get the elements
          std::string element1 = fp2nd_line.next_string();
          std::string element2 = fp2nd_line.next_string();
          std::string element3 = fp2nd_line.next_string();
          int itype = 0;
          int jtype = 0;
          int ktype = 0;
          for (int i=1; i<num_of_elements+1; i++) {
            if (std::string(elements[map[i]]) == element1) {
              itype = i;
              break;
            }
          }
          for (int i=1; i<num_of_elements+1; i++) {
            if (std::string(elements[map[i]]) == element2) {
              jtype = i;
              break;
            }
          }
          for (int i=1; i<num_of_elements+1; i++) {
            if (std::string(elements[map[i]]) == element3) {
              ktype = i;
              break;
            }
          }
        
          if ((itype != 0) && (jtype != 0) && (ktype!=0)) {
            //Trailing and leading trim check
            int leading_trim = fp2nd_line.next_int();
            int trailing_trim = fp2nd_line.next_int();
            if (leading_trim != 0)
              error->all(FLERR,
                         "UF3: Current implementation is throughly tested "
                         "only for leading_trim=0");
            if (trailing_trim != 3)
              error->all(FLERR,
                         "UF3: Current implementation is throughly tested "
                         "only for trailing_trim=3");

            //read next line, should contain cutoffs and size of knot vectors
            temp_line = txtfilereader.next_line(6);
            line_counter++;
            ValueTokenizer fp3rd_line(temp_line);
        
            if (fp3rd_line.count() != 6)
              error->all(FLERR,
                         "UF3: Expected only 6 numbers on 3rd line => "
                         "Rjk_CUTOFF Rik_CUTOFF Rij_CUTOFF NUM_OF_KNOTS_JK "
                         "NUM_OF_KNOTS_IK NUM_OF_KNOTS_IJ Found {} number/s",
                         fp3rd_line.count());
          
            double cut3b_rjk = fp3rd_line.next_double();
            double cut3b_rij = fp3rd_line.next_double();
            double cut3b_rik = fp3rd_line.next_double();
          
            if (cut3b_rij != cut3b_rik)
              error->all(FLERR,
                         "UF3: rij!=rik for {}-{}-{}. "
                         "Current implementation only works for rij=rik",
                         element1, element2, element3);
          
            if (2 * cut3b_rik != cut3b_rjk)
              error->all(FLERR,
                         "UF3: 2rij=2rik!=rik for {}-{}-{}. "
                         "Current implementation only works for 2rij=2rik!=rik",
                         element1, element2, element3);

            cut_3b_list[itype][jtype] = 
                std::max(cut3b_rij, cut_3b_list[itype][jtype]); //MPI_Bcast
            cut_3b_list[itype][ktype] = 
                std::max(cut_3b_list[itype][ktype], cut3b_rik); //MPI_Bcast

            cut_3b[itype][jtype][ktype] = cut3b_rij; //MPI_Bcast
            cut_3b[itype][ktype][jtype] = cut3b_rik; //MPI_Bcast

            int num_knots_3b_jk = fp3rd_line.next_int();
            int num_knots_3b_ik = fp3rd_line.next_int();
            int num_knots_3b_ij = fp3rd_line.next_int();

            n3b_knots_array_size[map_3b[itype][jtype][ktype]][0] = num_knots_3b_jk; //MPI_Bcast
            n3b_knots_array_size[map_3b[itype][jtype][ktype]][1] = num_knots_3b_ik;
            n3b_knots_array_size[map_3b[itype][jtype][ktype]][2] = num_knots_3b_ij;
            
            n3b_knots_array_size[map_3b[itype][ktype][jtype]][0] = num_knots_3b_jk; //MPI_Bcast
            n3b_knots_array_size[map_3b[itype][ktype][jtype]][1] = num_knots_3b_ij;
            n3b_knots_array_size[map_3b[itype][ktype][jtype]][2] = num_knots_3b_ik;

            max_num_knots_3b = std::max(max_num_knots_3b, num_knots_3b_jk);
            max_num_knots_3b = std::max(max_num_knots_3b, num_knots_3b_ik);
            max_num_knots_3b = std::max(max_num_knots_3b, num_knots_3b_ij); //MPI_Bcast
 
            //skip next 3 line
            txtfilereader.skip_line();
            line_counter++;
            txtfilereader.skip_line();
            line_counter++;
            txtfilereader.skip_line();
            line_counter++;

            //read number of coeff
            temp_line = txtfilereader.next_line(3);
            line_counter++;
            ValueTokenizer fp7th_line(temp_line);

            if (fp7th_line.count() != 3)
              error->all(FLERR,
                         "UF3: Expected 3 numbers on 7th line => "
                         "SHAPE_OF_COEFF_MATRIX[I][J][K] "
                         "found {} numbers",
                         fp7th_line.count());
        
            int coeff_matrix_dim1 = fp7th_line.next_int();
            int coeff_matrix_dim2 = fp7th_line.next_int();
            int coeff_matrix_dim3 = fp7th_line.next_int();

            n3b_coeff_array_size[map_3b[itype][jtype][ktype]][0] = coeff_matrix_dim1; //MPI_Bcast
            n3b_coeff_array_size[map_3b[itype][jtype][ktype]][1] = coeff_matrix_dim2;
            n3b_coeff_array_size[map_3b[itype][jtype][ktype]][2] = coeff_matrix_dim3;

            n3b_coeff_array_size[map_3b[itype][ktype][jtype]][0] = coeff_matrix_dim2;
            n3b_coeff_array_size[map_3b[itype][ktype][jtype]][1] = coeff_matrix_dim1;
            n3b_coeff_array_size[map_3b[itype][ktype][jtype]][2] = coeff_matrix_dim3;
            
            max_num_coeff_3b = std::max(max_num_coeff_3b,coeff_matrix_dim1);
            max_num_coeff_3b = std::max(max_num_coeff_3b,coeff_matrix_dim2);
            max_num_coeff_3b = std::max(max_num_coeff_3b,coeff_matrix_dim3);
          }
        }
      }
      else {
        if (!((nbody_on_file == "3B") && (!pot_3b)))
          error->all(FLERR,
                     "UF3: Expected either '2B' or '3B' word on line {} of {} file",
                     line_counter, potf_name);
      }
    } //if of #UF3 POT
    line_counter++;
  } // while
  //fclose(fp);

  //Create knot and coeff arrays
  if (max_num_knots_2b > 0) {
    //if (comm->me == 0)
    memory->create(n2b_knots_array, num_of_elements + 1, num_of_elements + 1,
                   max_num_knots_2b, "pair:n2b_knots_array");
  }
  else
    error->all(FLERR,
               "UF3: Error reading the size of 2B knot vector\n"
               "Possibly no 2B UF3 potential block detected in {} file",
               potf_name);

  if (max_num_coeff_2b > 0) {
    //if (comm->me == 0)
      /*utils::logmesg(lmp,
                     "max_num_coeff_2b = {}\n", max_num_coeff_2b);*/
    memory->create(n2b_coeff_array, num_of_elements + 1, num_of_elements + 1,
                   max_num_coeff_2b, "pair:n2b_coeff_array");
  }
  else
    error->all(FLERR,
               "UF3: Error reading the size of 2B coeff vector\n"
               "Possibly no 2B UF3 potential block detected in {} file",
               potf_name);

  if (pot_3b) {
    if (max_num_knots_3b > 0)
      memory->create(n3b_knots_array, tot_interaction_count_3b, 3,
                     max_num_knots_3b, "pair:n3b_knots_array");
  
    else
      error->all(FLERR,
                 "UF3: Error reading the size of 3B knot vector\n"
                 "Possibly no 3B UF3 potential block detected in {} file",
                 potf_name);
    
    if (max_num_coeff_3b > 0) 
      memory->create(n3b_coeff_array, tot_interaction_count_3b, max_num_coeff_3b,
                     max_num_coeff_3b, max_num_coeff_3b, "pair:n3b_coeff_array");
    else
      error->all(FLERR,
                 "UF3: Error reading the size of 3B coeff matrices\n"
                 "Possibly no 3B UF3 potential block detected in {} file",
                 potf_name);
  }

  //Go back to the begning of the file 
  txtfilereader.rewind();

  //Go through the file again and fill knot and coeff arrays
  //while loop to read the data
  //if (true) {
  /*FILE *fp = utils::open_potential(potf_name, lmp, nullptr);
  if (!fp)
    error->all(FLERR,
               "Cannot open UF3 potential file {}: {}",
               potf_name, utils::getsyserror());

  TextFileReader txtfilereader(fp, "UF3:POTFP");
  txtfilereader.ignore_comments = false;
  char *line;*/
  while((line = txtfilereader.next_line(1))){
    Tokenizer line_token(line);

    //Detect start of a block
    if (line_token.contains("#UF3 POT")) {
      //Block start detected
      //Read the 2nd line of the block
      std::string temp_line = txtfilereader.next_line(1);
      ValueTokenizer fp2nd_line(temp_line);
      std::string nbody_on_file = fp2nd_line.next_string();

      if (nbody_on_file == "2B") {
        //get the elements
        std::string element1 = fp2nd_line.next_string();
        std::string element2 = fp2nd_line.next_string();
        int itype = 0;
        int jtype = 0;
        for (int i=1; i<num_of_elements+1; i++){
          if (std::string(elements[map[i]]) == element1) {
            itype = i;
            break;
          }
        }
        for (int i=1; i<num_of_elements+1; i++){
          if (std::string(elements[map[i]]) == element2) {
            jtype = i;
            break;
          }
        }

        //skip the next two tokens
        fp2nd_line.skip(2);
       
        //uk or nk?
        std::string knot_type = fp2nd_line.next_string();
        if (knot_type == "uk") {
          knot_spacing_type_2b[itype][jtype] = 0; //MPI_Bcast
          knot_spacing_type_2b[jtype][itype] = 0;
        } else if (knot_type == "nk") {
          knot_spacing_type_2b[itype][jtype] = 1;
          knot_spacing_type_2b[jtype][itype] = 1;
        } else
          error->all(FLERR,
                     "UF3: Expected either 'uk'(uniform-knots) or 'nk'(non-uniform knots). "
                     "Found {} on the 2nd line of {}-{} interaction block",
                      knot_type, element1, element2);

        if ((itype != 0) && (jtype != 0)) {
          //skip line containing info of cutoff and knot vect size
          txtfilereader.skip_line();
          
          int num_knots_2b = n2b_knots_array_size[itype][jtype];

          temp_line = txtfilereader.next_line(num_knots_2b);
          ValueTokenizer fp4th_line(temp_line);

          if (fp4th_line.count() != num_knots_2b)
            error->all(FLERR, "UF3: Error readig the 2B potential block for {}-{}\n"
                       "Expecter {} numbers on 4th line of the block but found {} "
                       "numbers", num_knots_2b,fp4th_line.count());
       
          for (int k = 0; k < num_knots_2b; k++) {
            n2b_knots_array[itype][jtype][k] = fp4th_line.next_double(); //MPI_Bcast
            n2b_knots_array[jtype][itype][k] = n2b_knots_array[itype][jtype][k];
          }

          //skip next line
          txtfilereader.skip_line();

          int num_of_coeff_2b = n2b_coeff_array_size[itype][jtype];

          temp_line = txtfilereader.next_line(num_of_coeff_2b);
          ValueTokenizer fp6th_line(temp_line);

          if (fp6th_line.count() != num_of_coeff_2b)
            error->all(FLERR,
                       "UF3: Error readig the 2B potential block for {}-{}\n"
                       "Expecter {} numbers on 6th line of the block but found {} "
                       "numbers", num_knots_2b,fp4th_line.count());

          for (int k = 0; k < num_of_coeff_2b; k++) {
            n2b_coeff_array[itype][jtype][k] = fp6th_line.next_double(); //MPI_Bcast
            n2b_coeff_array[jtype][itype][k] = n2b_coeff_array[itype][jtype][k];
          }

          if (num_knots_2b != num_of_coeff_2b + 4)
            error->all(FLERR,
                       "UF3: {}-{} interaction block has incorrect knot and "
                       "coeff data nknots!=ncoeffs + 3 + 1",
                       element1, element2);

          setflag[itype][jtype] = 1; //MPI_Bcast
          setflag[jtype][itype] = 1;
        }
      }

      if ((nbody_on_file == "3B") && (pot_3b)) {
        //get the elements
        std::string element1 = fp2nd_line.next_string();
        std::string element2 = fp2nd_line.next_string();
        std::string element3 = fp2nd_line.next_string();
        int itype = 0;
        int jtype = 0;
        int ktype = 0;
        for (int i=1; i<num_of_elements+1; i++) {
          if (std::string(elements[map[i]]) == element1) {
            itype = i;
            break;
          }
        }
        for (int i=1; i<num_of_elements+1; i++) {
          if (std::string(elements[map[i]]) == element2) {
            jtype = i;
            break;
          }
        }
        for (int i=1; i<num_of_elements+1; i++) {
          if (std::string(elements[map[i]]) == element3) {
            ktype = i;
            break;
          }
        }
        
        //skip the next two tokens
        fp2nd_line.skip(2);

        //uk or nk?
        std::string knot_type = fp2nd_line.next_string();
        if (knot_type == "uk") {
          knot_spacing_type_3b[itype][jtype][ktype] = 0; //MPI_Bcast
          knot_spacing_type_3b[itype][ktype][jtype] = 0;
        } else if (knot_type == "nk") {
          knot_spacing_type_3b[itype][jtype][ktype] = 1;
          knot_spacing_type_3b[itype][ktype][jtype] = 1;
        } else
          error->all(FLERR,
                     "UF3: Expected either 'uk'(uniform-knots) or 'nk'(non-uniform knots) "
                     "Found {} on the 2nd line of {}-{}-{} interaction block",
                     knot_type, element1, element2, element3);

        if ((itype != 0) && (jtype != 0) && (ktype!=0)) { 
          //skip line containing info of cutoffs and knot vector sizes
          txtfilereader.skip_line();
          
          int num_knots_3b_jk = n3b_knots_array_size[map_3b[itype][jtype][ktype]][0];
          int num_knots_3b_ik = n3b_knots_array_size[map_3b[itype][jtype][ktype]][1];
          int num_knots_3b_ij = n3b_knots_array_size[map_3b[itype][jtype][ktype]][2];

          temp_line = txtfilereader.next_line(num_knots_3b_jk);
          ValueTokenizer fp4th_line(temp_line);
          if (fp4th_line.count() != num_knots_3b_jk)
            error->all(FLERR,
                       "UF3: Error readig the 3B potential block for {}-{}-{}\n"
                       "Expected {} numbers on 4th line of the block but found {} "
                       "numbers", element1, element2, element3,
                       num_knots_3b_jk, fp4th_line.count());

          for (int i = 0; i < num_knots_3b_jk; i++) {
            n3b_knots_array[map_3b[itype][jtype][ktype]][0][i] =
                fp4th_line.next_double();                     //MPI_Bcast
            n3b_knots_array[map_3b[itype][ktype][jtype]][0][i] =
                n3b_knots_array[map_3b[itype][jtype][ktype]][0][i];
          }

          min_cut_3b[itype][jtype][ktype][0] =
              n3b_knots_array[map_3b[itype][jtype][ktype]][0][0]; //MPI_Bcast
          min_cut_3b[itype][ktype][jtype][0] =
              n3b_knots_array[map_3b[itype][ktype][jtype]][0][0];

          temp_line = txtfilereader.next_line(num_knots_3b_ik);
          ValueTokenizer fp5th_line(temp_line);
          if (fp5th_line.count() != num_knots_3b_ik)
            error->all(FLERR,
                       "UF3: Error readig the 3B potential block for {}-{}-{}\n"
                       "Expected {} numbers on 5th line of the block but found {} "
                       "numbers",  element1, element2, element3,
                       num_knots_3b_ik, fp5th_line.count());

          for (int i = 0; i < num_knots_3b_ik; i++) {
            n3b_knots_array[map_3b[itype][jtype][ktype]][1][i] =
                fp5th_line.next_double();                         //MPI_Bcast
            n3b_knots_array[map_3b[itype][ktype][jtype]][2][i] =
                n3b_knots_array[map_3b[itype][jtype][ktype]][1][i];
          }

          min_cut_3b[itype][jtype][ktype][1] =
              n3b_knots_array[map_3b[itype][jtype][ktype]][1][0];
          min_cut_3b[itype][ktype][jtype][2] =
              n3b_knots_array[map_3b[itype][ktype][jtype]][2][0];

          temp_line = txtfilereader.next_line(num_knots_3b_ij);
          ValueTokenizer fp6th_line(temp_line);
          if (fp6th_line.count() != num_knots_3b_ij)
            error->all(FLERR,
                       "UF3: Error readig the 3B potential block for {}-{}-{}\n"
                       "Expected {} numbers on 6th line of the block but found {} "
                       "numbers",  element1, element2, element3,
                       num_knots_3b_ij, fp6th_line.count());

          for (int i = 0; i < num_knots_3b_ij; i++) {
            n3b_knots_array[map_3b[itype][jtype][ktype]][2][i] =
                fp6th_line.next_double();
            n3b_knots_array[map_3b[itype][ktype][jtype]][1][i] =
                n3b_knots_array[map_3b[itype][jtype][ktype]][2][i];
          }

          min_cut_3b[itype][jtype][ktype][2] = n3b_knots_array[map_3b[itype][jtype][ktype]][2][0];
          min_cut_3b[itype][ktype][jtype][1] = n3b_knots_array[map_3b[itype][ktype][jtype]][1][0];

          //skip next line
          txtfilereader.skip_line();

          int coeff_matrix_dim1 = n3b_coeff_array_size[map_3b[itype][jtype][ktype]][0];
          int coeff_matrix_dim2 = n3b_coeff_array_size[map_3b[itype][jtype][ktype]][1];
          int coeff_matrix_dim3 = n3b_coeff_array_size[map_3b[itype][jtype][ktype]][2];

          if (num_knots_3b_jk != coeff_matrix_dim3 + 3 + 1) 
            error->all(FLERR,
                       "UF3: {}-{}-{} interaction block has incorrect knot "
                       "(NUM_OF_KNOTS_JK) and coeff (coeff_matrix_dim3) data "
                       "nknots!=ncoeffs + 3 + 1",
                        element1, element2, element3);

          if (num_knots_3b_ik != coeff_matrix_dim2 + 3 + 1) 
            error->all(FLERR,
                       "UF3: {}-{}-{} interaction block has incorrect knot "
                       "(NUM_OF_KNOTS_IK) and coeff (coeff_matrix_dim2) data "
                       "nknots!=ncoeffs + 3 + 1",
                        element1, element2, element3);
  
          if (num_knots_3b_ij != coeff_matrix_dim1 + 3 + 1) 
            error->all(FLERR,
                       "UF3: {}-{}-{} interaction block has incorrect knot "
                       "(NUM_OF_KNOTS_IJ) and coeff (coeff_matrix_dim1) data "
                       "nknots!=ncoeffs + 3 + 1",
                        element1, element2, element3);

          int coeff_matrix_elements_len = coeff_matrix_dim3;
          int key1 = map_3b[itype][jtype][ktype];
          int key2 = map_3b[itype][ktype][jtype];

          int line_count = 0;
          for (int i = 0; i < coeff_matrix_dim1; i++) {
            for (int j = 0; j < coeff_matrix_dim2; j++) {
              temp_line = txtfilereader.next_line(coeff_matrix_elements_len);
              ValueTokenizer coeff_line(temp_line);
              if (coeff_line.count() != coeff_matrix_elements_len)
              error->all(FLERR,
                         "UF3: Error reading 3B potential block for {}-{}-{}\n"
                         "Expected {} numbers on {}th line of the block but found {} "
                         "numbers",  element1, element2, element3,
                         coeff_matrix_elements_len, line_count + 8,
                         coeff_line.count());

              for (int k = 0; k < coeff_matrix_dim3; k++) {
                n3b_coeff_array[key1][i][j][k] = coeff_line.next_double();
              }
              line_count += 1;
            }
          }

          for (int i = 0; i < coeff_matrix_dim1; i++) {
            for (int j = 0; j < coeff_matrix_dim2; j++) {
              for (int k = 0; k < coeff_matrix_dim3; k++) {
                n3b_coeff_array[key2][j][i][k] = n3b_coeff_array[key1][i][j][k]; //MPI_Bcast
              }
            }
          }

          setflag_3b[itype][jtype][ktype] = 1; //MPI_Bcast
          setflag_3b[itype][ktype][jtype] = 1;
        }
      }
    } // if #UF3 POT
  } //while
  fclose(fp);

  //Set interaction of atom types of the same elements
  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      if (setflag[i][j] != 1){
        //i-j interaction not set

        //maybe i-j is mapped to some other atom type interaction?
        int i_mapped_to = map[i]+1; //+1 as map starts from 0
        int j_mapped_to = map[j]+1; //+1 as map starts from 0

        if ((i_mapped_to == i) && (j_mapped_to == j))
          //i-j is not mapped to some other atom type ie interaction is missing on file
          error->all(FLERR,
                     "UF3: Potential for interaction {}-{} ie {}-{} not found "
                     "in {} file",
                     i, j, elements[i_mapped_to-1], elements[j_mapped_to-1],
                     potf_name);

        //utils::logmesg(lmp,"Setting stuff for {}-{} mapped to {}-{}\n",i,j,
        //                i_mapped_to, j_mapped_to);

        cut[i][j] = cut[i_mapped_to][j_mapped_to];

        n2b_knots_array_size[i][j] = n2b_knots_array_size[i_mapped_to][j_mapped_to];
        n2b_coeff_array_size[i][j] = n2b_coeff_array_size[i_mapped_to][j_mapped_to];

        knot_spacing_type_2b[i][j] = knot_spacing_type_2b[i_mapped_to][j_mapped_to];

        for (int knot_no = 0; knot_no < max_num_knots_2b; knot_no++)
          n2b_knots_array[i][j][knot_no] = 
              n2b_knots_array[i_mapped_to][j_mapped_to][knot_no];

        for (int coeff_no = 0; coeff_no < max_num_coeff_2b; coeff_no++)
          n2b_coeff_array[i][j][coeff_no] = 
              n2b_coeff_array[i_mapped_to][j_mapped_to][coeff_no];

        setflag[i][j] = 1;
      }
    }
  }

  if (pot_3b) {
    for (int i = 1; i < num_of_elements + 1; i++) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        //for (int k = j; k < num_of_elements + 1; k++) {
        for (int k = 1; k < num_of_elements + 1; k++) {
          if (setflag_3b[i][j][k] != 1) {
            //i-j-k interaction not set
            
            //maybe i-j-k is mapped to some other atom type interaction?
            int i_mapped_to = map[i]+1; //+1 as map starts from 0
            int j_mapped_to = map[j]+1; //+1 as map starts from 0
            int k_mapped_to = map[k]+1; //+1 as map starts from 0

            if ((i_mapped_to == i) && (j_mapped_to == j) && (k_mapped_to == k))
              error->all(FLERR,
                         "UF3: Potential for interaction {}-{}-{} ie {}-{}-{} "
                         " not found in {} file",
                         i, j, k, elements[i_mapped_to-1], elements[j_mapped_to-1],
                         elements[k_mapped_to-1], potf_name);
            if (setflag_3b[i_mapped_to][j_mapped_to][k_mapped_to] != 1)
              error->all(FLERR,
                         "UF3: Interaction {}-{}-{} was mapped to {}-{}-{}, but "
                         "potential interaction for {}-{}-{} was not found in "
                         "{} file",
                         i, j, k, i_mapped_to, j_mapped_to, k_mapped_to,
                         i_mapped_to, j_mapped_to, k_mapped_to, potf_name);
        

            cut_3b_list[i][j] = std::max(cut_3b_list[i_mapped_to][j_mapped_to],
                                         cut_3b_list[i][j]);


            cut_3b[i][j][k] = cut_3b[i_mapped_to][j_mapped_to][k_mapped_to];

            knot_spacing_type_3b[i][j][k] =
                knot_spacing_type_3b[i_mapped_to][j_mapped_to][k_mapped_to];

            int key = map_3b[i][j][k];
            int mapped_to_key = map_3b[i_mapped_to][j_mapped_to][k_mapped_to];

            n3b_knots_array_size[key][0] = n3b_knots_array_size[mapped_to_key][0];
            n3b_knots_array_size[key][1] = n3b_knots_array_size[mapped_to_key][1];
            n3b_knots_array_size[key][2] = n3b_knots_array_size[mapped_to_key][2];

            n3b_coeff_array_size[key][0] = n3b_coeff_array_size[mapped_to_key][0];
            n3b_coeff_array_size[key][1] = n3b_coeff_array_size[mapped_to_key][1];
            n3b_coeff_array_size[key][2] = n3b_coeff_array_size[mapped_to_key][2];

            min_cut_3b[i][j][k][0] = 
                min_cut_3b[i_mapped_to][j_mapped_to][k_mapped_to][0];//n3b_knots_array[key][0][0];
            
            min_cut_3b[i][j][k][1] = 
                min_cut_3b[i_mapped_to][j_mapped_to][k_mapped_to][1];//n3b_knots_array[key][1][0];
            
            min_cut_3b[i][j][k][2] = 
                min_cut_3b[i_mapped_to][j_mapped_to][k_mapped_to][2];//n3b_knots_array[key][2][0];
            
            for (int knot_no = 0; knot_no < n3b_knots_array_size[key][0]; knot_no++)
              n3b_knots_array[key][0][knot_no] = n3b_knots_array[mapped_to_key][0][knot_no];

            for (int knot_no = 0; knot_no < n3b_knots_array_size[key][1]; knot_no++)
              n3b_knots_array[key][1][knot_no] = n3b_knots_array[mapped_to_key][1][knot_no];

            for (int knot_no = 0; knot_no < n3b_knots_array_size[key][2]; knot_no++)
              n3b_knots_array[key][2][knot_no] = n3b_knots_array[mapped_to_key][2][knot_no];

            for (int coeff1 = 0; coeff1 < n3b_coeff_array_size[key][0]; coeff1++)
              for (int coeff2 = 0; coeff2 < n3b_coeff_array_size[key][1]; coeff2++)
                for (int coeff3 = 0; coeff3 < n3b_coeff_array_size[key][2]; coeff3++)
                  n3b_coeff_array[key][coeff1][coeff2][coeff3] = 
                      n3b_coeff_array[mapped_to_key][coeff1][coeff2][coeff3];
            setflag_3b[i][j][k] = 1;
          }
        }
      }
    }
  }
}

void PairUF3::uf3_read_pot_file(int itype, int jtype, char *potf_name)
{
  FILE *fp = utils::open_potential(potf_name, lmp, nullptr);
  if (!fp)
    error->all(FLERR, "Cannot open UF3 potential file {}: {}", potf_name, utils::getsyserror());

  TextFileReader txtfilereader(fp, "UF3:POTFP");
  txtfilereader.ignore_comments = false;

  std::string temp_line = txtfilereader.next_line(1);
  Tokenizer file_header(temp_line);

  if (file_header.count() <= 2)
    error->all(FLERR,
               "UF3: Expected more than two words on 1st line of {} \n"
               "but found {} word/s",
               potf_name, file_header.count());

  if (file_header.contains("#UF3 POT") == 0)
    error->all(FLERR,
               "UF3: {} file is not UF3 POT type, 1st line of UF3 POT \n"
               "files contain '#UF3 POT'. Found {} in the header",
               potf_name, temp_line);

  if (file_header.contains("UNITS:") == 0)
    error->all(FLERR,
               "UF3: {} file does not contain the 'UNITS:' metadata in \n"
               "the header",
               potf_name);

  temp_line = txtfilereader.next_line(1);
  ValueTokenizer fp2nd_line(temp_line);

  if (fp2nd_line.count() != 4)
    error->all(FLERR,
               "UF3: Expected 4 words on 2nd line =>\n"
               "  nBody leading_trim trailing_trim type_of_knot_spacing\n"
               "  Found {}",
               temp_line);

  std::string nbody_on_file = fp2nd_line.next_string();
  if (nbody_on_file != "2B")
    error->all(FLERR, "UF3: Expected a 2B UF3 file but found {}", nbody_on_file);

  int leading_trim = fp2nd_line.next_int();
  int trailing_trim = fp2nd_line.next_int();
  if (leading_trim != 0)
    error->all(FLERR,
               "UF3: Current implementation is throughly tested only for "
               "leading_trim=0");
  if (trailing_trim != 3)
    error->all(FLERR,
               "UF3: Current implementation is throughly tested only for "
               "trailing_trim=3");

  std::string knot_type = fp2nd_line.next_string();
  if (knot_type == "uk") {
    knot_spacing_type_2b[itype][jtype] = 0;
    knot_spacing_type_2b[jtype][itype] = 0;
  } else if (knot_type == "nk") {
    knot_spacing_type_2b[itype][jtype] = 1;
    knot_spacing_type_2b[jtype][itype] = 1;
  } else
    error->all(FLERR,
               "UF3: Expected either 'uk'(uniform-knots) or 'nk'(non-uniform knots). "
               "Found {} on the 2nd line of {} pot file",
               knot_type, potf_name);

  temp_line = txtfilereader.next_line(1);
  ValueTokenizer fp3rd_line(temp_line);
  if (fp3rd_line.count() != 2)
    error->all(FLERR,
               "UF3: Expected only 2 numbers on 3rd line => "
               "Rij_CUTOFF NUM_OF_KNOTS. Found {} number/s",
               fp3rd_line.count());

  //cut is used in init_one which is called by pair.cpp at line 267 where the return of init_one is squared
  cut[itype][jtype] = fp3rd_line.next_double();
  cut[jtype][itype] = cut[itype][jtype];

  int num_knots_2b = fp3rd_line.next_int();

  temp_line = txtfilereader.next_line(num_knots_2b);
  ValueTokenizer fp4th_line(temp_line);

  if (fp4th_line.count() != num_knots_2b)
    error->all(FLERR, "UF3: Expected {} numbers on 4th line but found {} numbers", num_knots_2b,
               fp4th_line.count());

  uf3_impl->n2b_knot[itype][jtype].resize(num_knots_2b);
  uf3_impl->n2b_knot[jtype][itype].resize(num_knots_2b);
  for (int k = 0; k < num_knots_2b; k++) {
    uf3_impl->n2b_knot[itype][jtype][k] = fp4th_line.next_double();
    uf3_impl->n2b_knot[jtype][itype][k] = uf3_impl->n2b_knot[itype][jtype][k];
  }

  temp_line = txtfilereader.next_line(1);
  ValueTokenizer fp5th_line(temp_line);
  int num_of_coeff_2b = fp5th_line.next_int();

  temp_line = txtfilereader.next_line(num_of_coeff_2b);
  ValueTokenizer fp6th_line(temp_line);

  if (fp6th_line.count() != num_of_coeff_2b)
    error->all(FLERR, "UF3: Expected {} numbers on 6th line but found {} numbers", num_of_coeff_2b,
               fp6th_line.count());

  uf3_impl->n2b_coeff[itype][jtype].resize(num_of_coeff_2b);
  uf3_impl->n2b_coeff[jtype][itype].resize(num_of_coeff_2b);
  for (int k = 0; k < num_of_coeff_2b; k++) {
    uf3_impl->n2b_coeff[itype][jtype][k] = fp6th_line.next_double();
    uf3_impl->n2b_coeff[jtype][itype][k] = uf3_impl->n2b_coeff[itype][jtype][k];
  }

  if (uf3_impl->n2b_knot[itype][jtype].size() != uf3_impl->n2b_coeff[itype][jtype].size() + 4) {
    error->all(FLERR, "UF3: {} has incorrect knot and coeff data nknots!=ncoeffs + 3 +1",
               potf_name);
  }
  setflag[itype][jtype] = 1;
  setflag[jtype][itype] = 1;
  fclose(fp);
}

void PairUF3::uf3_read_pot_file(int itype, int jtype, int ktype, char *potf_name)
{
  int coeff_matrix_dim1, coeff_matrix_dim2, coeff_matrix_dim3, coeff_matrix_elements_len;
  FILE *fp = utils::open_potential(potf_name, lmp, nullptr);
  if (!fp)
    error->all(FLERR, "Cannot open UF3 potential file {}: {}", potf_name, utils::getsyserror());

  TextFileReader txtfilereader(fp, "UF3:POTFP");
  txtfilereader.ignore_comments = false;

  std::string temp_line = txtfilereader.next_line(1);
  Tokenizer file_header(temp_line);

  if (file_header.count() <= 2)
    error->all(FLERR,
               "UF3: Expected more than two words on 1st line of {} \n"
               "but found {} word/s",
               potf_name, file_header.count());

  if (file_header.contains("#UF3 POT") == 0)
    error->all(FLERR,
               "UF3: {} file is not UF3 POT type, 1st line of UF3 POT \n"
               "files contain '#UF3 POT'. Found {} in the header",
               potf_name, temp_line);

  if (file_header.contains("UNITS:") == 0)
    error->all(FLERR,
               "UF3: {} file does not contain the 'UNITS:' metadata in \n"
               "the header",
               potf_name);

  temp_line = txtfilereader.next_line(1);
  ValueTokenizer fp2nd_line(temp_line);

  if (fp2nd_line.count() != 4)
    error->all(FLERR,
               "UF3: Expected 3 words on 2nd line => "
               "nBody leading_trim trailing_trim type_of_knot_spacing "
               "Found {}",
               temp_line);

  std::string nbody_on_file = fp2nd_line.next_string();
  if (nbody_on_file != "3B")
    error->all(FLERR, "UF3: Expected a 3B UF3 file but found {}", nbody_on_file);

  int leading_trim = fp2nd_line.next_int();
  int trailing_trim = fp2nd_line.next_int();
  if (leading_trim != 0)
    error->all(FLERR, "UF3: Current implementation is throughly tested only for leading_trim=0");
  if (trailing_trim != 3)
    error->all(FLERR, "UF3: Current implementation is throughly tested only for trailing_trim=3");

  std::string knot_type = fp2nd_line.next_string();
  if (knot_type == "uk") {
    knot_spacing_type_3b[itype][jtype][ktype] = 0;
    knot_spacing_type_3b[itype][ktype][jtype] = 0;
  } else if (knot_type == "nk") {
    knot_spacing_type_3b[itype][jtype][ktype] = 1;
    knot_spacing_type_3b[itype][ktype][jtype] = 1;
  } else
    error->all(FLERR,
               "UF3: Expected either 'uk'(uniform-knots) or 'nk'(non-uniform knots) "
               "Found {} on the 2nd line of {} pot file",
               knot_type, potf_name);

  temp_line = txtfilereader.next_line(6);
  ValueTokenizer fp3rd_line(temp_line);

  if (fp3rd_line.count() != 6)
    error->all(FLERR,
               "UF3: Expected only 6 numbers on 3rd line => "
               "Rjk_CUTOFF Rik_CUTOFF Rij_CUTOFF NUM_OF_KNOTS_JK NUM_OF_KNOTS_IK NUM_OF_KNOTS_IJ "
               "Found {} number/s",
               fp3rd_line.count());

  double cut3b_rjk = fp3rd_line.next_double();
  double cut3b_rij = fp3rd_line.next_double();
  double cut3b_rik = fp3rd_line.next_double();

  if (cut3b_rij != cut3b_rik) {
    error->all(FLERR, "UF3: rij!=rik, Current implementation only works for rij=rik");
  }

  if (2 * cut3b_rik != cut3b_rjk) {
    error->all(FLERR,
               "UF3: 2rij=2rik!=rik, Current implementation only works "
               "for 2rij=2rik!=rik");
  }

  cut_3b_list[itype][jtype] = std::max(cut3b_rij, cut_3b_list[itype][jtype]);
  cut_3b_list[itype][ktype] = std::max(cut_3b_list[itype][ktype], cut3b_rik);

  cut_3b[itype][jtype][ktype] = cut3b_rij;
  cut_3b[itype][ktype][jtype] = cut3b_rik;

  int num_knots_3b_jk = fp3rd_line.next_int();
  temp_line = txtfilereader.next_line(num_knots_3b_jk);
  ValueTokenizer fp4th_line(temp_line);

  if (fp4th_line.count() != num_knots_3b_jk)
    error->all(FLERR, "UF3: Expected {} numbers on 4th line but found {} numbers", num_knots_3b_jk,
               fp4th_line.count());

  uf3_impl->n3b_knot_matrix[itype][jtype][ktype].resize(3);
  uf3_impl->n3b_knot_matrix[itype][ktype][jtype].resize(3);

  uf3_impl->n3b_knot_matrix[itype][jtype][ktype][0].resize(num_knots_3b_jk);
  uf3_impl->n3b_knot_matrix[itype][ktype][jtype][0].resize(num_knots_3b_jk);

  for (int i = 0; i < num_knots_3b_jk; i++) {
    uf3_impl->n3b_knot_matrix[itype][jtype][ktype][0][i] = fp4th_line.next_double();
    uf3_impl->n3b_knot_matrix[itype][ktype][jtype][0][i] = uf3_impl->n3b_knot_matrix[itype][jtype][ktype][0][i];
  }

  min_cut_3b[itype][jtype][ktype][0] = uf3_impl->n3b_knot_matrix[itype][jtype][ktype][0][0];
  min_cut_3b[itype][ktype][jtype][0] = uf3_impl->n3b_knot_matrix[itype][ktype][jtype][0][0];

  int num_knots_3b_ik = fp3rd_line.next_int();
  temp_line = txtfilereader.next_line(num_knots_3b_ik);
  ValueTokenizer fp5th_line(temp_line);

  if (fp5th_line.count() != num_knots_3b_ik)
    error->all(FLERR, "UF3: Expected {} numbers on 5th line but found {} numbers", num_knots_3b_ik,
               fp5th_line.count());

  uf3_impl->n3b_knot_matrix[itype][jtype][ktype][1].resize(num_knots_3b_ik);
  uf3_impl->n3b_knot_matrix[itype][ktype][jtype][2].resize(num_knots_3b_ik);
  for (int i = 0; i < num_knots_3b_ik; i++) {
    uf3_impl->n3b_knot_matrix[itype][jtype][ktype][1][i] = fp5th_line.next_double();
    uf3_impl->n3b_knot_matrix[itype][ktype][jtype][2][i] = uf3_impl->n3b_knot_matrix[itype][jtype][ktype][1][i];
  }

  min_cut_3b[itype][jtype][ktype][1] = uf3_impl->n3b_knot_matrix[itype][jtype][ktype][1][0];
  min_cut_3b[itype][ktype][jtype][2] = uf3_impl->n3b_knot_matrix[itype][ktype][jtype][2][0];

  int num_knots_3b_ij = fp3rd_line.next_int();
  temp_line = txtfilereader.next_line(num_knots_3b_ij);
  ValueTokenizer fp6th_line(temp_line);

  if (fp6th_line.count() != num_knots_3b_ij)
    error->all(FLERR, "UF3: Expected {} numbers on 6th line but found {} numbers", num_knots_3b_ij,
               fp5th_line.count());

  uf3_impl->n3b_knot_matrix[itype][jtype][ktype][2].resize(num_knots_3b_ij);
  uf3_impl->n3b_knot_matrix[itype][ktype][jtype][1].resize(num_knots_3b_ij);
  for (int i = 0; i < num_knots_3b_ij; i++) {
    uf3_impl->n3b_knot_matrix[itype][jtype][ktype][2][i] = fp6th_line.next_double();
    uf3_impl->n3b_knot_matrix[itype][ktype][jtype][1][i] = uf3_impl->n3b_knot_matrix[itype][jtype][ktype][2][i];
  }

  min_cut_3b[itype][jtype][ktype][2] = uf3_impl->n3b_knot_matrix[itype][jtype][ktype][2][0];
  min_cut_3b[itype][ktype][jtype][1] = uf3_impl->n3b_knot_matrix[itype][ktype][jtype][1][0];

  temp_line = txtfilereader.next_line(3);
  ValueTokenizer fp7th_line(temp_line);

  if (fp7th_line.count() != 3)
    error->all(FLERR,
               "UF3: Expected 3 numbers on 7th line => "
               "SHAPE_OF_COEFF_MATRIX[I][J][K] "
               "found {} numbers",
               fp7th_line.count());

  coeff_matrix_dim1 = fp7th_line.next_int();
  coeff_matrix_dim2 = fp7th_line.next_int();
  coeff_matrix_dim3 = fp7th_line.next_int();

  if (uf3_impl->n3b_knot_matrix[itype][jtype][ktype][0].size() != coeff_matrix_dim3 + 3 + 1)
    error->all(FLERR,
               "UF3: {} has incorrect knot (NUM_OF_KNOTS_JK) and "
               "coeff (coeff_matrix_dim3) data nknots!=ncoeffs + 3 +1",
               potf_name);

  if (uf3_impl->n3b_knot_matrix[itype][jtype][ktype][1].size() != coeff_matrix_dim2 + 3 + 1)
    error->all(FLERR,
               "UF3: {} has incorrect knot (NUM_OF_KNOTS_IK) and "
               "coeff (coeff_matrix_dim2) data nknots!=ncoeffs + 3 +1",
               potf_name);

  if (uf3_impl->n3b_knot_matrix[itype][jtype][ktype][2].size() != coeff_matrix_dim1 + 3 + 1)
    error->all(FLERR,
               "UF3: {} has incorrect knot (NUM_OF_KNOTS_IJ) and "
               "coeff ()coeff_matrix_dim1 data nknots!=ncoeffs + 3 +1",
               potf_name);

  coeff_matrix_elements_len = coeff_matrix_dim3;

  std::string key = std::to_string(itype) + std::to_string(jtype) + std::to_string(ktype);
  uf3_impl->n3b_coeff_matrix[key].resize(coeff_matrix_dim1);

  int line_count = 0;
  for (int i = 0; i < coeff_matrix_dim1; i++) {
    uf3_impl->n3b_coeff_matrix[key][i].resize(coeff_matrix_dim2);
    for (int j = 0; j < coeff_matrix_dim2; j++) {
      temp_line = txtfilereader.next_line(coeff_matrix_elements_len);
      ValueTokenizer coeff_line(temp_line);
      uf3_impl->n3b_coeff_matrix[key][i][j].resize(coeff_matrix_dim3);

      if (coeff_line.count() != coeff_matrix_elements_len)
        error->all(FLERR, "UF3: Expected {} numbers on {}th line but found {} numbers",
                   coeff_matrix_elements_len, line_count + 8, coeff_line.count());
      for (int k = 0; k < coeff_matrix_dim3; k++) {
        uf3_impl->n3b_coeff_matrix[key][i][j][k] = coeff_line.next_double();
      }
      line_count += 1;
    }
  }

  std::string key2 = std::to_string(itype) + std::to_string(ktype) + std::to_string(jtype);
  uf3_impl->n3b_coeff_matrix[key2].resize(coeff_matrix_dim2);
  for (int j = 0; j < coeff_matrix_dim2; j++) {
    uf3_impl->n3b_coeff_matrix[key2][j].resize(coeff_matrix_dim1);
    for (int i = 0; i < coeff_matrix_dim1; i++) {
      uf3_impl->n3b_coeff_matrix[key2][j][i].resize(coeff_matrix_dim3);
    }
  }

  for (int i = 0; i < coeff_matrix_dim1; i++) {
    for (int j = 0; j < coeff_matrix_dim2; j++) {
      for (int k = 0; k < coeff_matrix_dim3; k++) {
        uf3_impl->n3b_coeff_matrix[key2][j][i][k] = uf3_impl->n3b_coeff_matrix[key][i][j][k];
      }
    }
  }

  setflag_3b[itype][jtype][ktype] = 1;
  setflag_3b[itype][ktype][jtype] = 1;
  fclose(fp);
}

void PairUF3::uf3_read_pot_file(char *potf_name)
{
  FILE *fp = utils::open_potential(potf_name, lmp, nullptr);
  if (!fp)
    error->all(FLERR, "Cannot open UF3 potential file {}: {}", potf_name, utils::getsyserror());

  TextFileReader txtfilereader(fp, "UF3:POTFP");
  txtfilereader.ignore_comments = false;

  std::string temp_line = txtfilereader.next_line(2);
  Tokenizer fp1st_line(temp_line);

  if (fp1st_line.count() <= 2)
    error->all(FLERR,
               "UF3: Expected more than two words on 1st line of {} \n"
               "but found {} word/s",
               potf_name, fp1st_line.count());

  if (fp1st_line.contains("#UF3 POT") == 0)
    error->all(FLERR,
               "UF3: {} file is not UF3 POT type, 1st line of UF3 POT \n"
               "files contain '#UF3 POT'. Found {} in the header",
               potf_name, temp_line);

  if (fp1st_line.contains("UNITS:") == 0)
    error->all(FLERR,
               "UF3: {} file does not contain the 'UNITS:' metadata in \n"
               "the header",
               potf_name);

  temp_line = txtfilereader.next_line(1);
  Tokenizer fp2nd_line(temp_line);
  if (fp2nd_line.contains("2B") == 1) {
    temp_line = txtfilereader.next_line(4);
    ValueTokenizer fp3rd_line(temp_line);
    int temp_type1 = fp3rd_line.next_int();
    int temp_type2 = fp3rd_line.next_int();

    //cut is used in init_one which is called by pair.cpp at line 267 where the return of init_one is squared
    cut[temp_type1][temp_type2] = fp3rd_line.next_double();
    cut[temp_type2][temp_type1] = cut[temp_type1][temp_type2];

    int temp_line_len = fp3rd_line.next_int();

    temp_line = txtfilereader.next_line(temp_line_len);
    ValueTokenizer fp4th_line(temp_line);

    uf3_impl->n2b_knot[temp_type1][temp_type2].resize(temp_line_len);
    uf3_impl->n2b_knot[temp_type2][temp_type1].resize(temp_line_len);
    for (int k = 0; k < temp_line_len; k++) {
      uf3_impl->n2b_knot[temp_type1][temp_type2][k] = fp4th_line.next_double();
      uf3_impl->n2b_knot[temp_type2][temp_type1][k] = uf3_impl->n2b_knot[temp_type1][temp_type2][k];
    }

    temp_line = txtfilereader.next_line(1);
    ValueTokenizer fp5th_line(temp_line);

    temp_line_len = fp5th_line.next_int();

    temp_line = txtfilereader.next_line(temp_line_len);
    ValueTokenizer fp6th_line(temp_line);
    uf3_impl->n2b_coeff[temp_type1][temp_type2].resize(temp_line_len);
    uf3_impl->n2b_coeff[temp_type2][temp_type1].resize(temp_line_len);

    for (int k = 0; k < temp_line_len; k++) {
      uf3_impl->n2b_coeff[temp_type1][temp_type2][k] = fp6th_line.next_double();
      uf3_impl->n2b_coeff[temp_type2][temp_type1][k] = uf3_impl->n2b_coeff[temp_type1][temp_type2][k];
    }
    if (uf3_impl->n2b_knot[temp_type1][temp_type2].size() != uf3_impl->n2b_coeff[temp_type1][temp_type2].size() + 4) {
      error->all(FLERR, "UF3: {} has incorrect knot and coeff data nknots!=ncoeffs + 3 +1",
                 potf_name);
    }
    setflag[temp_type1][temp_type2] = 1;
    setflag[temp_type2][temp_type1] = 1;
  } else if (fp2nd_line.contains("3B") == 1) {
    int coeff_matrix_dim1, coeff_matrix_dim2, coeff_matrix_dim3, coeff_matrix_elements_len;
    temp_line = txtfilereader.next_line(9);
    ValueTokenizer fp3rd_line(temp_line);
    int temp_type1 = fp3rd_line.next_int();
    int temp_type2 = fp3rd_line.next_int();
    int temp_type3 = fp3rd_line.next_int();

    double cut3b_rjk = fp3rd_line.next_double();
    double cut3b_rij = fp3rd_line.next_double();
    // cut_3b[temp_type1][temp_type2] = std::max(cut3b_rij,
    // cut_3b[temp_type1][temp_type2]);
    cut_3b_list[temp_type1][temp_type2] =
        std::max(cut3b_rij, cut_3b_list[temp_type1][temp_type2]);

    double cut3b_rik = fp3rd_line.next_double();
    if (cut3b_rij != cut3b_rik) {
      error->all(FLERR, "UF3: rij!=rik, Current implementation only works for rij=rik");
    }
    if (2 * cut3b_rik != cut3b_rjk) {
      error->all(FLERR,
                 "UF3: 2rij=2rik!=rik, Current implementation only works for 2rij=2rik!=rik");
    }
    // cut_3b[temp_type1][temp_type3] = std::max(cut_3b[temp_type1][temp_type3],cut3b_rik);
    cut_3b_list[temp_type1][temp_type3] =
        std::max(cut_3b_list[temp_type1][temp_type3], cut3b_rik);

    cut_3b[temp_type1][temp_type2][temp_type3] = cut3b_rij;
    cut_3b[temp_type1][temp_type3][temp_type2] = cut3b_rik;

    int temp_line_len = fp3rd_line.next_int();
    temp_line = txtfilereader.next_line(temp_line_len);
    ValueTokenizer fp4th_line(temp_line);

    uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3].resize(3);
    uf3_impl->n3b_knot_matrix[temp_type1][temp_type3][temp_type2].resize(3);

    uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][0].resize(temp_line_len);
    uf3_impl->n3b_knot_matrix[temp_type1][temp_type3][temp_type2][0].resize(temp_line_len);
    for (int i = 0; i < temp_line_len; i++) {
      uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][0][i] = fp4th_line.next_double();
      uf3_impl->n3b_knot_matrix[temp_type1][temp_type3][temp_type2][0][i] =
          uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][0][i];
    }

    min_cut_3b[temp_type1][temp_type2][temp_type3][0] =
        uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][0][0];
    min_cut_3b[temp_type1][temp_type3][temp_type2][0] =
        uf3_impl->n3b_knot_matrix[temp_type1][temp_type3][temp_type2][0][0];

    temp_line_len = fp3rd_line.next_int();
    temp_line = txtfilereader.next_line(temp_line_len);
    ValueTokenizer fp5th_line(temp_line);
    uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][1].resize(temp_line_len);
    uf3_impl->n3b_knot_matrix[temp_type1][temp_type3][temp_type2][2].resize(temp_line_len);
    for (int i = 0; i < temp_line_len; i++) {
      uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][1][i] = fp5th_line.next_double();
      uf3_impl->n3b_knot_matrix[temp_type1][temp_type3][temp_type2][2][i] =
          uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][1][i];
    }

    min_cut_3b[temp_type1][temp_type2][temp_type3][1] =
        uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][1][0];
    min_cut_3b[temp_type1][temp_type3][temp_type2][2] =
        uf3_impl->n3b_knot_matrix[temp_type1][temp_type3][temp_type2][2][0];

    temp_line_len = fp3rd_line.next_int();
    temp_line = txtfilereader.next_line(temp_line_len);
    ValueTokenizer fp6th_line(temp_line);
    uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][2].resize(temp_line_len);
    uf3_impl->n3b_knot_matrix[temp_type1][temp_type3][temp_type2][1].resize(temp_line_len);
    for (int i = 0; i < temp_line_len; i++) {
      uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][2][i] = fp6th_line.next_double();
      uf3_impl->n3b_knot_matrix[temp_type1][temp_type3][temp_type2][1][i] =
          uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][2][i];
    }

    min_cut_3b[temp_type1][temp_type2][temp_type3][2] =
        uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][2][0];
    min_cut_3b[temp_type1][temp_type3][temp_type2][1] =
        uf3_impl->n3b_knot_matrix[temp_type1][temp_type3][temp_type2][1][0];

    temp_line = txtfilereader.next_line(3);
    ValueTokenizer fp7th_line(temp_line);

    coeff_matrix_dim1 = fp7th_line.next_int();
    coeff_matrix_dim2 = fp7th_line.next_int();
    coeff_matrix_dim3 = fp7th_line.next_int();
    if (uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][0].size() !=
        coeff_matrix_dim3 + 3 + 1) {
      error->all(FLERR, "UF3: {} has incorrect knot and coeff data nknots!=ncoeffs + 3 +1",
                 potf_name);
    }
    if (uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][1].size() !=
        coeff_matrix_dim2 + 3 + 1) {
      error->all(FLERR, "UF3: {} has incorrect knot and coeff data nknots!=ncoeffs + 3 +1",
                 potf_name);
    }
    if (uf3_impl->n3b_knot_matrix[temp_type1][temp_type2][temp_type3][2].size() !=
        coeff_matrix_dim1 + 3 + 1) {
      error->all(FLERR, "UF3: {} has incorrect knot and coeff data nknots!=ncoeffs + 3 +1",
                 potf_name);
    }

    coeff_matrix_elements_len = coeff_matrix_dim3;

    std::string key =
        std::to_string(temp_type1) + std::to_string(temp_type2) + std::to_string(temp_type3);
    uf3_impl->n3b_coeff_matrix[key].resize(coeff_matrix_dim1);
    for (int i = 0; i < coeff_matrix_dim1; i++) {
      uf3_impl->n3b_coeff_matrix[key][i].resize(coeff_matrix_dim2);
      for (int j = 0; j < coeff_matrix_dim2; j++) {
        temp_line = txtfilereader.next_line(coeff_matrix_elements_len);
        ValueTokenizer coeff_line(temp_line);
        uf3_impl->n3b_coeff_matrix[key][i][j].resize(coeff_matrix_dim3);
        for (int k = 0; k < coeff_matrix_dim3; k++) {
          uf3_impl->n3b_coeff_matrix[key][i][j][k] = coeff_line.next_double();
        }
      }
    }

    key = std::to_string(temp_type1) + std::to_string(temp_type3) + std::to_string(temp_type2);
    uf3_impl->n3b_coeff_matrix[key] =
        uf3_impl->n3b_coeff_matrix[std::to_string(temp_type1) +
                std::to_string(temp_type2) +
                std::to_string(temp_type3)];
    setflag_3b[temp_type1][temp_type2][temp_type3] = 1;
    setflag_3b[temp_type1][temp_type3][temp_type2] = 1;
  } else
    error->all(
        FLERR,
        "UF3: {} file does not contain right words indicating whether it is 2 or 3 body potential",
        potf_name);
  fclose(fp);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
void PairUF3::init_style()
{
  if (force->newton_pair == 0) error->all(FLERR, "UF3: Pair style requires newton pair on");
  // request a default neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init list sets the pointer to full neighbour list requested in previous function
------------------------------------------------------------------------- */

void PairUF3::init_list(int /*id*/, class NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
double PairUF3::init_one(int i /*i*/, int /*j*/ j)
{

  if (!bsplines_created) create_bsplines();

  //init_one is called by pair.cpp at line 267 where it is squred
  //at line 268
  return cut[i][j];
}

void PairUF3::create_bsplines()
{
  const int num_of_elements = atom->ntypes;
  bsplines_created = 1;
  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      if (setflag[i][j] != 1)
        error->all(FLERR,
                   "UF3: Not all 2-body UF potentials are set, "
                   "missing potential for {}-{} interaction",
                   i, j);
    }
  }
  if (pot_3b) {
    for (int i = 1; i < num_of_elements + 1; i++) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        for (int k = 1; k < num_of_elements + 1; k++) {
          if (setflag_3b[i][j][k] != 1)
            error->all(FLERR,
                       "UF3: Not all 3-body UF potentials are set, "
                       "missing potential for {}-{}-{} interaction",
                       i, j, k);
        }
      }
    }
  }

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = i; j < num_of_elements + 1; j++) {
      /*uf3_impl->UFBS2b[i][j] =
          uf3_pair_bspline(lmp, uf3_impl->n2b_knot[i][j], uf3_impl->n2b_coeff[i][j],
                  knot_spacing_type_2b[i][j]);*/

        uf3_impl->UFBS2b[i][j] = uf3_pair_bspline(lmp, n2b_knots_array[i][j],
                                                  n2b_knots_array_size[i][j],
                                                  n2b_coeff_array[i][j],
                                                  n2b_coeff_array_size[i][j],
                                                  knot_spacing_type_2b[i][j]);
      uf3_impl->UFBS2b[j][i] = uf3_impl->UFBS2b[i][j];
    }
    if (pot_3b) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        for (int k = j; k < num_of_elements + 1; k++) {
          /*std::string key = std::to_string(i) + std::to_string(j) + std::to_string(k);
          uf3_impl->UFBS3b[i][j][k] = uf3_triplet_bspline(
              lmp, uf3_impl->n3b_knot_matrix[i][j][k], uf3_impl->n3b_coeff_matrix[key], knot_spacing_type_3b[i][j][k]);*/
          int key = map_3b[i][j][k];
          int key2 = map_3b[i][k][j];
          /*utils::logmesg(lmp, "Setting UFBS3b for {}-{}-{} map_3b={} and for {}-{}-{} "
                         "map_3b={}\n", i, j, k, key, i, k, j,
                         key2);*/
          uf3_impl->UFBS3b[i][j][k] = uf3_triplet_bspline(
              lmp, n3b_knots_array[key], n3b_knots_array_size[key],
              n3b_coeff_array[key], n3b_coeff_array_size[key],
              knot_spacing_type_3b[i][j][k]);

          /*std::string key2 = std::to_string(i) + std::to_string(k) + std::to_string(j);
          uf3_impl->UFBS3b[i][k][j] = uf3_triplet_bspline(
              lmp, uf3_impl->n3b_knot_matrix[i][k][j], uf3_impl->n3b_coeff_matrix[key2], knot_spacing_type_3b[i][k][j]);*/
          //int key2 = map_3b[i][k][j];
          uf3_impl->UFBS3b[i][k][j] = uf3_triplet_bspline(
              lmp, n3b_knots_array[key2], n3b_knots_array_size[key2],
              n3b_coeff_array[key2], n3b_coeff_array_size[key2],
              knot_spacing_type_3b[i][k][j]);
        }
      }
    }
  }
}

void PairUF3::compute(int eflag, int vflag)
{
  int i, j, k, ii, jj, kk, inum, jnum, itype, jtype, ktype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair, fx, fy, fz;
  double del_rji[3], del_rki[3], del_rkj[3];
  double fij[3], fik[3], fjk[3];
  double fji[3], fki[3], fkj[3];
  double Fi[3], Fj[3], Fk[3];
  double rsq, rij, rik, rjk;
  int *ilist, *jlist, *numneigh, **firstneigh;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    evdwl = 0;
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;
    for (jj = 0; jj < jnum; jj++) {
      fx = 0;
      fy = 0;
      fz = 0;
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      if (rsq < cutsq[itype][jtype]) {
        rij = sqrt(rsq);

        if (pot_3b) {
          if (rij <= cut_3b_list[itype][jtype]) {
            neighshort[numshort] = j;
            if (numshort >= maxshort - 1) {
              maxshort += maxshort / 2;
              memory->grow(neighshort, maxshort, "pair:neighshort");
            }
            numshort = numshort + 1;
          }
        }

        double *pair_eval = uf3_impl->UFBS2b[itype][jtype].eval(rij);

        fpair = -1 * pair_eval[1] / rij;

        fx = delx * fpair;
        fy = dely * fpair;
        fz = delz * fpair;

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;

        if (eflag) evdwl = pair_eval[0];

        if (evflag) {
          ev_tally_xyz(i, j, nlocal, newton_pair, evdwl, 0.0, fx, fy, fz, delx, dely, delz);

          // Centroid Stress
          if (vflag_either && cvflag_atom) {
            double v[6];

            v[0] = delx * fx;
            v[1] = dely * fy;
            v[2] = delz * fz;
            v[3] = delx * fy;
            v[4] = delx * fz;
            v[5] = dely * fz;

            cvatom[i][0] += 0.5 * v[0];
            cvatom[i][1] += 0.5 * v[1];
            cvatom[i][2] += 0.5 * v[2];
            cvatom[i][3] += 0.5 * v[3];
            cvatom[i][4] += 0.5 * v[4];
            cvatom[i][5] += 0.5 * v[5];
            cvatom[i][6] += 0.5 * v[3];
            cvatom[i][7] += 0.5 * v[4];
            cvatom[i][8] += 0.5 * v[5];

            cvatom[j][0] += 0.5 * v[0];
            cvatom[j][1] += 0.5 * v[1];
            cvatom[j][2] += 0.5 * v[2];
            cvatom[j][3] += 0.5 * v[3];
            cvatom[j][4] += 0.5 * v[4];
            cvatom[j][5] += 0.5 * v[5];
            cvatom[j][6] += 0.5 * v[3];
            cvatom[j][7] += 0.5 * v[4];
            cvatom[j][8] += 0.5 * v[5];
          }
        }
      }
    }

    // 3-body interaction
    // jth atom
    jnum = numshort - 1;
    for (jj = 0; jj < jnum; jj++) {
      fij[0] = fji[0] = 0;
      fij[1] = fji[1] = 0;
      fij[2] = fji[2] = 0;
      j = neighshort[jj];
      jtype = type[j];
      del_rji[0] = x[j][0] - xtmp;
      del_rji[1] = x[j][1] - ytmp;
      del_rji[2] = x[j][2] - ztmp;
      rij =
          sqrt(((del_rji[0] * del_rji[0]) + (del_rji[1] * del_rji[1]) + (del_rji[2] * del_rji[2])));

      // kth atom
      for (kk = jj + 1; kk < numshort; kk++) {

        fik[0] = fki[0] = 0;
        fik[1] = fki[1] = 0;
        fik[2] = fki[2] = 0;

        fjk[0] = fkj[0] = 0;
        fjk[1] = fkj[1] = 0;
        fjk[2] = fkj[2] = 0;

        k = neighshort[kk];
        ktype = type[k];
        del_rki[0] = x[k][0] - xtmp;
        del_rki[1] = x[k][1] - ytmp;
        del_rki[2] = x[k][2] - ztmp;
        rik = sqrt(
            ((del_rki[0] * del_rki[0]) + (del_rki[1] * del_rki[1]) + (del_rki[2] * del_rki[2])));

        if ((rij <= cut_3b[itype][jtype][ktype]) &&
                (rik <= cut_3b[itype][ktype][jtype]) &&
                (rij >= min_cut_3b[itype][jtype][ktype][2]) &&
                (rik >= min_cut_3b[itype][jtype][ktype][1])) {

          del_rkj[0] = x[k][0] - x[j][0];
          del_rkj[1] = x[k][1] - x[j][1];
          del_rkj[2] = x[k][2] - x[j][2];
          rjk = sqrt(
              ((del_rkj[0] * del_rkj[0]) + (del_rkj[1] * del_rkj[1]) + (del_rkj[2] * del_rkj[2])));

          if (rjk >= min_cut_3b[itype][jtype][ktype][0]) {
            double *triangle_eval = uf3_impl->UFBS3b[itype][jtype][ktype].eval(rij, rik, rjk);

            fij[0] = *(triangle_eval + 1) * (del_rji[0] / rij);
            fji[0] = -fij[0];
            fik[0] = *(triangle_eval + 2) * (del_rki[0] / rik);
            fki[0] = -fik[0];
            fjk[0] = *(triangle_eval + 3) * (del_rkj[0] / rjk);
            fkj[0] = -fjk[0];

            fij[1] = *(triangle_eval + 1) * (del_rji[1] / rij);
            fji[1] = -fij[1];
            fik[1] = *(triangle_eval + 2) * (del_rki[1] / rik);
            fki[1] = -fik[1];
            fjk[1] = *(triangle_eval + 3) * (del_rkj[1] / rjk);
            fkj[1] = -fjk[1];

            fij[2] = *(triangle_eval + 1) * (del_rji[2] / rij);
            fji[2] = -fij[2];
            fik[2] = *(triangle_eval + 2) * (del_rki[2] / rik);
            fki[2] = -fik[2];
            fjk[2] = *(triangle_eval + 3) * (del_rkj[2] / rjk);
            fkj[2] = -fjk[2];

            Fi[0] = fij[0] + fik[0];
            Fi[1] = fij[1] + fik[1];
            Fi[2] = fij[2] + fik[2];
            f[i][0] += Fi[0];
            f[i][1] += Fi[1];
            f[i][2] += Fi[2];

            Fj[0] = fji[0] + fjk[0];
            Fj[1] = fji[1] + fjk[1];
            Fj[2] = fji[2] + fjk[2];
            f[j][0] += Fj[0];
            f[j][1] += Fj[1];
            f[j][2] += Fj[2];

            Fk[0] = fki[0] + fkj[0];
            Fk[1] = fki[1] + fkj[1];
            Fk[2] = fki[2] + fkj[2];
            f[k][0] += Fk[0];
            f[k][1] += Fk[1];
            f[k][2] += Fk[2];

            if (eflag) evdwl = *triangle_eval;

            if (evflag) {
              ev_tally3(i, j, k, evdwl, 0, Fj, Fk, del_rji, del_rki);
              // Centroid stress 3-body term
              if (vflag_either && cvflag_atom) {
                double ric[3];
                ric[0] = THIRD * (-del_rji[0] - del_rki[0]);
                ric[1] = THIRD * (-del_rji[1] - del_rki[1]);
                ric[2] = THIRD * (-del_rji[2] - del_rki[2]);

                cvatom[i][0] += ric[0] * Fi[0];
                cvatom[i][1] += ric[1] * Fi[1];
                cvatom[i][2] += ric[2] * Fi[2];
                cvatom[i][3] += ric[0] * Fi[1];
                cvatom[i][4] += ric[0] * Fi[2];
                cvatom[i][5] += ric[1] * Fi[2];
                cvatom[i][6] += ric[1] * Fi[0];
                cvatom[i][7] += ric[2] * Fi[0];
                cvatom[i][8] += ric[2] * Fi[1];

                double rjc[3];
                rjc[0] = THIRD * (del_rji[0] - del_rkj[0]);
                rjc[1] = THIRD * (del_rji[1] - del_rkj[1]);
                rjc[2] = THIRD * (del_rji[2] - del_rkj[2]);

                cvatom[j][0] += rjc[0] * Fj[0];
                cvatom[j][1] += rjc[1] * Fj[1];
                cvatom[j][2] += rjc[2] * Fj[2];
                cvatom[j][3] += rjc[0] * Fj[1];
                cvatom[j][4] += rjc[0] * Fj[2];
                cvatom[j][5] += rjc[1] * Fj[2];
                cvatom[j][6] += rjc[1] * Fj[0];
                cvatom[j][7] += rjc[2] * Fj[0];
                cvatom[j][8] += rjc[2] * Fj[1];

                double rkc[3];
                rkc[0] = THIRD * (del_rki[0] + del_rkj[0]);
                rkc[1] = THIRD * (del_rki[1] + del_rkj[1]);
                rkc[2] = THIRD * (del_rki[2] + del_rkj[2]);

                cvatom[k][0] += rkc[0] * Fk[0];
                cvatom[k][1] += rkc[1] * Fk[1];
                cvatom[k][2] += rkc[2] * Fk[2];
                cvatom[k][3] += rkc[0] * Fk[1];
                cvatom[k][4] += rkc[0] * Fk[2];
                cvatom[k][5] += rkc[1] * Fk[2];
                cvatom[k][6] += rkc[1] * Fk[0];
                cvatom[k][7] += rkc[2] * Fk[0];
                cvatom[k][8] += rkc[2] * Fk[1];
              }
            }
          }
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

double PairUF3::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                       double /*factor_coul*/, double factor_lj, double &fforce)
{
  double value = 0.0;
  double r = sqrt(rsq);

  if (r < cut[itype][jtype]) {
    double *pair_eval = uf3_impl->UFBS2b[itype][jtype].eval(r);
    value = pair_eval[0];
    fforce = factor_lj * pair_eval[1];
  }

  return factor_lj * value;
}

double PairUF3::memory_usage()
{
  const int num_of_elements = atom->ntypes;
  double bytes = Pair::memory_usage();

  bytes += (double) 5 * sizeof(double);    //num_of_elements, nbody_flag,
                                           //n2body_pot_files, n3body_pot_files,
                                           //tot_pot_files;

  bytes += (double) 5 * sizeof(double);    //bsplines_created, coeff_matrix_dim1,
                                           //coeff_matrix_dim2, coeff_matrix_dim3,
                                           //coeff_matrix_elements_len
  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) *
      sizeof(double);    //***setflag_3b

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * sizeof(double);    //cut

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) *
      sizeof(double);    //***cut_3b

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * sizeof(double);    //cut_3b_list

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) * 3 *
      sizeof(double);    //min_cut_3b

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * sizeof(double);   //n2b_knots_array_size
  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * sizeof(double);   //n2b_coeff_array_size
  
  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * max_num_knots_2b *
                    sizeof(double);         //n2b_knots_array

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * max_num_coeff_2b *
                    sizeof(double);         //n2b_coeff_array

  if (pot_3b) {
    bytes += (double) tot_interaction_count_3b * 3 * sizeof(double);       //n3b_knots_array_size
    bytes += (double) tot_interaction_count_3b * 3 * sizeof(double);       //n3b_coeff_array_size
    bytes += (double) tot_interaction_count_3b * 3 * max_num_knots_3b * sizeof(double); //n3b_knots_array
    bytes += (double) tot_interaction_count_3b * max_num_coeff_3b * max_num_coeff_3b *
                      max_num_coeff_3b * sizeof(double); //n3b_coeff_array
  }
  /*for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = i; j < num_of_elements + 1; j++) {
      bytes += (double) 2 * uf3_impl->n2b_knot[i][j].size() * sizeof(double);     //n2b_knot
      bytes += (double) 2 * uf3_impl->n2b_coeff[i][j].size() * sizeof(double);    //n2b_coeff
    }
    if (pot_3b) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        for (int k = j; k < num_of_elements + 1; k++) {
          bytes += (double) 2 * uf3_impl->n3b_knot_matrix[i][j][k][0].size() * sizeof(double);
          bytes += (double) 2 * uf3_impl->n3b_knot_matrix[i][j][k][1].size() * sizeof(double);
          bytes += (double) 2 * uf3_impl->n3b_knot_matrix[i][j][k][2].size() * sizeof(double);

          std::string key = std::to_string(i) + std::to_string(j) + std::to_string(k);

          for (int l = 0; l < uf3_impl->n3b_coeff_matrix[key].size(); l++) {
            for (int m = 0; m < uf3_impl->n3b_coeff_matrix[key][l].size(); m++) {
              bytes += (double) 2 * uf3_impl->n3b_coeff_matrix[key][l][m].size() * sizeof(double);
              //key = ijk
              //key = ikj
            }
          }
        }
      }
    }
  }*/

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = i; j < num_of_elements + 1; j++) {
      bytes += (double) 2 * uf3_impl->UFBS2b[i][j].memory_usage();    //UFBS2b[i][j] UFBS2b[j][1]
    }
    if (pot_3b) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        for (int k = j; k < num_of_elements + 1; k++) {
          bytes += (double) 2 * uf3_impl->UFBS3b[i][j][k].memory_usage();
        }
      }
    }
  }

  bytes += (double) (maxshort + 1) * sizeof(int);    //neighshort, maxshort

  return bytes;
}

//Accessor function called by pair_uf3_kokkos.cpp
//Will probably be removed once std::vector are converted to arrays
std::vector<std::vector<std::vector<double>>>& PairUF3::get_n2b_knot()
{
  return uf3_impl->n2b_knot;
}

std::vector<std::vector<std::vector<double>>>& PairUF3::get_n2b_coeff()
{
  return uf3_impl->n2b_coeff;
}
//Accessor function called by pair_uf3_kokkos.cpp
//Will probably be removed once std::vector are converted to arrays
std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>&
    PairUF3::get_n3b_knot_matrix()
{
  return uf3_impl->n3b_knot_matrix;
}

//Accessor function called by pair_uf3_kokkos.cpp
//Will probably be removed once std::vector are converted to arrays
std::vector<std::vector<std::vector<double>>>&
    PairUF3::get_n3b_coeff_matrix_key(std::string key)
{
  return uf3_impl->n3b_coeff_matrix[key];
}

double PairUF3::get_knot_spacing_2b(int i, int j)
{
  return uf3_impl->UFBS2b[i][j].knot_spacing;
}

double PairUF3::get_knot_spacing_3b_ij(int i, int j, int k)
{
  return uf3_impl->UFBS3b[i][j][k].knot_spacing_ij;
}
double PairUF3::get_knot_spacing_3b_ik(int i, int j, int k)
{
  return uf3_impl->UFBS3b[i][j][k].knot_spacing_ik;
}
double PairUF3::get_knot_spacing_3b_jk(int i, int j, int k)
{
  return uf3_impl->UFBS3b[i][j][k].knot_spacing_jk;
}

