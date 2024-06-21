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

#include "uf3_bspline_basis2.h"
#include "uf3_bspline_basis3.h"

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

using namespace LAMMPS_NS;
using MathConst::THIRD;

/* ---------------------------------------------------------------------- */

PairUF3::PairUF3(LAMMPS *lmp) :
    Pair(lmp), setflag_3b(nullptr), knot_spacing_type_2b(nullptr), knot_spacing_type_3b(nullptr),
    cut(nullptr), cut_3b(nullptr), cut_3b_list(nullptr), min_cut_3b(nullptr),
    knot_spacing_2b(nullptr), knot_spacing_3b(nullptr), n2b_knots_array(nullptr),
    n2b_coeff_array(nullptr), n2b_knots_array_size(nullptr), n2b_coeff_array_size(nullptr),
    cached_constants_2b(nullptr), cached_constants_2b_deri(nullptr), map_3b(nullptr),
    n3b_knots_array(nullptr), n3b_coeff_array(nullptr), n3b_knots_array_size(nullptr),
    n3b_coeff_array_size(nullptr), coeff_for_der_jk(nullptr), coeff_for_der_ik(nullptr),
    coeff_for_der_ij(nullptr), cached_constants_3b(nullptr), cached_constants_3b_deri(nullptr),
    neighshort(nullptr), get_starting_index_2b(nullptr), get_starting_index_3b(nullptr)
{
  single_enable = 1;    // 1 if single() routine exists
  one_coeff = 1;        // 1 if allows only one coeff * * call
  restartinfo = 0;      // 1 if pair style writes restart info
  maxshort = 20;
  centroidstressflag = CENTROID_AVAIL;
  manybody_flag = 1;
  bsplines_created = 0;
  pot_3b = false;
  nbody_flag = 3;
  max_num_knots_2b = 0;
  max_num_coeff_2b = 0;
  max_num_knots_3b = 0;
  max_num_coeff_3b = 0;
  tot_interaction_count_3b = 0;
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
    memory->destroy(knot_spacing_2b);
    memory->destroy(n2b_knots_array_size);
    memory->destroy(n2b_coeff_array_size);
    memory->destroy(n2b_knots_array);
    memory->destroy(n2b_coeff_array);
    memory->destroy(cached_constants_2b);
    memory->destroy(cached_constants_2b_deri);

    if (pot_3b) {
      memory->destroy(setflag_3b);
      memory->destroy(cut_3b);
      memory->destroy(cut_3b_list);
      memory->destroy(min_cut_3b);
      memory->destroy(neighshort);
      memory->destroy(knot_spacing_type_3b);
      memory->destroy(knot_spacing_3b);
      memory->destroy(map_3b);
      memory->destroy(n3b_knots_array_size);
      memory->destroy(n3b_coeff_array_size);
      memory->destroy(n3b_knots_array);
      memory->destroy(n3b_coeff_array);
      memory->destroy(coeff_for_der_jk);
      memory->destroy(coeff_for_der_ik);
      memory->destroy(coeff_for_der_ij);
      memory->destroy(cached_constants_3b);
      memory->destroy(cached_constants_3b_deri);
    }
  }
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
  } else if (nbody_flag == 3) {
    pot_3b = true;
    single_enable = 0;
  } else
    error->all(FLERR, "Pair style uf3 not (yet) implemented for {}-body terms", nbody_flag);
}

/* ----------------------------------------------------------------------
 *    set coeffs for one or more type pairs
 * ---------------------------------------------------------------------- */
void PairUF3::coeff(int narg, char **arg)
{
  if (narg != 3 + atom->ntypes)
    error->all(FLERR, "Invalid number of arguments uf3 in pair coeffs.");

  if (!allocated) allocate();

  map_element2type(narg - 3, arg + 3, false);

  if (comm->me == 0) uf3_read_unified_pot_file(arg[2]);
  communicate();
}

void PairUF3::allocate()
{
  allocated = 1;
  const int num_of_elements = atom->ntypes;

  map = new int[num_of_elements + 1];    //No need to delete map as ~Pair deletes map

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
                 "pair:knot_spacing_type_2b");
  memory->create(knot_spacing_2b, num_of_elements + 1, num_of_elements + 1, "pair:knot_spacing_2b");

  //Contains size of 2b knots vectors and 2b coeff matrices
  memory->create(n2b_knots_array_size, num_of_elements + 1, num_of_elements + 1,
                 "pair:n2b_knots_array_size");
  memory->create(n2b_coeff_array_size, num_of_elements + 1, num_of_elements + 1,
                 "pair:n2b_coeff_array_size");

  if (pot_3b) {
    // Contains info about wether UF potential were found for type i, j and k
    memory->create(setflag_3b, num_of_elements + 1, num_of_elements + 1, num_of_elements + 1,
                   "pair:setflag_3b");
    // Contains info about 3-body cutoff distance for type i, j and k
    memory->create(cut_3b, num_of_elements + 1, num_of_elements + 1, num_of_elements + 1,
                   "pair:cut_3b");
    // Contains info about 3-body cutoff distance for type i, j and k
    // for constructing 3-body list
    memory->create(cut_3b_list, num_of_elements + 1, num_of_elements + 1, "pair:cut_3b_list");
    // Contains info about minimum 3-body cutoff distance for type i, j and k
    memory->create(min_cut_3b, num_of_elements + 1, num_of_elements + 1, num_of_elements + 1, 3,
                   "pair:min_cut_3b");
    //Contains info about type of knot_spacing--> 0 = uniform knot spacing (default)
    //1 = non-uniform knot spacing
    memory->create(knot_spacing_type_3b, num_of_elements + 1, num_of_elements + 1,
                   num_of_elements + 1, "pair:knot_spacing_type_3b");
    memory->create(knot_spacing_3b, num_of_elements + 1, num_of_elements + 1, num_of_elements + 1,
                   3, "pair:knot_spacing_3b");

    tot_interaction_count_3b = 0;
    //conatins map of I-J-K interaction
    memory->create(map_3b, num_of_elements + 1, num_of_elements + 1, num_of_elements + 1,
                   "pair:map_3b");

    // setting cut_3b, setflag = 0 and map_3b
    for (int i = 1; i < num_of_elements + 1; i++) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        cut_3b_list[i][j] = 0;
        setflag[i][j] = 0;
        n2b_coeff_array_size[i][j] = 0;
        n2b_knots_array_size[i][j] = 0;

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
    memory->create(n3b_knots_array_size, tot_interaction_count_3b, 3, "pair:n3b_knots_array_size");
    memory->create(n3b_coeff_array_size, tot_interaction_count_3b, 3, "pair:n3b_coeff_array_size");
    for (int i = 0; i < tot_interaction_count_3b; i++) {
      n3b_coeff_array_size[i][0] = 0;
      n3b_coeff_array_size[i][1] = 0;
      n3b_coeff_array_size[i][2] = 0;

      n3b_knots_array_size[i][0] = 0;
      n3b_knots_array_size[i][1] = 0;
      n3b_knots_array_size[i][2] = 0;
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
    error->all(FLERR, "Cannot open UF3 potential file {}: {}", potf_name, utils::getsyserror());

  TextFileReader txtfilereader(fp, "UF3:POTFP");
  txtfilereader.ignore_comments = false;

  //while loop over the entire file, find blocks starting with #UF3 POT
  //if block found read the very next line to determine 2B or 3B block
  //if 2B read the knot vector and coeff vector size
  //if 3B read the knot vectors and coeff matrix size
  int line_counter = 1;
  char *line;
  while ((line = txtfilereader.next_line(1))) {
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
          error->all(FLERR,
                     "UF3: Expected 6 words on line {} of {} file but found {} word/s",
                     line_counter, potf_name, fp2nd_line.count());

        //get the elements
        std::string element1 = fp2nd_line.next_string();
        std::string element2 = fp2nd_line.next_string();
        int itype = 0;
        int jtype = 0;
        for (int i = 1; i < num_of_elements + 1; i++) {
          if (std::string(elements[map[i]]) == element1) {
            itype = i;
            break;
          }
        }
        for (int i = 1; i < num_of_elements + 1; i++) {
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
          cut[itype][jtype] = fp3rd_line.next_double();
          cut[jtype][itype] = cut[itype][jtype];

          int num_knots_2b = fp3rd_line.next_int();
          n2b_knots_array_size[itype][jtype] = num_knots_2b;
          n2b_knots_array_size[jtype][itype] = num_knots_2b;
          max_num_knots_2b = std::max(max_num_knots_2b, num_knots_2b);

          //skip next line
          txtfilereader.skip_line();
          line_counter++;

          //read number of coeff
          temp_line = txtfilereader.next_line(1);
          line_counter++;
          ValueTokenizer fp5th_line(temp_line);

          int num_coeff_2b = fp5th_line.next_int();
          if (num_coeff_2b <= 0)
            error->all(FLERR,
                       "UF3: 0 or negative number found for num_coeff_2b"
                       " on line {} of the potential file",
                       line_counter);
          n2b_coeff_array_size[itype][jtype] = num_coeff_2b;
          n2b_coeff_array_size[jtype][itype] = num_coeff_2b;
          max_num_coeff_2b = std::max(max_num_coeff_2b, num_coeff_2b);
        }
      } else if ((nbody_on_file == "3B") && (pot_3b)) {
        //3B block
        if (fp2nd_line.count() != 7)
          error->all(FLERR,
                     "UF3: Expected 7 words on line {} of {} file"
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
          for (int i = 1; i < num_of_elements + 1; i++) {
            if (std::string(elements[map[i]]) == element1) {
              itype = i;
              break;
            }
          }
          for (int i = 1; i < num_of_elements + 1; i++) {
            if (std::string(elements[map[i]]) == element2) {
              jtype = i;
              break;
            }
          }
          for (int i = 1; i < num_of_elements + 1; i++) {
            if (std::string(elements[map[i]]) == element3) {
              ktype = i;
              break;
            }
          }

          if ((itype != 0) && (jtype != 0) && (ktype != 0)) {
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

            cut_3b_list[itype][jtype] = std::max(cut3b_rij, cut_3b_list[itype][jtype]);
            cut_3b_list[itype][ktype] = std::max(cut_3b_list[itype][ktype], cut3b_rik);

            cut_3b[itype][jtype][ktype] = cut3b_rij;
            cut_3b[itype][ktype][jtype] = cut3b_rik;

            int num_knots_3b_jk = fp3rd_line.next_int();
            int num_knots_3b_ik = fp3rd_line.next_int();
            int num_knots_3b_ij = fp3rd_line.next_int();

            n3b_knots_array_size[map_3b[itype][jtype][ktype]][0] = num_knots_3b_jk;
            n3b_knots_array_size[map_3b[itype][jtype][ktype]][1] = num_knots_3b_ik;
            n3b_knots_array_size[map_3b[itype][jtype][ktype]][2] = num_knots_3b_ij;

            n3b_knots_array_size[map_3b[itype][ktype][jtype]][0] = num_knots_3b_jk;
            n3b_knots_array_size[map_3b[itype][ktype][jtype]][1] = num_knots_3b_ij;
            n3b_knots_array_size[map_3b[itype][ktype][jtype]][2] = num_knots_3b_ik;

            max_num_knots_3b = std::max(max_num_knots_3b, num_knots_3b_jk);
            max_num_knots_3b = std::max(max_num_knots_3b, num_knots_3b_ik);
            max_num_knots_3b = std::max(max_num_knots_3b, num_knots_3b_ij);

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
                         "SHAPE_OF_COEFF_MATRIX[I][J][K] found {} numbers",
                         fp7th_line.count());

            int coeff_matrix_dim1 = fp7th_line.next_int();
            int coeff_matrix_dim2 = fp7th_line.next_int();
            int coeff_matrix_dim3 = fp7th_line.next_int();

            n3b_coeff_array_size[map_3b[itype][jtype][ktype]][0] = coeff_matrix_dim1;
            n3b_coeff_array_size[map_3b[itype][jtype][ktype]][1] = coeff_matrix_dim2;
            n3b_coeff_array_size[map_3b[itype][jtype][ktype]][2] = coeff_matrix_dim3;

            n3b_coeff_array_size[map_3b[itype][ktype][jtype]][0] = coeff_matrix_dim2;
            n3b_coeff_array_size[map_3b[itype][ktype][jtype]][1] = coeff_matrix_dim1;
            n3b_coeff_array_size[map_3b[itype][ktype][jtype]][2] = coeff_matrix_dim3;

            max_num_coeff_3b = std::max(max_num_coeff_3b, coeff_matrix_dim1);
            max_num_coeff_3b = std::max(max_num_coeff_3b, coeff_matrix_dim2);
            max_num_coeff_3b = std::max(max_num_coeff_3b, coeff_matrix_dim3);
          }
        }
      } else {
        if (!((nbody_on_file == "3B") && (!pot_3b)))
          error->all(FLERR, "UF3: Expected either '2B' or '3B' word on line {} of {} file",
                     line_counter, potf_name);
      }
    }    //if of #UF3 POT
    line_counter++;
  }    // while

  //Create knot and coeff arrays
  if (max_num_knots_2b <= 0)
    error->all(FLERR,
               "UF3: Error reading the size of 2B knot vector\n"
               "Possibly no 2B UF3 potential block detected in {} file",
               potf_name);
  memory->destroy(n2b_knots_array);
  memory->create(n2b_knots_array, num_of_elements + 1, num_of_elements + 1, max_num_knots_2b,
                 "pair:n2b_knots_array");

  if (max_num_coeff_2b <= 0)
    error->all(FLERR,
               "UF3: Error reading the size of 2B coeff vector\n"
               "Possibly no 2B UF3 potential block detected in {} file",
               potf_name);

  memory->destroy(n2b_coeff_array);
  memory->create(n2b_coeff_array, num_of_elements + 1, num_of_elements + 1, max_num_coeff_2b,
                 "pair:n2b_coeff_array");

  if (pot_3b) {
    if (max_num_knots_3b <= 0)
      error->all(FLERR,
                 "UF3: Error reading the size of 3B knot vector\n"
                 "Possibly no 3B UF3 potential block detected in {} file",
                 potf_name);
    memory->destroy(n3b_knots_array);
    memory->create(n3b_knots_array, tot_interaction_count_3b, 3, max_num_knots_3b,
                   "pair:n3b_knots_array");

    if (max_num_coeff_3b <= 0)
      error->all(FLERR,
                 "UF3: Error reading the size of 3B coeff matrices\n"
                 "Possibly no 3B UF3 potential block detected in {} file",
                 potf_name);
    memory->destroy(n3b_coeff_array);
    memory->create(n3b_coeff_array, tot_interaction_count_3b, max_num_coeff_3b, max_num_coeff_3b,
                   max_num_coeff_3b, "pair:n3b_coeff_array");
  }

  //Go back to the begning of the file
  txtfilereader.rewind();

  //Go through the file again and fill knot and coeff arrays
  //while loop to read the data
  while ((line = txtfilereader.next_line(1))) {
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
        for (int i = 1; i < num_of_elements + 1; i++) {
          if (std::string(elements[map[i]]) == element1) {
            itype = i;
            break;
          }
        }
        for (int i = 1; i < num_of_elements + 1; i++) {
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
          knot_spacing_type_2b[itype][jtype] = 0;
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
            error->all(FLERR,
                       "UF3: Error reading the 2B potential block for {}-{}\n"
                       "Expected {} numbers on 4th line of the block but found {} numbers",
                       element1, element2, num_knots_2b, fp4th_line.count());

          for (int k = 0; k < num_knots_2b; k++) {
            n2b_knots_array[itype][jtype][k] = fp4th_line.next_double();
            n2b_knots_array[jtype][itype][k] = n2b_knots_array[itype][jtype][k];
          }

          knot_spacing_2b[itype][jtype] =
              n2b_knots_array[itype][jtype][4] - n2b_knots_array[itype][jtype][3];
          knot_spacing_2b[jtype][itype] = knot_spacing_2b[itype][jtype];

          //skip next line
          txtfilereader.skip_line();

          int num_of_coeff_2b = n2b_coeff_array_size[itype][jtype];

          temp_line = txtfilereader.next_line(num_of_coeff_2b);
          ValueTokenizer fp6th_line(temp_line);

          if (fp6th_line.count() != num_of_coeff_2b)
            error->all(FLERR,
                       "UF3: Error reading the 2B potential block for {}-{}\n"
                       "Expected {} numbers on 6th line of the block but found {} numbers",
                       element1, element2, num_of_coeff_2b, fp6th_line.count());

          for (int k = 0; k < num_of_coeff_2b; k++) {
            n2b_coeff_array[itype][jtype][k] = fp6th_line.next_double();
            n2b_coeff_array[jtype][itype][k] = n2b_coeff_array[itype][jtype][k];
          }

          if (num_knots_2b != num_of_coeff_2b + 4)
            error->all(FLERR,
                       "UF3: {}-{} interaction block has incorrect knot and "
                       "coeff data nknots (={}) != ncoeffs (={}) + 3 + 1",
                       element1, element2, num_knots_2b, num_of_coeff_2b);

          setflag[itype][jtype] = 1;
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
        for (int i = 1; i < num_of_elements + 1; i++) {
          if (std::string(elements[map[i]]) == element1) {
            itype = i;
            break;
          }
        }
        for (int i = 1; i < num_of_elements + 1; i++) {
          if (std::string(elements[map[i]]) == element2) {
            jtype = i;
            break;
          }
        }
        for (int i = 1; i < num_of_elements + 1; i++) {
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
          knot_spacing_type_3b[itype][jtype][ktype] = 0;
          knot_spacing_type_3b[itype][ktype][jtype] = 0;
        } else if (knot_type == "nk") {
          knot_spacing_type_3b[itype][jtype][ktype] = 1;
          knot_spacing_type_3b[itype][ktype][jtype] = 1;
        } else
          error->all(FLERR,
                     "UF3: Expected either 'uk'(uniform-knots) or 'nk'(non-uniform knots) "
                     "Found {} on the 2nd line of {}-{}-{} interaction block",
                     knot_type, element1, element2, element3);

        if ((itype != 0) && (jtype != 0) && (ktype != 0)) {
          //skip line containing info of cutoffs and knot vector sizes
          txtfilereader.skip_line();

          int num_knots_3b_jk = n3b_knots_array_size[map_3b[itype][jtype][ktype]][0];
          int num_knots_3b_ik = n3b_knots_array_size[map_3b[itype][jtype][ktype]][1];
          int num_knots_3b_ij = n3b_knots_array_size[map_3b[itype][jtype][ktype]][2];

          temp_line = txtfilereader.next_line(num_knots_3b_jk);
          ValueTokenizer fp4th_line(temp_line);
          if (fp4th_line.count() != num_knots_3b_jk)
            error->all(FLERR,
                       "UF3: Error reading the 3B potential block for {}-{}-{}\n"
                       "Expected {} numbers on 4th line of the block but found {} "
                       "numbers",
                       element1, element2, element3, num_knots_3b_jk, fp4th_line.count());

          for (int i = 0; i < num_knots_3b_jk; i++) {
            n3b_knots_array[map_3b[itype][jtype][ktype]][0][i] = fp4th_line.next_double();
            n3b_knots_array[map_3b[itype][ktype][jtype]][0][i] =
                n3b_knots_array[map_3b[itype][jtype][ktype]][0][i];
          }

          min_cut_3b[itype][jtype][ktype][0] = n3b_knots_array[map_3b[itype][jtype][ktype]][0][0];
          min_cut_3b[itype][ktype][jtype][0] = n3b_knots_array[map_3b[itype][ktype][jtype]][0][0];

          knot_spacing_3b[itype][jtype][ktype][0] =
              n3b_knots_array[map_3b[itype][jtype][ktype]][0][4] -
              n3b_knots_array[map_3b[itype][jtype][ktype]][0][3];
          knot_spacing_3b[itype][ktype][jtype][0] = knot_spacing_3b[itype][jtype][ktype][0];

          temp_line = txtfilereader.next_line(num_knots_3b_ik);
          ValueTokenizer fp5th_line(temp_line);
          if (fp5th_line.count() != num_knots_3b_ik)
            error->all(FLERR,
                       "UF3: Error reading the 3B potential block for {}-{}-{}\n"
                       "Expected {} numbers on 5th line of the block but found {} "
                       "numbers",
                       element1, element2, element3, num_knots_3b_ik, fp5th_line.count());

          for (int i = 0; i < num_knots_3b_ik; i++) {
            n3b_knots_array[map_3b[itype][jtype][ktype]][1][i] = fp5th_line.next_double();
            n3b_knots_array[map_3b[itype][ktype][jtype]][2][i] =
                n3b_knots_array[map_3b[itype][jtype][ktype]][1][i];
          }

          min_cut_3b[itype][jtype][ktype][1] = n3b_knots_array[map_3b[itype][jtype][ktype]][1][0];
          min_cut_3b[itype][ktype][jtype][2] = n3b_knots_array[map_3b[itype][ktype][jtype]][2][0];

          knot_spacing_3b[itype][jtype][ktype][1] =
              n3b_knots_array[map_3b[itype][jtype][ktype]][1][4] -
              n3b_knots_array[map_3b[itype][jtype][ktype]][1][3];
          knot_spacing_3b[itype][ktype][jtype][2] = knot_spacing_3b[itype][jtype][ktype][1];

          temp_line = txtfilereader.next_line(num_knots_3b_ij);
          ValueTokenizer fp6th_line(temp_line);
          if (fp6th_line.count() != num_knots_3b_ij)
            error->all(FLERR,
                       "UF3: Error reading the 3B potential block for {}-{}-{}\n"
                       "Expected {} numbers on 6th line of the block but found {} "
                       "numbers",
                       element1, element2, element3, num_knots_3b_ij, fp6th_line.count());

          for (int i = 0; i < num_knots_3b_ij; i++) {
            n3b_knots_array[map_3b[itype][jtype][ktype]][2][i] = fp6th_line.next_double();
            n3b_knots_array[map_3b[itype][ktype][jtype]][1][i] =
                n3b_knots_array[map_3b[itype][jtype][ktype]][2][i];
          }

          min_cut_3b[itype][jtype][ktype][2] = n3b_knots_array[map_3b[itype][jtype][ktype]][2][0];
          min_cut_3b[itype][ktype][jtype][1] = n3b_knots_array[map_3b[itype][ktype][jtype]][1][0];

          knot_spacing_3b[itype][jtype][ktype][2] =
              n3b_knots_array[map_3b[itype][jtype][ktype]][2][4] -
              n3b_knots_array[map_3b[itype][jtype][ktype]][2][3];
          knot_spacing_3b[itype][ktype][jtype][1] = knot_spacing_3b[itype][jtype][ktype][2];

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
                           "numbers",
                           element1, element2, element3, coeff_matrix_elements_len, line_count + 8,
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
                n3b_coeff_array[key2][j][i][k] = n3b_coeff_array[key1][i][j][k];
              }
            }
          }

          setflag_3b[itype][jtype][ktype] = 1;
          setflag_3b[itype][ktype][jtype] = 1;
        }
      }
    }    // if #UF3 POT
  }      //while
  fclose(fp);

  //Set interaction of atom types of the same elements
  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      if (setflag[i][j] != 1) {
        //i-j interaction not set

        //maybe i-j is mapped to some other atom type interaction?
        int i_mapped_to = map[i] + 1;    //+1 as map starts from 0
        int j_mapped_to = map[j] + 1;    //+1 as map starts from 0

        if ((i_mapped_to == i) && (j_mapped_to == j))
          //i-j is not mapped to some other atom type ie interaction is missing on file
          error->all(FLERR,
                     "UF3: Potential for interaction {}-{} ie {}-{} not found "
                     "in {} file",
                     i, j, elements[i_mapped_to - 1], elements[j_mapped_to - 1], potf_name);

        cut[i][j] = cut[i_mapped_to][j_mapped_to];

        n2b_knots_array_size[i][j] = n2b_knots_array_size[i_mapped_to][j_mapped_to];
        n2b_coeff_array_size[i][j] = n2b_coeff_array_size[i_mapped_to][j_mapped_to];

        knot_spacing_type_2b[i][j] = knot_spacing_type_2b[i_mapped_to][j_mapped_to];
        knot_spacing_2b[i][j] = knot_spacing_2b[i_mapped_to][j_mapped_to];

        for (int knot_no = 0; knot_no < max_num_knots_2b; knot_no++)
          n2b_knots_array[i][j][knot_no] = n2b_knots_array[i_mapped_to][j_mapped_to][knot_no];

        for (int coeff_no = 0; coeff_no < max_num_coeff_2b; coeff_no++)
          n2b_coeff_array[i][j][coeff_no] = n2b_coeff_array[i_mapped_to][j_mapped_to][coeff_no];

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
            int i_mapped_to = map[i] + 1;    //+1 as map starts from 0
            int j_mapped_to = map[j] + 1;    //+1 as map starts from 0
            int k_mapped_to = map[k] + 1;    //+1 as map starts from 0

            if ((i_mapped_to == i) && (j_mapped_to == j) && (k_mapped_to == k))
              error->all(FLERR,
                         "UF3: Potential for interaction {}-{}-{} ie {}-{}-{} "
                         " not found in {} file",
                         i, j, k, elements[i_mapped_to - 1], elements[j_mapped_to - 1],
                         elements[k_mapped_to - 1], potf_name);
            if (setflag_3b[i_mapped_to][j_mapped_to][k_mapped_to] != 1)
              error->all(FLERR,
                         "UF3: Interaction {}-{}-{} was mapped to {}-{}-{}, but "
                         "potential interaction for {}-{}-{} was not found in "
                         "{} file",
                         i, j, k, i_mapped_to, j_mapped_to, k_mapped_to, i_mapped_to, j_mapped_to,
                         k_mapped_to, potf_name);

            cut_3b_list[i][j] = std::max(cut_3b_list[i_mapped_to][j_mapped_to], cut_3b_list[i][j]);

            cut_3b[i][j][k] = cut_3b[i_mapped_to][j_mapped_to][k_mapped_to];

            knot_spacing_type_3b[i][j][k] =
                knot_spacing_type_3b[i_mapped_to][j_mapped_to][k_mapped_to];
            knot_spacing_3b[i][j][k][0] = knot_spacing_3b[i_mapped_to][j_mapped_to][k_mapped_to][0];
            knot_spacing_3b[i][j][k][1] = knot_spacing_3b[i_mapped_to][j_mapped_to][k_mapped_to][1];
            knot_spacing_3b[i][j][k][2] = knot_spacing_3b[i_mapped_to][j_mapped_to][k_mapped_to][2];

            int key = map_3b[i][j][k];
            int mapped_to_key = map_3b[i_mapped_to][j_mapped_to][k_mapped_to];

            n3b_knots_array_size[key][0] = n3b_knots_array_size[mapped_to_key][0];
            n3b_knots_array_size[key][1] = n3b_knots_array_size[mapped_to_key][1];
            n3b_knots_array_size[key][2] = n3b_knots_array_size[mapped_to_key][2];

            n3b_coeff_array_size[key][0] = n3b_coeff_array_size[mapped_to_key][0];
            n3b_coeff_array_size[key][1] = n3b_coeff_array_size[mapped_to_key][1];
            n3b_coeff_array_size[key][2] = n3b_coeff_array_size[mapped_to_key][2];

            min_cut_3b[i][j][k][0] = min_cut_3b[i_mapped_to][j_mapped_to][k_mapped_to][0];

            min_cut_3b[i][j][k][1] = min_cut_3b[i_mapped_to][j_mapped_to][k_mapped_to][1];

            min_cut_3b[i][j][k][2] = min_cut_3b[i_mapped_to][j_mapped_to][k_mapped_to][2];

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

//Broadcast data read from potential file to all processors
void PairUF3::communicate()
{
  const int num_of_elements = atom->ntypes;
  MPI_Bcast(&cut[0][0], (num_of_elements + 1) * (num_of_elements + 1), MPI_DOUBLE, 0, world);

  MPI_Bcast(&n2b_knots_array_size[0][0], (num_of_elements + 1) * (num_of_elements + 1), MPI_INT, 0,
            world);
  MPI_Bcast(&n2b_coeff_array_size[0][0], (num_of_elements + 1) * (num_of_elements + 1), MPI_INT, 0,
            world);

  MPI_Bcast(&max_num_knots_2b, 1, MPI_INT, 0, world);
  MPI_Bcast(&max_num_coeff_2b, 1, MPI_INT, 0, world);

  if (pot_3b) {
    MPI_Bcast(&cut_3b_list[0][0], (num_of_elements + 1) * (num_of_elements + 1), MPI_DOUBLE, 0,
              world);

    MPI_Bcast(&cut_3b[0][0][0],
              (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1), MPI_DOUBLE, 0,
              world);

    MPI_Bcast(&n3b_knots_array_size[0][0], tot_interaction_count_3b * 3, MPI_INT, 0, world);
    MPI_Bcast(&n3b_coeff_array_size[0][0], tot_interaction_count_3b * 3, MPI_INT, 0, world);

    MPI_Bcast(&max_num_knots_3b, 1, MPI_INT, 0, world);
    MPI_Bcast(&max_num_coeff_3b, 1, MPI_INT, 0, world);
  }

  if (comm->me != 0) {
    memory->destroy(n2b_knots_array);
    memory->destroy(n2b_coeff_array);

    memory->create(n2b_knots_array, num_of_elements + 1, num_of_elements + 1, max_num_knots_2b,
                   "pair:n2b_knots_array");
    memory->create(n2b_coeff_array, num_of_elements + 1, num_of_elements + 1, max_num_coeff_2b,
                   "pair:n2b_coeff_array");
    if (pot_3b) {
      memory->destroy(n3b_knots_array);
      memory->destroy(n3b_coeff_array);

      memory->create(n3b_knots_array, tot_interaction_count_3b, 3, max_num_knots_3b,
                     "pair:n3b_knots_array");

      memory->create(n3b_coeff_array, tot_interaction_count_3b, max_num_coeff_3b, max_num_coeff_3b,
                     max_num_coeff_3b, "pair:n3b_coeff_array");
    }
  }

  MPI_Bcast(&knot_spacing_type_2b[0][0], (num_of_elements + 1) * (num_of_elements + 1), MPI_INT, 0,
            world);

  MPI_Bcast(&knot_spacing_2b[0][0], (num_of_elements + 1) * (num_of_elements + 1), MPI_DOUBLE, 0,
            world);

  MPI_Bcast(&n2b_knots_array[0][0][0],
            (num_of_elements + 1) * (num_of_elements + 1) * max_num_knots_2b, MPI_DOUBLE, 0, world);
  MPI_Bcast(&n2b_coeff_array[0][0][0],
            (num_of_elements + 1) * (num_of_elements + 1) * max_num_coeff_2b, MPI_DOUBLE, 0, world);

  MPI_Bcast(&setflag[0][0], (num_of_elements + 1) * (num_of_elements + 1), MPI_INT, 0, world);

  if (pot_3b) {
    MPI_Bcast(&knot_spacing_type_3b[0][0][0],
              (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1), MPI_INT, 0,
              world);

    MPI_Bcast(&knot_spacing_3b[0][0][0][0],
              (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) * 3, MPI_DOUBLE,
              0, world);
    MPI_Bcast(&n3b_knots_array[0][0][0], tot_interaction_count_3b * 3 * max_num_knots_3b,
              MPI_DOUBLE, 0, world);
    MPI_Bcast(&n3b_coeff_array[0][0][0][0],
              tot_interaction_count_3b * max_num_coeff_3b * max_num_coeff_3b * max_num_coeff_3b,
              MPI_DOUBLE, 0, world);
    MPI_Bcast(&setflag_3b[0][0][0],
              (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1), MPI_INT, 0,
              world);
    MPI_Bcast(&min_cut_3b[0][0][0][0],
              (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) * 3, MPI_DOUBLE,
              0, world);
  }
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
  int spacing_type = knot_spacing_type_2b[1][1];
  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      if (setflag[i][j] != 1)
        error->all(FLERR,
                   "UF3: Not all 2-body UF potentials are set, "
                   "missing potential for {}-{} interaction",
                   i, j);
      /*if (spacing_type != knot_spacing_type_2b[i][j])
        error->all(FLERR,
                   "UF3: In the current version the knot spacing type, "
                   "for all interactions needs to be same. For {}-{} "
                   "i.e. {}-{} interaction expected {}, but found {}",
                   i,j,elements[map[i]],elements[map[j]],spacing_type,
                   knot_spacing_type_2b[i][j]);*/
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
          if (spacing_type != knot_spacing_type_3b[i][j][k])
            error->all(FLERR,
                       "UF3: In the current version the knot spacing type, "
                       "for all interactions needs to be same. For {}-{}-{} "
                       "i.e. {}-{}-{} interaction expected{}, but found {}",
                       i, j, k, elements[map[i]], elements[map[j]], elements[map[k]], spacing_type,
                       knot_spacing_type_3b[i][j][k]);
        }
      }
    }
  }

  if (spacing_type) {
    get_starting_index_2b = &PairUF3::get_starting_index_nonuniform_2b;
    if (pot_3b) get_starting_index_3b = &PairUF3::get_starting_index_nonuniform_3b;
  } else {
    get_starting_index_2b = &PairUF3::get_starting_index_uniform_2b;
    if (pot_3b) get_starting_index_3b = &PairUF3::get_starting_index_uniform_3b;
  }

  create_cached_constants_2b();
  if (pot_3b) create_cached_constants_3b();
}

int PairUF3::get_starting_index_uniform_2b(int i, int j, double r)
{
  return 3 + (int) ((r - n2b_knots_array[i][j][0]) / (knot_spacing_2b[i][j]));
}

int PairUF3::get_starting_index_uniform_3b(int i, int j, int k, double r, int knot_dim)
{
  return 3 +
      (int) (((r - n3b_knots_array[map_3b[i][j][k]][knot_dim][0]) /
              knot_spacing_3b[i][j][k][knot_dim]));
}

int PairUF3::get_starting_index_nonuniform_2b(int i, int j, double r)
{
  for (int l = 3; l < n2b_knots_array_size[i][j] - 1; ++l) {
    if ((n2b_knots_array[i][j][l] <= r) && (r < n2b_knots_array[i][j][l + 1])) return l;
  }
  return -1;
}

int PairUF3::get_starting_index_nonuniform_3b(int i, int j, int k, double r, int knot_dim)
{
  for (int l = 3; l < n3b_knots_array_size[map_3b[i][j][k]][knot_dim] - 1; ++l) {
    if ((n3b_knots_array[map_3b[i][j][k]][knot_dim][l] <= r) &&
        (r < n3b_knots_array[map_3b[i][j][k]][knot_dim][l + 1]))
      return l;
  }
  return -1;
}

void PairUF3::create_cached_constants_2b()
{
  const int num_of_elements = atom->ntypes;
  memory->destroy(cached_constants_2b);
  memory->destroy(cached_constants_2b_deri);
  memory->create(cached_constants_2b, num_of_elements + 1, num_of_elements + 1, max_num_coeff_2b,
                 16, "pair:cached_constants_2b");

  memory->create(cached_constants_2b_deri, num_of_elements + 1, num_of_elements + 1,
                 max_num_coeff_2b - 1, 9, "pair:cached_constants_2b_deri");

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      for (int l = 0; l < n2b_coeff_array_size[i][j]; l++) {
        uf3_bspline_basis3 bspline_basis(lmp, &n2b_knots_array[i][j][l], n2b_coeff_array[i][j][l]);
        for (int cc = 0; cc < 16; cc++) {
          cached_constants_2b[i][j][l][cc] = bspline_basis.constants[cc];
        }
      }
    }
  }

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      //initialize coeff and knots for derivative
      //double* knots_for_deri = new double[n2b_knots_array_size[i][j]-2];
      double *knots_for_deri = nullptr;
      memory->create(knots_for_deri, n2b_knots_array_size[i][j] - 2, "pair:knots_for_deri");

      for (int l = 1; l < n2b_knots_array_size[i][j] - 1; l++)
        knots_for_deri[l - 1] = n2b_knots_array[i][j][l];

      //double* coeff_for_deri = new double[n2b_coeff_array_size[i][j]-1];
      double *coeff_for_deri = nullptr;
      memory->create(coeff_for_deri, n2b_coeff_array_size[i][j] - 1, "pair:coeff_for_deri");
      for (int l = 0; l < n2b_coeff_array_size[i][j] - 1; l++) {
        double dntemp = 3 / (n2b_knots_array[i][j][l + 4] - n2b_knots_array[i][j][l + 1]);
        coeff_for_deri[l] = (n2b_coeff_array[i][j][l + 1] - n2b_coeff_array[i][j][l]) * dntemp;
      }

      for (int l = 0; l < n2b_coeff_array_size[i][j] - 1; l++) {
        uf3_bspline_basis2 bspline_basis_deri(lmp, &knots_for_deri[l], coeff_for_deri[l]);
        for (int cc = 0; cc < 9; cc++) {
          cached_constants_2b_deri[i][j][l][cc] = bspline_basis_deri.constants[cc];
        }
      }
      memory->destroy(knots_for_deri);
      memory->destroy(coeff_for_deri);
      //delete[] knots_for_deri;
      //delete[] coeff_for_deri;
    }
  }
}

void PairUF3::create_cached_constants_3b()
{
  const int num_of_elements = atom->ntypes;
  memory->destroy(coeff_for_der_jk);
  memory->destroy(coeff_for_der_ik);
  memory->destroy(coeff_for_der_ij);
  memory->destroy(cached_constants_3b);
  memory->destroy(cached_constants_3b_deri);

  memory->create(coeff_for_der_jk, tot_interaction_count_3b, max_num_coeff_3b, max_num_coeff_3b,
                 max_num_coeff_3b, "pair:coeff_for_der_jk");

  memory->create(coeff_for_der_ik, tot_interaction_count_3b, max_num_coeff_3b, max_num_coeff_3b,
                 max_num_coeff_3b, "pair:coeff_for_der_ik");

  memory->create(coeff_for_der_ij, tot_interaction_count_3b, max_num_coeff_3b, max_num_coeff_3b,
                 max_num_coeff_3b, "pair:coeff_for_der_ij");

  memory->create(cached_constants_3b, tot_interaction_count_3b, 3, max_num_coeff_3b, 16,
                 "pair:cached_constants_3b");

  memory->create(cached_constants_3b_deri, tot_interaction_count_3b, 3, max_num_coeff_3b - 1, 9,
                 "pair:cached_constants_3b_deri");

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      for (int k = 1; k < num_of_elements + 1; k++) {
        int map_to = map_3b[i][j][k];

        for (int l = 0; l < n3b_knots_array_size[map_to][2] - 4; l++) {
          uf3_bspline_basis3 bspline_basis_ij(lmp, &n3b_knots_array[map_to][2][l], 1);
          for (int cc = 0; cc < 16; cc++)
            cached_constants_3b[map_to][0][l][cc] = bspline_basis_ij.constants[cc];
        }

        for (int l = 0; l < n3b_knots_array_size[map_to][1] - 4; l++) {
          uf3_bspline_basis3 bspline_basis_ik(lmp, &n3b_knots_array[map_to][1][l], 1);
          for (int cc = 0; cc < 16; cc++)
            cached_constants_3b[map_to][1][l][cc] = bspline_basis_ik.constants[cc];
        }

        for (int l = 0; l < n3b_knots_array_size[map_to][0] - 4; l++) {
          uf3_bspline_basis3 bspline_basis_jk(lmp, &n3b_knots_array[map_to][0][l], 1);
          for (int cc = 0; cc < 16; cc++)
            cached_constants_3b[map_to][2][l][cc] = bspline_basis_jk.constants[cc];
        }
      }
    }
  }

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      for (int k = 1; k < num_of_elements + 1; k++) {
        int map_to = map_3b[i][j][k];
        double **knots_for_der = nullptr;    //new double*[3];

        //n3b_knots_array_size[map_to][0] for jk knot vector --> always largest
        memory->create(knots_for_der, 3, n3b_knots_array_size[map_to][0] - 1, "pair:knots_for_der");

        //--deri_basis_jk
        for (int l = 1; l < n3b_knots_array_size[map_to][0] - 1; l++)
          knots_for_der[0][l - 1] = n3b_knots_array[map_to][0][l];

        for (int l = 0; l < n3b_coeff_array_size[map_to][0]; l++) {
          for (int m = 0; m < n3b_coeff_array_size[map_to][1]; m++) {
            for (int n = 0; n < n3b_coeff_array_size[map_to][2] - 1; n++) {
              double dntemp =
                  3 / (n3b_knots_array[map_to][0][n + 4] - n3b_knots_array[map_to][0][n + 1]);
              coeff_for_der_jk[map_to][l][m][n] =
                  ((n3b_coeff_array[map_to][l][m][n + 1] - n3b_coeff_array[map_to][l][m][n]) *
                   dntemp);
            }
          }
        }

        //--deri_basis_ik
        for (int l = 1; l < n3b_knots_array_size[map_to][1] - 1; l++)
          knots_for_der[1][l - 1] = n3b_knots_array[map_to][1][l];

        for (int l = 0; l < n3b_coeff_array_size[map_to][0]; l++) {
          for (int m = 0; m < n3b_coeff_array_size[map_to][1] - 1; m++) {
            double dntemp =
                3 / (n3b_knots_array[map_to][1][m + 4] - n3b_knots_array[map_to][1][m + 1]);
            for (int n = 0; n < n3b_coeff_array_size[map_to][2]; n++) {
              coeff_for_der_ik[map_to][l][m][n] =
                  ((n3b_coeff_array[map_to][l][m + 1][n] - n3b_coeff_array[map_to][l][m][n]) *
                   dntemp);
            }
          }
        }

        //--deri_basis_ij
        for (int l = 1; l < n3b_knots_array_size[map_to][2] - 1; l++)
          knots_for_der[2][l - 1] = n3b_knots_array[map_to][2][l];

        for (int l = 0; l < n3b_coeff_array_size[map_to][0] - 1; l++) {
          double dntemp =
              3 / (n3b_knots_array[map_to][2][l + 4] - n3b_knots_array[map_to][2][l + 1]);
          for (int m = 0; m < n3b_coeff_array_size[map_to][1]; m++) {
            for (int n = 0; n < n3b_coeff_array_size[map_to][2]; n++) {
              coeff_for_der_ij[map_to][l][m][n] =
                  ((n3b_coeff_array[map_to][l + 1][m][n] - n3b_coeff_array[map_to][l][m][n]) *
                   dntemp);
            }
          }
        }

        for (int l = 0; l < n3b_coeff_array_size[map_to][0] - 1; l++) {
          uf3_bspline_basis2 bspline_basis_deri_ij(lmp, &knots_for_der[2][l], 1);
          for (int cc = 0; cc < 9; cc++) {
            cached_constants_3b_deri[map_to][0][l][cc] = bspline_basis_deri_ij.constants[cc];
          }
        }

        for (int l = 0; l < n3b_coeff_array_size[map_to][1] - 1; l++) {
          uf3_bspline_basis2 bspline_basis_deri_ik(lmp, &knots_for_der[1][l], 1);
          for (int cc = 0; cc < 9; cc++) {
            cached_constants_3b_deri[map_to][1][l][cc] = bspline_basis_deri_ik.constants[cc];
          }
        }

        for (int l = 0; l < n3b_coeff_array_size[map_to][2] - 1; l++) {
          uf3_bspline_basis2 bspline_basis_deri_jk(lmp, &knots_for_der[0][l], 1);
          for (int cc = 0; cc < 9; cc++) {
            cached_constants_3b_deri[map_to][2][l][cc] = bspline_basis_deri_jk.constants[cc];
          }
        }

        memory->destroy(knots_for_der);
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
  double rij_sq, rik_sq, rjk_sq;
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

        int knot_start_index = (this->*get_starting_index_2b)(itype, jtype, rij);

        double force_2b = cached_constants_2b_deri[itype][jtype][knot_start_index - 1][0];
        force_2b += rij * cached_constants_2b_deri[itype][jtype][knot_start_index - 1][1];
        force_2b += rsq * cached_constants_2b_deri[itype][jtype][knot_start_index - 1][2];
        force_2b += cached_constants_2b_deri[itype][jtype][knot_start_index - 2][3];
        force_2b += rij * cached_constants_2b_deri[itype][jtype][knot_start_index - 2][4];
        force_2b += rsq * cached_constants_2b_deri[itype][jtype][knot_start_index - 2][5];
        force_2b += cached_constants_2b_deri[itype][jtype][knot_start_index - 3][6];
        force_2b += rij * cached_constants_2b_deri[itype][jtype][knot_start_index - 3][7];
        force_2b += rsq * cached_constants_2b_deri[itype][jtype][knot_start_index - 3][8];

        fpair = -1 * force_2b / rij;

        fx = delx * fpair;
        fy = dely * fpair;
        fz = delz * fpair;

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;

        if (eflag) {
          double rth = rsq * rij;
          evdwl = cached_constants_2b[itype][jtype][knot_start_index][0];
          evdwl += rij * cached_constants_2b[itype][jtype][knot_start_index][1];
          evdwl += rsq * cached_constants_2b[itype][jtype][knot_start_index][2];
          evdwl += rth * cached_constants_2b[itype][jtype][knot_start_index][3];
          evdwl += cached_constants_2b[itype][jtype][knot_start_index - 1][4];
          evdwl += rij * cached_constants_2b[itype][jtype][knot_start_index - 1][5];
          evdwl += rsq * cached_constants_2b[itype][jtype][knot_start_index - 1][6];
          evdwl += rth * cached_constants_2b[itype][jtype][knot_start_index - 1][7];
          evdwl += cached_constants_2b[itype][jtype][knot_start_index - 2][8];
          evdwl += rij * cached_constants_2b[itype][jtype][knot_start_index - 2][9];
          evdwl += rsq * cached_constants_2b[itype][jtype][knot_start_index - 2][10];
          evdwl += rth * cached_constants_2b[itype][jtype][knot_start_index - 2][11];
          evdwl += cached_constants_2b[itype][jtype][knot_start_index - 3][12];
          evdwl += rij * cached_constants_2b[itype][jtype][knot_start_index - 3][13];
          evdwl += rsq * cached_constants_2b[itype][jtype][knot_start_index - 3][14];
          evdwl += rth * cached_constants_2b[itype][jtype][knot_start_index - 3][15];
        };

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
      rij_sq = (del_rji[0] * del_rji[0]) + (del_rji[1] * del_rji[1]) + (del_rji[2] * del_rji[2]);
      rij = sqrt(rij_sq);

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
        rik_sq = (del_rki[0] * del_rki[0]) + (del_rki[1] * del_rki[1]) + (del_rki[2] * del_rki[2]);
        rik = sqrt(rik_sq);

        if ((rij <= cut_3b[itype][jtype][ktype]) && (rik <= cut_3b[itype][ktype][jtype]) &&
            (rij >= min_cut_3b[itype][jtype][ktype][2]) &&
            (rik >= min_cut_3b[itype][jtype][ktype][1])) {

          del_rkj[0] = x[k][0] - x[j][0];
          del_rkj[1] = x[k][1] - x[j][1];
          del_rkj[2] = x[k][2] - x[j][2];

          rjk_sq =
              (del_rkj[0] * del_rkj[0]) + (del_rkj[1] * del_rkj[1]) + (del_rkj[2] * del_rkj[2]);
          rjk = sqrt(rjk_sq);

          if (rjk >= min_cut_3b[itype][jtype][ktype][0]) {
            double rij_th = rij * rij_sq;
            double rik_th = rik * rik_sq;
            double rjk_th = rjk * rjk_sq;

            int map_to = map_3b[itype][jtype][ktype];
            int knot_start_index_ij = (this->*get_starting_index_3b)(itype, jtype, ktype, rij, 2);
            int knot_start_index_ik = (this->*get_starting_index_3b)(itype, jtype, ktype, rik, 1);
            int knot_start_index_jk = (this->*get_starting_index_3b)(itype, jtype, ktype, rjk, 0);
            double basis_ij[4];
            double basis_ik[4];
            double basis_jk[4];
            double basis_ij_der[3];
            double basis_ik_der[3];
            double basis_jk_der[3];

            //--------------basis_ij
            basis_ij[0] = cached_constants_3b[map_to][0][knot_start_index_ij - 3][12];
            basis_ij[0] += rij * cached_constants_3b[map_to][0][knot_start_index_ij - 3][13];
            basis_ij[0] += rij_sq * cached_constants_3b[map_to][0][knot_start_index_ij - 3][14];
            basis_ij[0] += rij_th * cached_constants_3b[map_to][0][knot_start_index_ij - 3][15];

            basis_ij[1] = cached_constants_3b[map_to][0][knot_start_index_ij - 2][8];
            basis_ij[1] += rij * cached_constants_3b[map_to][0][knot_start_index_ij - 2][9];
            basis_ij[1] += rij_sq * cached_constants_3b[map_to][0][knot_start_index_ij - 2][10];
            basis_ij[1] += rij_th * cached_constants_3b[map_to][0][knot_start_index_ij - 2][11];

            basis_ij[2] = cached_constants_3b[map_to][0][knot_start_index_ij - 1][4];
            basis_ij[2] += rij * cached_constants_3b[map_to][0][knot_start_index_ij - 1][5];
            basis_ij[2] += rij_sq * cached_constants_3b[map_to][0][knot_start_index_ij - 1][6];
            basis_ij[2] += rij_th * cached_constants_3b[map_to][0][knot_start_index_ij - 1][7];

            basis_ij[3] = cached_constants_3b[map_to][0][knot_start_index_ij][0];
            basis_ij[3] += rij * cached_constants_3b[map_to][0][knot_start_index_ij][1];
            basis_ij[3] += rij_sq * cached_constants_3b[map_to][0][knot_start_index_ij][2];
            basis_ij[3] += rij_th * cached_constants_3b[map_to][0][knot_start_index_ij][3];

            //--------------basis_ik
            basis_ik[0] = cached_constants_3b[map_to][1][knot_start_index_ik - 3][12];
            basis_ik[0] += rik * cached_constants_3b[map_to][1][knot_start_index_ik - 3][13];
            basis_ik[0] += rik_sq * cached_constants_3b[map_to][1][knot_start_index_ik - 3][14];
            basis_ik[0] += rik_th * cached_constants_3b[map_to][1][knot_start_index_ik - 3][15];

            basis_ik[1] = cached_constants_3b[map_to][1][knot_start_index_ik - 2][8];
            basis_ik[1] += rik * cached_constants_3b[map_to][1][knot_start_index_ik - 2][9];
            basis_ik[1] += rik_sq * cached_constants_3b[map_to][1][knot_start_index_ik - 2][10];
            basis_ik[1] += rik_th * cached_constants_3b[map_to][1][knot_start_index_ik - 2][11];

            basis_ik[2] = cached_constants_3b[map_to][1][knot_start_index_ik - 1][4];
            basis_ik[2] += rik * cached_constants_3b[map_to][1][knot_start_index_ik - 1][5];
            basis_ik[2] += rik_sq * cached_constants_3b[map_to][1][knot_start_index_ik - 1][6];
            basis_ik[2] += rik_th * cached_constants_3b[map_to][1][knot_start_index_ik - 1][7];

            basis_ik[3] = cached_constants_3b[map_to][1][knot_start_index_ik][0];
            basis_ik[3] += rik * cached_constants_3b[map_to][1][knot_start_index_ik][1];
            basis_ik[3] += rik_sq * cached_constants_3b[map_to][1][knot_start_index_ik][2];
            basis_ik[3] += rik_th * cached_constants_3b[map_to][1][knot_start_index_ik][3];

            //--------------basis_jk
            basis_jk[0] = cached_constants_3b[map_to][2][knot_start_index_jk - 3][12];
            basis_jk[0] += rjk * cached_constants_3b[map_to][2][knot_start_index_jk - 3][13];
            basis_jk[0] += rjk_sq * cached_constants_3b[map_to][2][knot_start_index_jk - 3][14];
            basis_jk[0] += rjk_th * cached_constants_3b[map_to][2][knot_start_index_jk - 3][15];

            basis_jk[1] = cached_constants_3b[map_to][2][knot_start_index_jk - 2][8];
            basis_jk[1] += rjk * cached_constants_3b[map_to][2][knot_start_index_jk - 2][9];
            basis_jk[1] += rjk_sq * cached_constants_3b[map_to][2][knot_start_index_jk - 2][10];
            basis_jk[1] += rjk_th * cached_constants_3b[map_to][2][knot_start_index_jk - 2][11];

            basis_jk[2] = cached_constants_3b[map_to][2][knot_start_index_jk - 1][4];
            basis_jk[2] += rjk * cached_constants_3b[map_to][2][knot_start_index_jk - 1][5];
            basis_jk[2] += rjk_sq * cached_constants_3b[map_to][2][knot_start_index_jk - 1][6];
            basis_jk[2] += rjk_th * cached_constants_3b[map_to][2][knot_start_index_jk - 1][7];

            basis_jk[3] = cached_constants_3b[map_to][2][knot_start_index_jk][0];
            basis_jk[3] += rjk * cached_constants_3b[map_to][2][knot_start_index_jk][1];
            basis_jk[3] += rjk_sq * cached_constants_3b[map_to][2][knot_start_index_jk][2];
            basis_jk[3] += rjk_th * cached_constants_3b[map_to][2][knot_start_index_jk][3];

            //----------------basis_ij_der
            basis_ij_der[0] = cached_constants_3b_deri[map_to][0][knot_start_index_ij - 3][6];
            basis_ij_der[0] +=
                rij * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 3][7];
            basis_ij_der[0] +=
                rij_sq * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 3][8];

            basis_ij_der[1] = cached_constants_3b_deri[map_to][0][knot_start_index_ij - 2][3];
            basis_ij_der[1] +=
                rij * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 2][4];
            basis_ij_der[1] +=
                rij_sq * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 2][5];

            basis_ij_der[2] = cached_constants_3b_deri[map_to][0][knot_start_index_ij - 1][0];
            basis_ij_der[2] +=
                rij * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 1][1];
            basis_ij_der[2] +=
                rij_sq * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 1][2];

            //----------------basis_ik_der
            basis_ik_der[0] = cached_constants_3b_deri[map_to][1][knot_start_index_ik - 3][6];
            basis_ik_der[0] +=
                rik * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 3][7];
            basis_ik_der[0] +=
                rik_sq * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 3][8];

            basis_ik_der[1] = cached_constants_3b_deri[map_to][1][knot_start_index_ik - 2][3];
            basis_ik_der[1] +=
                rik * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 2][4];
            basis_ik_der[1] +=
                rik_sq * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 2][5];

            basis_ik_der[2] = cached_constants_3b_deri[map_to][1][knot_start_index_ik - 1][0];
            basis_ik_der[2] +=
                rik * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 1][1];
            basis_ik_der[2] +=
                rik_sq * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 1][2];

            //----------------basis_jk_der
            basis_jk_der[0] = cached_constants_3b_deri[map_to][2][knot_start_index_jk - 3][6];
            basis_jk_der[0] +=
                rjk * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 3][7];
            basis_jk_der[0] +=
                rjk_sq * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 3][8];

            basis_jk_der[1] = cached_constants_3b_deri[map_to][2][knot_start_index_jk - 2][3];
            basis_jk_der[1] +=
                rjk * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 2][4];
            basis_jk_der[1] +=
                rjk_sq * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 2][5];

            basis_jk_der[2] = cached_constants_3b_deri[map_to][2][knot_start_index_jk - 1][0];
            basis_jk_der[2] +=
                rjk * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 1][1];
            basis_jk_der[2] +=
                rjk_sq * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 1][2];

            double triangle_eval[4] = {0, 0, 0, 0};

            int iknot_ij = knot_start_index_ij - 3;
            int iknot_ik = knot_start_index_ik - 3;
            int iknot_jk = knot_start_index_jk - 3;

            for (int l = 0; l < 3; l++) {
              const double basis_ij_der_i = basis_ij_der[l];
              for (int m = 0; m < 4; m++) {
                const double factor = basis_ij_der_i * basis_ik[m];
                const double *slice =
                    &coeff_for_der_ij[map_to][iknot_ij + l][iknot_ik + m][iknot_jk];
                double tmp[4];
                tmp[0] = slice[0] * basis_jk[0];
                tmp[1] = slice[1] * basis_jk[1];
                tmp[2] = slice[2] * basis_jk[2];
                tmp[3] = slice[3] * basis_jk[3];
                double sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];
                triangle_eval[1] += factor * sum;
              }
            }

            for (int l = 0; l < 4; l++) {
              const double basis_ij_i = basis_ij[l];
              for (int m = 0; m < 3; m++) {
                const double factor = basis_ij_i * basis_ik_der[m];
                const double *slice =
                    &coeff_for_der_ik[map_to][iknot_ij + l][iknot_ik + m][iknot_jk];
                double tmp[4];
                tmp[0] = slice[0] * basis_jk[0];
                tmp[1] = slice[1] * basis_jk[1];
                tmp[2] = slice[2] * basis_jk[2];
                tmp[3] = slice[3] * basis_jk[3];
                double sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];
                triangle_eval[2] += factor * sum;
              }
            }

            for (int l = 0; l < 4; l++) {
              const double basis_ij_i = basis_ij[l];
              for (int m = 0; m < 4; m++) {
                const double factor = basis_ij_i * basis_ik[m];
                const double *slice =
                    &coeff_for_der_jk[map_to][iknot_ij + l][iknot_ik + m][iknot_jk];
                double tmp[3];
                tmp[0] = slice[0] * basis_jk_der[0];
                tmp[1] = slice[1] * basis_jk_der[1];
                tmp[2] = slice[2] * basis_jk_der[2];
                double sum = tmp[0] + tmp[1] + tmp[2];
                triangle_eval[3] += factor * sum;
              }
            }

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

            if (eflag) {
              for (int l = 0; l < 4; l++) {
                const double basis_ij_i = basis_ij[l];
                for (int m = 0; m < 4; m++) {
                  const double factor = basis_ij_i * basis_ik[m];
                  const double *slice =
                      &n3b_coeff_array[map_to][iknot_ij + l][iknot_ik + m][iknot_jk];
                  double tmp[4];
                  tmp[0] = slice[0] * basis_jk[0];
                  tmp[1] = slice[1] * basis_jk[1];
                  tmp[2] = slice[2] * basis_jk[2];
                  tmp[3] = slice[3] * basis_jk[3];
                  double sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];
                  triangle_eval[0] += factor * sum;
                }
              }
              evdwl = *triangle_eval;
            }

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
    int knot_start_index = (this->*get_starting_index_2b)(itype, jtype, r);

    double force_2b = cached_constants_2b_deri[itype][jtype][knot_start_index - 1][0];
    force_2b += r * cached_constants_2b_deri[itype][jtype][knot_start_index - 1][1];
    force_2b += rsq * cached_constants_2b_deri[itype][jtype][knot_start_index - 1][2];
    force_2b += cached_constants_2b_deri[itype][jtype][knot_start_index - 2][3];
    force_2b += r * cached_constants_2b_deri[itype][jtype][knot_start_index - 2][4];
    force_2b += rsq * cached_constants_2b_deri[itype][jtype][knot_start_index - 2][5];
    force_2b += cached_constants_2b_deri[itype][jtype][knot_start_index - 3][6];
    force_2b += r * cached_constants_2b_deri[itype][jtype][knot_start_index - 3][7];
    force_2b += rsq * cached_constants_2b_deri[itype][jtype][knot_start_index - 3][8];
    fforce = factor_lj * force_2b;

    double rth = rsq * r;
    value = cached_constants_2b[itype][jtype][knot_start_index][0];
    value += r * cached_constants_2b[itype][jtype][knot_start_index][1];
    value += rsq * cached_constants_2b[itype][jtype][knot_start_index][2];
    value += rth * cached_constants_2b[itype][jtype][knot_start_index][3];
    value += cached_constants_2b[itype][jtype][knot_start_index - 1][4];
    value += r * cached_constants_2b[itype][jtype][knot_start_index - 1][5];
    value += rsq * cached_constants_2b[itype][jtype][knot_start_index - 1][6];
    value += rth * cached_constants_2b[itype][jtype][knot_start_index - 1][7];
    value += cached_constants_2b[itype][jtype][knot_start_index - 2][8];
    value += r * cached_constants_2b[itype][jtype][knot_start_index - 2][9];
    value += rsq * cached_constants_2b[itype][jtype][knot_start_index - 2][10];
    value += rth * cached_constants_2b[itype][jtype][knot_start_index - 2][11];
    value += cached_constants_2b[itype][jtype][knot_start_index - 3][12];
    value += r * cached_constants_2b[itype][jtype][knot_start_index - 3][13];
    value += rsq * cached_constants_2b[itype][jtype][knot_start_index - 3][14];
    value += rth * cached_constants_2b[itype][jtype][knot_start_index - 3][15];
  }

  return factor_lj * value;
}

double PairUF3::memory_usage()
{
  const int num_of_elements = atom->ntypes;
  double bytes = Pair::memory_usage();

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) *
      sizeof(int);    //***setflag_3b

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) *
      sizeof(int);    //knot_spacing_type_2b
  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) *
      sizeof(int);    //knot_spacing_type_3b

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * sizeof(double);    //cut

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) *
      sizeof(double);    //***cut_3b

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * sizeof(double);    //cut_3b_list

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) * 3 *
      sizeof(double);    //min_cut_3b

  bytes +=
      (double) (num_of_elements + 1) * (num_of_elements + 1) * sizeof(double);    //knot_spacing_2b
  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) *
      sizeof(double);    //knot_spacing_3b

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * max_num_knots_2b *
      sizeof(double);    //n2b_knots_array

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * max_num_coeff_2b *
      sizeof(double);    //n2b_coeff_array

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) *
      sizeof(int);    //n2b_knots_array_size
  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) *
      sizeof(int);    //n2b_coeff_array_size

  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * max_num_coeff_2b * 16 *
      sizeof(double);    //cached_constants_2b,
  bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * (max_num_coeff_2b - 1) * 9 *
      sizeof(double);    //cached_constants_2b_deri

  if (pot_3b) {
    bytes += (double) (num_of_elements + 1) * (num_of_elements + 1) * (num_of_elements + 1) *
        sizeof(int);    //map_3b

    bytes += (double) tot_interaction_count_3b * 3 * max_num_knots_3b *
        sizeof(double);    //n3b_knots_array
    bytes += (double) tot_interaction_count_3b * max_num_coeff_3b * max_num_coeff_3b *
        max_num_coeff_3b * sizeof(double);    //n3b_coeff_array

    bytes += (double) tot_interaction_count_3b * 3 * sizeof(int);    //n3b_knots_array_size
    bytes += (double) tot_interaction_count_3b * 3 * sizeof(int);    //n3b_coeff_array_size

    bytes += (double) tot_interaction_count_3b * max_num_coeff_3b * max_num_coeff_3b *
        max_num_coeff_3b * 3 *
        sizeof(double);    //coeff_for_der_jk coeff_for_der_ik coeff_for_der_ij

    bytes += (double) tot_interaction_count_3b * 3 * max_num_coeff_3b * 16 *
        sizeof(double);    //cached_constants_3b
    bytes += (double) tot_interaction_count_3b * 3 * (max_num_coeff_3b - 1) * 16 *
        sizeof(double);    //cached_constants_3b_deri
  }

  bytes += (double) maxshort * sizeof(int);    //neighshort

  bytes += (double) 6 * sizeof(int);     //maxshort, bsplines_created, nbody_flag,
                                         //max_num_knots_2b, max_num_coeff_2b,
                                         //max_num_knots_3b, max_num_coeff_3b
  bytes += (double) 1 * sizeof(bool);    //pot_3b

  return bytes;
}
