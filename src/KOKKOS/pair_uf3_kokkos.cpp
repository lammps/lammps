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
#include "pair_uf3_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "kokkos_type.h"
#include "math_const.h"
#include "math_special_kokkos.h"
#include "memory.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_kokkos.h"
#include "text_file_reader.h"

#include <algorithm>
#include <cmath>
#include <utility>

using namespace LAMMPS_NS;
using namespace MathConst;
using MathSpecialKokkos::cube;
using MathSpecialKokkos::square;

template <class DeviceType> PairUF3Kokkos<DeviceType>::PairUF3Kokkos(LAMMPS *lmp) : PairUF3(lmp)
{
  respa_enable = 0;

  //kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TAG_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

}

template <class DeviceType> PairUF3Kokkos<DeviceType>::~PairUF3Kokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom, eatom); //destory eatom from host, set it to nullptr
                                              //Also set k_eatom to empty View
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
    destroy_3d(k_cut_3b,cut_3b);
    destroy_4d(k_min_cut_3b,min_cut_3b);
    eatom = NULL;
    vatom = NULL;
    cvatom = NULL;
  }
}

template <class DeviceType>
template <typename TYPE>
void PairUF3Kokkos<DeviceType>::destroy_3d(TYPE data, typename TYPE::value_type*** &array)
{
  if (array == nullptr) return;
  data = TYPE();
  memory->sfree(array);
  array = nullptr;
}

template <class DeviceType>
template <typename TYPE>
void PairUF3Kokkos<DeviceType>::destroy_4d(TYPE data, typename TYPE::value_type**** &array)
{
  if (array == nullptr) return;
  data = TYPE();
  memory->sfree(array);
  array = nullptr;
}

/* ----------------------------------------------------------------------
 *     global settings
 * ---------------------------------------------------------------------- */

template <class DeviceType> void PairUF3Kokkos<DeviceType>::settings(int narg, char **arg)
{
  PairUF3::settings(narg, arg);
  //1. Determines whether the simulation is 2-body or 2 and 3-body
  //2. Set nbody_flag, num_of_elements, pot_3b
}

/* ----------------------------------------------------------------------
 *    set coeffs for one or more type pairs
 * ---------------------------------------------------------------------- */
template <class DeviceType> void PairUF3Kokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairUF3::coeff(narg,arg);
  //Also calls allocate internally
  //Grows arrays to the right dimensions --> setflag, cutsq, cut, knot_spacing_type_2b,
  //knot_spacing_2b, n2b_knots_array_size, n2b_coeff_array_size, setflag_3b,
  //cut_3b, cut_3b_list, min_cut_3b, knot_spacing_type_3b, knot_spacing_3b,
  //tot_interaction_count_3b, map_3b, n3b_knots_array_size, n3b_coeff_array_size,
  //neighshort
  //
  //Also reads the pot_files_internally

}

template<class DeviceType>
void PairUF3Kokkos<DeviceType>::allocate()
{
  if (!allocated) PairUF3::allocate();
  int n = atom->ntypes;
  memory->destroy(cutsq); //Why are we destroying cutsq? cutsq is allocated when
  //PairUF3::coeff or PairUF3::allocate is called; in the next step when k_cutsq
  //is created cutsq is set to point to the host array of k_cutsq
  //memory->destroy(cut_3b);

  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>(); //assignment; get the device
  //view of k_cutsq and assign it to d_cutsq; in the header file we just
  //decleared d_cutsq's type
  memoryKK->create_kokkos(k_cut_3b,n+1,n+1,n+1,"threebody:cut");
  memoryKK->create_kokkos(k_min_cut_3b,n+1,n+1,n+1,3,"threebody:cut");
  d_cut_3b = k_cut_3b.template view<DeviceType>();
  d_min_cut_3b = k_min_cut_3b.template view<DeviceType>();
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template <class DeviceType> void PairUF3Kokkos<DeviceType>::init_style()
{

  PairUF3::init_style();

  neighflag = lmp->kokkos->neighflag;

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same<DeviceType, LMPHostType>::value &&
                           !std::is_same<DeviceType, LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType, LMPDeviceType>::value);

  request->enable_full();
  request->enable_ghost();
}

/* ----------------------------------------------------------------------
   init list sets the pointer to full neighbour list requested in previous function
------------------------------------------------------------------------- */

template <class DeviceType>
void PairUF3Kokkos<DeviceType>::init_list(int /*id*/, class NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
template <class DeviceType> double PairUF3Kokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairUF3::init_one(i, j);

  if (!coefficients_created) create_coefficients();

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone; //Update the k_cutsq's
  //host memory
  k_cutsq.template modify<LMPHostType>(); //Record that k_cutsq's host memory has
  //been modified

  return cutone;
}

template <class DeviceType> void PairUF3Kokkos<DeviceType>::create_coefficients()
{
  const int num_of_elements = atom->ntypes;
  coefficients_created = 1;

  if (pot_3b) {
    for (int i = 1; i < num_of_elements + 1; i++) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        for (int k = 1; k < num_of_elements + 1; k++) {
          k_cut_3b.h_view(i,j,k) = cut_3b[i][j][k];

          // Notice the order of min_cut_3b[i][j][k]
          //In min_cut_3b[i][j][k],
          //min_cut_3b[i][j][k][0] is the knot_vector along jk,
          //min_cut_3b[i][j][k][1] is the knot_vector along ik,
          //min_cut_3b[i][j][k][2] is the knot_vector along ij,
          //see pair_uf3.cpp for more details
          k_min_cut_3b.h_view(i,j,k,0) = min_cut_3b[i][j][k][0];
          k_min_cut_3b.h_view(i,j,k,1) = min_cut_3b[i][j][k][1];
          k_min_cut_3b.h_view(i,j,k,2) = min_cut_3b[i][j][k][2];
        }
      }
    }
    k_cut_3b.template modify<LMPHostType>();
    k_min_cut_3b.template modify<LMPHostType>();
  }

  //No allocation on device for --> setflag, cut, knot_spacing_type_2b,
  //n2b_knot, n2b_coeff, n2b_knot[i], n2b_coeff[i], setflag_3b, cut_3b,
  //cut_3b_list, min_cut_3b, knot_spacing_type_3b, cut_3b_list, n3b_knot_matrix,
  //neighshort

  //UFBS2b and UFBS3b are array of objects. Bad idea to use kokkos view(array)
  //for it
  create_2b_coefficients();
  if (pot_3b) create_3b_coefficients();

}

template <class DeviceType> void PairUF3Kokkos<DeviceType>::create_2b_coefficients()
{
  const int num_of_elements = atom->ntypes;

  // Setup interaction pair map
  //TODO: Instead of using map2b and map3b use simple indexing
  Kokkos::realloc(map2b, num_of_elements + 1, num_of_elements + 1);
  auto map2b_view = Kokkos::create_mirror(map2b);

  int interaction_count = 0;
  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = i; j < num_of_elements + 1; j++) {
      map2b_view(i, j) = interaction_count;
      map2b_view(j, i) = interaction_count++;
    }
  }
  Kokkos::deep_copy(map2b, map2b_view);

  // Count max knots for array size

  int max_knots = max_num_knots_2b;

  // Copy coefficients to view

  Kokkos::realloc(d_coefficients_2b, interaction_count, max_knots - 4);
  auto d_coefficients_2b_view = Kokkos::create_mirror(d_coefficients_2b);

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = i; j < num_of_elements + 1; j++) {
      for (int k = 0; k < max_num_coeff_2b; k++) {
        d_coefficients_2b_view(map2b_view(i, j), k) = n2b_coeff_array[i][j][k];
      }
    }
  }
  Kokkos::deep_copy(d_coefficients_2b, d_coefficients_2b_view);

  // Copy knots from array to view

  Kokkos::realloc(d_n2b_knot, interaction_count, max_knots);
  Kokkos::realloc(d_n2b_knot_spacings, interaction_count);
  auto d_n2b_knot_view = Kokkos::create_mirror(d_n2b_knot);
  auto d_n2b_knot_spacings_view = Kokkos::create_mirror(d_n2b_knot_spacings);

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = i; j < num_of_elements + 1; j++) {
      for (int k = 0; k < max_num_knots_2b; k++) {
        d_n2b_knot_view(map2b_view(i, j), k) = n2b_knots_array[i][j][k];
      }
      d_n2b_knot_spacings_view(map2b_view(i, j)) = n2b_knots_array[i][j][4] - n2b_knots_array[i][j][3];
    }
  }

  Kokkos::deep_copy(d_n2b_knot, d_n2b_knot_view);
  Kokkos::deep_copy(d_n2b_knot_spacings, d_n2b_knot_spacings_view);
  // Set spline constants

  Kokkos::realloc(constants_2b, interaction_count, max_knots - 4);
  auto constants_2b_view = Kokkos::create_mirror(constants_2b);

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = i; j < num_of_elements + 1; j++) {
      for (int l = 0; l < n2b_knots_array_size[i][j] - 4; l++) {
        //n2b_knot[i][j].size() - 4; l++) {
        auto c = get_constants(&n2b_knots_array[i][j][l], n2b_coeff_array[i][j][l]);
        for (int k = 0; k < 16; k++)
          constants_2b_view(map2b_view(i, j), l, k) = (std::isinf(c[k]) || std::isnan(c[k])) ? 0
                                                                                             : c[k];
      }
    }
  }
  Kokkos::deep_copy(constants_2b, constants_2b_view);

  Kokkos::realloc(dnconstants_2b, interaction_count, max_knots - 5);
  auto dnconstants_2b_view = Kokkos::create_mirror(dnconstants_2b);

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = i; j < num_of_elements + 1; j++) {
      for (int l = 0; l < n2b_knots_array_size[i][j] - 5; l++) {
        double dntemp4 = 3 / (n2b_knots_array[i][j][l + 4] - n2b_knots_array[i][j][l + 1]);
        double coeff = (n2b_coeff_array[i][j][l + 1] - n2b_coeff_array[i][j][l]) * dntemp4;
        auto c = get_dnconstants(&n2b_knots_array[i][j][l + 1], coeff);
        for (int k = 0; k < 9; k++)
          dnconstants_2b_view(map2b_view(i, j), l, k) =
              (std::isinf(c[k]) || std::isnan(c[k])) ? 0 : c[k];
      }
    }
  }
  Kokkos::deep_copy(dnconstants_2b, dnconstants_2b_view);
}

template <class DeviceType> void PairUF3Kokkos<DeviceType>::create_3b_coefficients()
{
  const int num_of_elements = atom->ntypes;
  // Init interaction map for 3B

  Kokkos::realloc(map3b, num_of_elements + 1, num_of_elements + 1, num_of_elements + 1);
  auto map3b_view = Kokkos::create_mirror(map3b);

  int interaction_count = 0;
  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      for (int k = 1; k < num_of_elements + 1; k++) {
        map3b_view(i, j, k) = interaction_count;
        interaction_count++;
      }
    }
  }
  Kokkos::deep_copy(map3b, map3b_view);

  // Count max knots for view

  int max_knots = max_num_knots_3b;
  //In n3b_knot_matrix[i][j][k],
  //n3b_knot_matrix[i][j][k][0] is the knot_vector along jk,
  //n3b_knot_matrix[i][j][k][1] is the knot_vector along ik,
  //n3b_knot_matrix[i][j][k][2] is the knot_vector along ij,
  //see pair_uf3.cpp for more details


  // Init knot matrix view

  Kokkos::realloc(d_n3b_knot_matrix, interaction_count, 3, max_knots);
  Kokkos::realloc(d_n3b_knot_matrix_spacings, interaction_count, 3);
  auto d_n3b_knot_matrix_view = Kokkos::create_mirror(d_n3b_knot_matrix);
  auto d_n3b_knot_matrix_spacings_view = Kokkos::create_mirror(d_n3b_knot_matrix_spacings);

  for (int i = 1; i < num_of_elements + 1; i++)
    for (int j = 1; j < num_of_elements + 1; j++)
      for (int k = 1; k < num_of_elements + 1; k++) {
        for (int m = 0; m < n3b_knots_array_size[map_3b[i][j][k]][0]; m++)
          d_n3b_knot_matrix_view(map3b_view(i, j, k), 0, m) =
              n3b_knots_array[map_3b[i][j][k]][0][m];
        for (int m = 0; m < n3b_knots_array_size[map_3b[i][j][k]][1]; m++)
          d_n3b_knot_matrix_view(map3b_view(i, j, k), 1, m) =
              n3b_knots_array[map_3b[i][j][k]][1][m];
        for (int m = 0; m < n3b_knots_array_size[map_3b[i][j][k]][2]; m++)
          d_n3b_knot_matrix_view(map3b_view(i, j, k), 2, m) =
              n3b_knots_array[map_3b[i][j][k]][2][m];

        d_n3b_knot_matrix_spacings_view(map3b_view(i, j, k),2) =
            n3b_knots_array[map_3b[i][j][k]][2][4] - n3b_knots_array[map_3b[i][j][k]][2][3];

        d_n3b_knot_matrix_spacings_view(map3b_view(i, j, k),1) =
            n3b_knots_array[map_3b[i][j][k]][1][4] - n3b_knots_array[map_3b[i][j][k]][1][3];

        d_n3b_knot_matrix_spacings_view(map3b_view(i, j, k),0) =
            n3b_knots_array[map_3b[i][j][k]][0][4] - n3b_knots_array[map_3b[i][j][k]][0][3];
      }
  Kokkos::deep_copy(d_n3b_knot_matrix, d_n3b_knot_matrix_view);
  Kokkos::deep_copy(d_n3b_knot_matrix_spacings, d_n3b_knot_matrix_spacings_view);

  // Set knots spacings

  Kokkos::realloc(d_n3b_knot_spacings, interaction_count, 3);
  auto d_n3b_knot_spacings_view = Kokkos::create_mirror(d_n3b_knot_spacings);

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      for (int k = 1; k < num_of_elements + 1; k++) {
        d_n3b_knot_spacings_view(map3b_view(i, j, k), 0) =
            1 / (n3b_knots_array[map_3b[i][j][k]][0][5] - n3b_knots_array[map_3b[i][j][k]][0][4]);

        d_n3b_knot_spacings_view(map3b_view(i, j, k), 1) =
            1 / (n3b_knots_array[map_3b[i][j][k]][1][5] - n3b_knots_array[map_3b[i][j][k]][1][4]);

        d_n3b_knot_spacings_view(map3b_view(i, j, k), 2) =
            1 / (n3b_knots_array[map_3b[i][j][k]][2][5] - n3b_knots_array[map_3b[i][j][k]][2][4]);
      }
    }
  }
  Kokkos::deep_copy(d_n3b_knot_spacings, d_n3b_knot_spacings_view);

  // Copy coefficients

  Kokkos::realloc(d_coefficients_3b, interaction_count, max_knots - 4, max_knots - 4,
                  max_knots - 4);
  auto d_coefficients_3b_view = Kokkos::create_mirror(d_coefficients_3b);

  for (int n = 1; n < num_of_elements + 1; n++) {
    for (int m = 1; m < num_of_elements + 1; m++) {
      for (int o = 1; o < num_of_elements + 1; o++) {
        for (int i = 0; i < n3b_coeff_array_size[map_3b[n][m][o]][0]; i++) {
          for (int j = 0; j < n3b_coeff_array_size[map_3b[n][m][o]][1]; j++) {
            for (int k = 0; k < n3b_coeff_array_size[map_3b[n][m][o]][2]; k++) {
              d_coefficients_3b_view(map3b_view(n, m, o), i, j, k) =
                  n3b_coeff_array[map_3b[n][m][o]][i][j][k];
            }
          }
        }
      }
    }
  }
  Kokkos::deep_copy(d_coefficients_3b, d_coefficients_3b_view);
  //
  // Create derivative coefficients

  // TODO: Shrink size
  Kokkos::realloc(d_dncoefficients_3b, interaction_count, 3, max_knots - 4, max_knots - 4,
                  max_knots - 4);
  auto d_dncoefficients_3b_view = Kokkos::create_mirror(d_dncoefficients_3b);

  //Notice the order for d_dncoefficients_3b_view(map3b_view(n, m, o), X, i, j, k)
  //d_dncoefficients_3b_view(map3b_view(n, m, o), 2, i, j, k) --> coeff for rjk
  //d_dncoefficients_3b_view(map3b_view(n, m, o), 1, i, j, k) --> coeff for rik
  //d_dncoefficients_3b_view(map3b_view(n, m, o), 0, i, j, k) --> coeff for rij
  //
  //This is because-
  //In n3b_knot_matrix[i][j][k],
  //n3b_knot_matrix[i][j][k][0] is the knot_vector along jk,
  //n3b_knot_matrix[i][j][k][1] is the knot_vector along ik,
  //n3b_knot_matrix[i][j][k][2] is the knot_vector along ij,
  //see pair_uf3.cpp for more details

  for (int n = 1; n < num_of_elements + 1; n++) {
    for (int m = 1; m < num_of_elements + 1; m++) {
      for (int o = 1; o < num_of_elements + 1; o++) {
        int coeff_dim1 = n3b_coeff_array_size[map_3b[n][m][o]][0];
        int coeff_dim2 = n3b_coeff_array_size[map_3b[n][m][o]][1];
        int coeff_dim3 = n3b_coeff_array_size[map_3b[n][m][o]][2];
        for (int i = 0; i < coeff_dim1; i++) {
          for (int j = 0; j < coeff_dim2; j++) {
            for (int k = 0; k < coeff_dim3; k++) {
              F_FLOAT dntemp4 =
                  3 / (n3b_knots_array[map_3b[n][m][o]][0][k + 4] - n3b_knots_array[map_3b[n][m][o]][0][k + 1]);

              d_dncoefficients_3b_view(map3b_view(n, m, o), 2, i, j, k) =
                  (n3b_coeff_array[map_3b[n][m][o]][i][j][k + 1] - n3b_coeff_array[map_3b[n][m][o]][i][j][k]) * dntemp4;
            }
          }
        }

        for (int i = 0; i < coeff_dim1; i++) {
          std::vector<std::vector<F_FLOAT>> dncoeff_vect2;
          for (int j = 0; j < coeff_dim2; j++) {
            F_FLOAT dntemp4 =
                3 / (n3b_knots_array[map_3b[n][m][o]][1][j + 4] - n3b_knots_array[map_3b[n][m][o]][1][j + 1]);

            std::vector<F_FLOAT> dncoeff_vect;
            for (int k = 0; k < coeff_dim3; k++) {
              d_dncoefficients_3b_view(map3b_view(n, m, o), 1, i, j, k) =
                  (n3b_coeff_array[map_3b[n][m][o]][i][j + 1][k] - n3b_coeff_array[map_3b[n][m][o]][i][j][k]) * dntemp4;
            }
          }
        }

        for (int i = 0; i < coeff_dim1; i++) {
          F_FLOAT dntemp4 =
              3 / (n3b_knots_array[map_3b[n][m][o]][2][i + 4] - n3b_knots_array[map_3b[n][m][o]][2][i + 1]);
          for (int j = 0; j < coeff_dim2; j++) {
            for (int k = 0; k < coeff_dim3; k++) {
              d_dncoefficients_3b_view(map3b_view(n, m, o), 0, i, j, k) =
                  (n3b_coeff_array[map_3b[n][m][o]][i + 1][j][k] - n3b_coeff_array[map_3b[n][m][o]][i][j][k]) * dntemp4;
            }
          }
        }
      }
    }
  }
  Kokkos::deep_copy(d_dncoefficients_3b, d_dncoefficients_3b_view);

  // Set spline constants

  Kokkos::realloc(constants_3b, interaction_count, 3, max_knots - 4);
  auto constants_3b_view = Kokkos::create_mirror(constants_3b);

  //In n3b_knot_matrix[i][j][k],
  //n3b_knot_matrix[i][j][k][0] is the knot_vector along jk,
  //n3b_knot_matrix[i][j][k][1] is the knot_vector along ik,
  //n3b_knot_matrix[i][j][k][2] is the knot_vector along ij,
  //see pair_uf3.cpp for more details
  for (int n = 1; n < num_of_elements + 1; n++) {
    for (int m = 1; m < num_of_elements + 1; m++) {
      for (int o = 1; o < num_of_elements + 1; o++) {
        for (int l = 0; l < n3b_knots_array_size[map_3b[n][m][o]][2] - 4; l++) {
          auto c = get_constants(&n3b_knots_array[map_3b[n][m][o]][2][l], 1);
          for (int k = 0; k < 16; k++)
            constants_3b_view(map3b_view(n, m, o), 0, l, k) =
                (std::isinf(c[k]) || std::isnan(c[k])) ? 0 : c[k];
        }
        for (int l = 0; l < n3b_knots_array_size[map_3b[n][m][o]][1] - 4; l++) {
          auto c = get_constants(&n3b_knots_array[map_3b[n][m][o]][1][l], 1);
          for (int k = 0; k < 16; k++)
            constants_3b_view(map3b_view(n, m, o), 1, l, k) =
                (std::isinf(c[k]) || std::isnan(c[k])) ? 0 : c[k];
        }
        for (int l = 0; l < n3b_knots_array_size[map_3b[n][m][o]][0] -4; l++) {
          auto c = get_constants(&n3b_knots_array[map_3b[n][m][o]][0][l], 1);
          for (int k = 0; k < 16; k++)
            constants_3b_view(map3b_view(n, m, o), 2, l, k) =
                (std::isinf(c[k]) || std::isnan(c[k])) ? 0 : c[k];
        }
      }
    }
  }
  Kokkos::deep_copy(constants_3b, constants_3b_view);

  Kokkos::realloc(dnconstants_3b, interaction_count, 3, max_knots - 6);
  auto dnconstants_3b_view = Kokkos::create_mirror(dnconstants_3b);

  for (int n = 1; n < num_of_elements + 1; n++) {
    for (int m = 1; m < num_of_elements + 1; m++) {
      for (int o = 1; o < num_of_elements + 1; o++) {
        for (int l = 1; l < n3b_knots_array_size[map_3b[n][m][o]][2] - 5; l++) {
          auto c = get_dnconstants(&n3b_knots_array[map_3b[n][m][o]][2][l], 1);
          for (int k = 0; k < 9; k++)
            dnconstants_3b_view(map3b_view(n, m, o), 0, l - 1, k) =
                (std::isinf(c[k]) || std::isnan(c[k])) ? 0 : c[k];
        }
        for (int l = 1; l < n3b_knots_array_size[map_3b[n][m][o]][1] - 5; l++) {
          auto c = get_dnconstants(&n3b_knots_array[map_3b[n][m][o]][1][l], 1);
          for (int k = 0; k < 9; k++)
            dnconstants_3b_view(map3b_view(n, m, o), 1, l - 1, k) =
                (std::isinf(c[k]) || std::isnan(c[k])) ? 0 : c[k];
        }
        for (int l = 1; l < n3b_knots_array_size[map_3b[n][m][o]][0] - 5; l++) {
          auto c = get_dnconstants(&n3b_knots_array[map_3b[n][m][o]][0][l], 1);
          for (int k = 0; k < 9; k++)
            dnconstants_3b_view(map3b_view(n, m, o), 2, l - 1, k) =
                (std::isinf(c[k]) || std::isnan(c[k])) ? 0 : c[k];
        }
      }
    }
  }
  Kokkos::deep_copy(dnconstants_3b, dnconstants_3b_view);
}

template <class DeviceType>
template <int EVFLAG>
KOKKOS_INLINE_FUNCTION void PairUF3Kokkos<DeviceType>::twobody(const int itype, const int jtype,
                                                               const F_FLOAT r, F_FLOAT &evdwl,
                                                               F_FLOAT &fpair) const
{

  // Find knot starting position
  int interaction_id = map2b(itype, jtype);
  int start_index = 3;
  while (r > d_n2b_knot(interaction_id, start_index + 1)) start_index++;
  //int start_index = 3+(int)((r-d_n2b_knot(interaction_id,0))/d_n2b_knot_spacings(interaction_id));

  F_FLOAT r_values[4];
  r_values[0] = 1;
  r_values[1] = r;
  r_values[2] = r_values[1] * r_values[1];

  if (EVFLAG) {
    r_values[3] = r_values[2] * r_values[1];
    // Calculate energy
    evdwl = constants_2b(interaction_id, start_index, 0);
    evdwl += r_values[1] * constants_2b(interaction_id, start_index, 1);
    evdwl += r_values[2] * constants_2b(interaction_id, start_index, 2);
    evdwl += r_values[3] * constants_2b(interaction_id, start_index, 3);
    evdwl += constants_2b(interaction_id, start_index - 1, 4);
    evdwl += r_values[1] * constants_2b(interaction_id, start_index - 1, 5);
    evdwl += r_values[2] * constants_2b(interaction_id, start_index - 1, 6);
    evdwl += r_values[3] * constants_2b(interaction_id, start_index - 1, 7);
    evdwl += constants_2b(interaction_id, start_index - 2, 8);
    evdwl += r_values[1] * constants_2b(interaction_id, start_index - 2, 9);
    evdwl += r_values[2] * constants_2b(interaction_id, start_index - 2, 10);
    evdwl += r_values[3] * constants_2b(interaction_id, start_index - 2, 11);
    evdwl += constants_2b(interaction_id, start_index - 3, 12);
    evdwl += r_values[1] * constants_2b(interaction_id, start_index - 3, 13);
    evdwl += r_values[2] * constants_2b(interaction_id, start_index - 3, 14);
    evdwl += r_values[3] * constants_2b(interaction_id, start_index - 3, 15);
  }

  // Calculate force
  fpair = dnconstants_2b(interaction_id, start_index - 1, 0);
  fpair += r_values[1] * dnconstants_2b(interaction_id, start_index - 1, 1);
  fpair += r_values[2] * dnconstants_2b(interaction_id, start_index - 1, 2);
  fpair += dnconstants_2b(interaction_id, start_index - 2, 3);
  fpair += r_values[1] * dnconstants_2b(interaction_id, start_index - 2, 4);
  fpair += r_values[2] * dnconstants_2b(interaction_id, start_index - 2, 5);
  fpair += dnconstants_2b(interaction_id, start_index - 3, 6);
  fpair += r_values[1] * dnconstants_2b(interaction_id, start_index - 3, 7);
  fpair += r_values[2] * dnconstants_2b(interaction_id, start_index - 3, 8);
}

template <class DeviceType>
template <int EVFLAG>
KOKKOS_INLINE_FUNCTION void PairUF3Kokkos<DeviceType>::threebody(
    const int itype, const int jtype, const int ktype, const F_FLOAT value_rij,
    const F_FLOAT value_rik, const F_FLOAT value_rjk, F_FLOAT &evdwl, F_FLOAT (&fforce)[3]) const
{
  evdwl = 0;
  fforce[0] = 0;
  fforce[1] = 0;
  fforce[2] = 0;

  F_FLOAT evals[3][4];
  F_FLOAT dnevals[3][4];
  int start_indices[3];
  F_FLOAT r[3] = {value_rij, value_rik, value_rjk};
  int interaction_id = map3b(itype, jtype, ktype);

  auto coefficients =
      Kokkos::subview(d_coefficients_3b, interaction_id, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
  auto dncoefficients = Kokkos::subview(d_dncoefficients_3b, interaction_id, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
  //Notice the 2-d in d_n3b_knot_matrix
  //
  //In d_n3b_knot_matrix[i][j][k],
  //d_n3b_knot_matrix[i][j][k][0] is the knot_vector along jk,
  //d_n3b_knot_matrix[i][j][k][1] is the knot_vector along ik,
  //d_n3b_knot_matrix[i][j][k][2] is the knot_vector along ij,
  //
  //and r[0] = rij, r[1] = rik and r[2] = rjk
  //see n3b_knot_matrix and pair_uf3.cpp for more details
  for (int d = 0; d < 3; d++) {
    start_indices[d] = 3;
    while (r[d] > d_n3b_knot_matrix(interaction_id, 2-d, start_indices[d] + 1)) start_indices[d]++;
    //start_indices[d] = 3+(int)((r[d]-d_n3b_knot_matrix(interaction_id, 2-d, 0))/d_n3b_knot_matrix_spacings(interaction_id, 2-d));

    F_FLOAT r_values[4];
    r_values[0] = 1;
    r_values[1] = r[d];
    r_values[2] = r_values[1] * r_values[1];

    r_values[3] = r_values[2] * r_values[1];

    // Calculate energy
    evals[d][0] = constants_3b(interaction_id, d, start_indices[d], 0);
    evals[d][0] += r_values[1] * constants_3b(interaction_id, d, start_indices[d], 1);
    evals[d][0] += r_values[2] * constants_3b(interaction_id, d, start_indices[d], 2);
    evals[d][0] += r_values[3] * constants_3b(interaction_id, d, start_indices[d], 3);
    evals[d][1] = constants_3b(interaction_id, d, start_indices[d] - 1, 4);
    evals[d][1] += r_values[1] * constants_3b(interaction_id, d, start_indices[d] - 1, 5);
    evals[d][1] += r_values[2] * constants_3b(interaction_id, d, start_indices[d] - 1, 6);
    evals[d][1] += r_values[3] * constants_3b(interaction_id, d, start_indices[d] - 1, 7);
    evals[d][2] = constants_3b(interaction_id, d, start_indices[d] - 2, 8);
    evals[d][2] += r_values[1] * constants_3b(interaction_id, d, start_indices[d] - 2, 9);
    evals[d][2] += r_values[2] * constants_3b(interaction_id, d, start_indices[d] - 2, 10);
    evals[d][2] += r_values[3] * constants_3b(interaction_id, d, start_indices[d] - 2, 11);
    evals[d][3] = constants_3b(interaction_id, d, start_indices[d] - 3, 12);
    evals[d][3] += r_values[1] * constants_3b(interaction_id, d, start_indices[d] - 3, 13);
    evals[d][3] += r_values[2] * constants_3b(interaction_id, d, start_indices[d] - 3, 14);
    evals[d][3] += r_values[3] * constants_3b(interaction_id, d, start_indices[d] - 3, 15);

    dnevals[d][0] = dnconstants_3b(interaction_id, d, start_indices[d] - 1, 0);
    dnevals[d][0] += r_values[1] * dnconstants_3b(interaction_id, d, start_indices[d] - 1, 1);
    dnevals[d][0] += r_values[2] * dnconstants_3b(interaction_id, d, start_indices[d] - 1, 2);
    dnevals[d][1] = dnconstants_3b(interaction_id, d, start_indices[d] - 2, 3);
    dnevals[d][1] += r_values[1] * dnconstants_3b(interaction_id, d, start_indices[d] - 2, 4);
    dnevals[d][1] += r_values[2] * dnconstants_3b(interaction_id, d, start_indices[d] - 2, 5);
    dnevals[d][2] = dnconstants_3b(interaction_id, d, start_indices[d] - 3, 6);
    dnevals[d][2] += r_values[1] * dnconstants_3b(interaction_id, d, start_indices[d] - 3, 7);
    dnevals[d][2] += r_values[2] * dnconstants_3b(interaction_id, d, start_indices[d] - 3, 8);
    dnevals[d][3] = 0;
  }

  if (EVFLAG) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          evdwl += coefficients(start_indices[0] - i, start_indices[1] - j, start_indices[2] - k) *
              evals[0][i] * evals[1][j] * evals[2][k];
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        fforce[0] += dncoefficients(0, start_indices[0] - 3 + i, start_indices[1] - 3 + j,
                                    start_indices[2] - 3 + k) *
            dnevals[0][2 - i] * evals[1][3 - j] * evals[2][3 - k];
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 4; k++) {
        fforce[1] += dncoefficients(1, start_indices[0] - 3 + i, start_indices[1] - 3 + j,
                                    start_indices[2] - 3 + k) *
            evals[0][3 - i] * dnevals[1][2 - j] * evals[2][3 - k];
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 3; k++) {
        fforce[2] += dncoefficients(2, start_indices[0] - 3 + i, start_indices[1] - 3 + j,
                                    start_indices[2] - 3 + k) *
            evals[0][3 - i] * evals[1][3 - j] * dnevals[2][2 - k];
      }
    }
  }
}

template <class DeviceType> void PairUF3Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag, vflag, 0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->create_kokkos(k_eatom, eatom, maxeatom, "pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space, datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  tag = atomKK->k_tag.template view<DeviceType>();
  type = atomKK->k_type.template view<DeviceType>();
  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  nall = atom->nlocal + atom->nghost;
  k_cutsq.template sync<DeviceType>(); //Sync the device memory of k_cutsq with
  //the array from the host memory; this updates d_cutsq also
  k_cut_3b.template sync<DeviceType>();
  k_min_cut_3b.template sync<DeviceType>();

  inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType> *k_list = static_cast<NeighListKokkos<DeviceType> *>(list);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

  copymode = 1;

  escatter = ScatterEType(d_eatom);
  fscatter = ScatterFType(f);
  vscatter = ScatterVType(d_vatom);
  //cvscatter = ScatterCVType(d_cvatom);

  EV_FLOAT ev;
  EV_FLOAT ev_all;

  // build short neighbor list

  int max_neighs = d_neighbors.extent(1);

  if (((int)d_neighbors_short.extent(1) != max_neighs) ||
          ((int)d_neighbors_short.extent(0) != ignum)) {
    d_neighbors_short = Kokkos::View<int **, DeviceType>("UF3::neighbors_short", ignum, max_neighs);
  }
  if (d_numneigh_short.extent(0) != ignum)
    d_numneigh_short = Kokkos::View<int *, DeviceType>("UF3::numneighs_short", ignum);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagPairUF3ComputeShortNeigh>(0, ignum), *this);

  // loop over neighbor list of my atoms

  if (evflag) {
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<DeviceType, TagPairUF3ComputeFullA<FULL, 1>>(0, inum), *this, ev);
  }
  else{
    Kokkos::parallel_for(
        Kokkos::RangePolicy<DeviceType, TagPairUF3ComputeFullA<FULL, 0>>(0, inum), *this);
  }
  ev_all += ev;

  Kokkos::Experimental::contribute(d_eatom, escatter);
  Kokkos::Experimental::contribute(d_vatom, vscatter);
  //Kokkos::Experimental::contribute(d_cvatom, cvscatter);
  Kokkos::Experimental::contribute(f, fscatter);

  if (eflag_global) eng_vdwl += ev_all.evdwl;
  if (vflag_global) {
    virial[0] += ev_all.v[0];
    virial[1] += ev_all.v[1];
    virial[2] += ev_all.v[2];
    virial[3] += ev_all.v[3];
    virial[4] += ev_all.v[4];
    virial[5] += ev_all.v[5];
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (cvflag_atom) {
    //k_cvatom.template modify<DeviceType>();
    //k_cvatom.template sync<LMPHostType>();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION void PairUF3Kokkos<DeviceType>::operator()(TagPairUF3ComputeShortNeigh,
                                                                  const int &ii) const
{
  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i, 0);
  const X_FLOAT ytmp = x(i, 1);
  const X_FLOAT ztmp = x(i, 2);

  const int jnum = d_numneigh[i];
  int inside = 0;
  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i, jj);
    j &= NEIGHMASK;

    const X_FLOAT delx = xtmp - x(j, 0);
    const X_FLOAT dely = ytmp - x(j, 1);
    const X_FLOAT delz = ztmp - x(j, 2);
    const F_FLOAT rsq = delx * delx + dely * dely + delz * delz;

    const int itype = type[i];
    const int jtype = type[j];

    if (rsq <= d_cutsq(itype, jtype)) {
      d_neighbors_short(i, inside) = j;
      inside++;
    }
  }
  d_numneigh_short(i) = inside;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION void
PairUF3Kokkos<DeviceType>::operator()(TagPairUF3ComputeFullA<NEIGHFLAG, EVFLAG>, const int &ii,
                                      EV_FLOAT &ev) const
{
  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = vscatter.access();
  auto a_f = fscatter.access();
  auto a_cvatom = cvscatter.access();

  F_FLOAT del_rji[3], del_rki[3], del_rkj[3], triangle_eval[3];
  F_FLOAT fij[3], fik[3], fjk[3];
  F_FLOAT fji[3], fki[3], fkj[3];
  F_FLOAT Fj[3], Fk[3];
  F_FLOAT evdwl = 0, evdwl3 = 0;
  F_FLOAT fpair = 0;

  const int i = d_ilist[ii];
  const int itype = type[i];
  const X_FLOAT xtmp = x(i, 0);
  const X_FLOAT ytmp = x(i, 1);
  const X_FLOAT ztmp = x(i, 2);

  // two-body interactions

  const int jnum = d_numneigh_short[i];

  F_FLOAT fxtmpi = 0.0;
  F_FLOAT fytmpi = 0.0;
  F_FLOAT fztmpi = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short(i, jj);
    j &= NEIGHMASK;
    const int jtype = type[j];

    const X_FLOAT delx = xtmp - x(j, 0);
    const X_FLOAT dely = ytmp - x(j, 1);
    const X_FLOAT delz = ztmp - x(j, 2);
    const F_FLOAT rsq = delx * delx + dely * dely + delz * delz;

    if (rsq >= d_cutsq(itype, jtype)) continue;

    const F_FLOAT rij = sqrt(rsq);
    this->template twobody<EVFLAG>(itype, jtype, rij, evdwl, fpair);

    fpair = -fpair / rij;

    fxtmpi += delx * fpair;
    fytmpi += dely * fpair;
    fztmpi += delz * fpair;
    a_f(j, 0) -= delx * fpair;
    a_f(j, 1) -= dely * fpair;
    a_f(j, 2) -= delz * fpair;

    if (EVFLAG) {
      if (eflag) ev.evdwl += evdwl;
      if (vflag_either || eflag_atom)
        this->template ev_tally<NEIGHFLAG>(ev, i, j, evdwl, fpair, delx, dely, delz);
    }
  }

  // 3-body interaction
  // jth atom
  const int jnumm1 = jnum - 1;
  for (int jj = 0; jj < jnumm1; jj++) {
    int j = d_neighbors_short(i, jj);
    j &= NEIGHMASK;
    const int jtype = type[j];
    del_rji[0] = x(j, 0) - xtmp;
    del_rji[1] = x(j, 1) - ytmp;
    del_rji[2] = x(j, 2) - ztmp;
    F_FLOAT rij = sqrt(del_rji[0] * del_rji[0] + del_rji[1] * del_rji[1] + del_rji[2] * del_rji[2]);

    F_FLOAT fxtmpj = 0.0;
    F_FLOAT fytmpj = 0.0;
    F_FLOAT fztmpj = 0.0;

    for (int kk = jj + 1; kk < jnum; kk++) {
      int k = d_neighbors_short(i, kk);
      k &= NEIGHMASK;
      const int ktype = type[k];

      // Notice the order of d_min_cut_3b[i][j][k]
      //In d_min_cut_3b[i][j][k],
      //d_min_cut_3b[i][j][k][0] is the knot_vector along jk,
      //d_min_cut_3b[i][j][k][1] is the knot_vector along ik,
      //d_min_cut_3b[i][j][k][2] is the knot_vector along ij,
      //see pair_uf3.cpp for more details
      if (rij < d_min_cut_3b(itype, jtype, ktype, 2)) continue;
      if (rij > d_cut_3b(itype, jtype, ktype)) continue;

      del_rki[0] = x(k, 0) - xtmp;
      del_rki[1] = x(k, 1) - ytmp;
      del_rki[2] = x(k, 2) - ztmp;
      F_FLOAT rik =
          sqrt(del_rki[0] * del_rki[0] + del_rki[1] * del_rki[1] + del_rki[2] * del_rki[2]);

      if (rik < d_min_cut_3b(itype, jtype, ktype, 1)) continue;
      if (rik > d_cut_3b(itype, ktype, jtype)) continue;

      del_rkj[0] = x(k, 0) - x(j, 0);
      del_rkj[1] = x(k, 1) - x(j, 1);
      del_rkj[2] = x(k, 2) - x(j, 2);
      F_FLOAT rjk =
          sqrt(del_rkj[0] * del_rkj[0] + del_rkj[1] * del_rkj[1] + del_rkj[2] * del_rkj[2]);
      if (rjk < d_min_cut_3b(itype, jtype, ktype, 0)) continue;
      this->template threebody<EVFLAG>(itype, jtype, ktype, rij, rik, rjk, evdwl3, triangle_eval);

      fij[0] = *(triangle_eval + 0) * (del_rji[0] / rij);
      fji[0] = -fij[0];
      fik[0] = *(triangle_eval + 1) * (del_rki[0] / rik);
      fki[0] = -fik[0];
      fjk[0] = *(triangle_eval + 2) * (del_rkj[0] / rjk);
      fkj[0] = -fjk[0];

      fij[1] = *(triangle_eval + 0) * (del_rji[1] / rij);
      fji[1] = -fij[1];
      fik[1] = *(triangle_eval + 1) * (del_rki[1] / rik);
      fki[1] = -fik[1];
      fjk[1] = *(triangle_eval + 2) * (del_rkj[1] / rjk);
      fkj[1] = -fjk[1];

      fij[2] = *(triangle_eval + 0) * (del_rji[2] / rij);
      fji[2] = -fij[2];
      fik[2] = *(triangle_eval + 1) * (del_rki[2] / rik);
      fki[2] = -fik[2];
      fjk[2] = *(triangle_eval + 2) * (del_rkj[2] / rjk);
      fkj[2] = -fjk[2];

      Fj[0] = fji[0] + fjk[0];
      Fj[1] = fji[1] + fjk[1];
      Fj[2] = fji[2] + fjk[2];

      Fk[0] = fki[0] + fkj[0];
      Fk[1] = fki[1] + fkj[1];
      Fk[2] = fki[2] + fkj[2];

      fxtmpi += (fij[0] + fik[0]);
      fytmpi += (fij[1] + fik[1]);
      fztmpi += (fij[2] + fik[2]);
      fxtmpj += Fj[0];
      fytmpj += Fj[1];
      fztmpj += Fj[2];
      a_f(k, 0) += Fk[0];
      a_f(k, 1) += Fk[1];
      a_f(k, 2) += Fk[2];

      if (EVFLAG) {
        if (eflag) { ev.evdwl += evdwl3; }
        if (vflag_either || eflag_atom) {
          this->template ev_tally3<NEIGHFLAG>(ev, i, j, k, evdwl3, 0.0, Fj, Fk, del_rji, del_rki);
          if (cvflag_atom) {

            F_FLOAT ric[3];
            ric[0] = THIRD * (-del_rji[0] - del_rki[0]);
            ric[1] = THIRD * (-del_rji[1] - del_rki[1]);
            ric[2] = THIRD * (-del_rji[2] - del_rki[2]);
            a_cvatom(i, 0) += ric[0] * (-Fj[0] - Fk[0]);
            a_cvatom(i, 1) += ric[1] * (-Fj[1] - Fk[1]);
            a_cvatom(i, 2) += ric[2] * (-Fj[2] - Fk[2]);
            a_cvatom(i, 3) += ric[0] * (-Fj[1] - Fk[1]);
            a_cvatom(i, 4) += ric[0] * (-Fj[2] - Fk[2]);
            a_cvatom(i, 5) += ric[1] * (-Fj[2] - Fk[2]);
            a_cvatom(i, 6) += ric[1] * (-Fj[0] - Fk[0]);
            a_cvatom(i, 7) += ric[2] * (-Fj[0] - Fk[0]);
            a_cvatom(i, 8) += ric[2] * (-Fj[1] - Fk[1]);

            double rjc[3];
            rjc[0] = THIRD * (del_rji[0] - del_rkj[0]);
            rjc[1] = THIRD * (del_rji[1] - del_rkj[1]);
            rjc[2] = THIRD * (del_rji[2] - del_rkj[2]);

            a_cvatom(j, 0) += rjc[0] * Fj[0];
            a_cvatom(j, 1) += rjc[1] * Fj[1];
            a_cvatom(j, 2) += rjc[2] * Fj[2];
            a_cvatom(j, 3) += rjc[0] * Fj[1];
            a_cvatom(j, 4) += rjc[0] * Fj[2];
            a_cvatom(j, 5) += rjc[1] * Fj[2];
            a_cvatom(j, 6) += rjc[1] * Fj[0];
            a_cvatom(j, 7) += rjc[2] * Fj[0];
            a_cvatom(j, 8) += rjc[2] * Fj[1];

            double rkc[3];
            rkc[0] = THIRD * (del_rki[0] + del_rkj[0]);
            rkc[1] = THIRD * (del_rki[1] + del_rkj[1]);
            rkc[2] = THIRD * (del_rki[2] + del_rkj[2]);

            a_cvatom(k, 0) += rkc[0] * Fk[0];
            a_cvatom(k, 1) += rkc[1] * Fk[1];
            a_cvatom(k, 2) += rkc[2] * Fk[2];
            a_cvatom(k, 3) += rkc[0] * Fk[1];
            a_cvatom(k, 4) += rkc[0] * Fk[2];
            a_cvatom(k, 5) += rkc[1] * Fk[2];
            a_cvatom(k, 6) += rkc[1] * Fk[0];
            a_cvatom(k, 7) += rkc[2] * Fk[0];
            a_cvatom(k, 8) += rkc[2] * Fk[1];
          }
        }
      }
    }
    a_f(j, 0) += fxtmpj;
    a_f(j, 1) += fytmpj;
    a_f(j, 2) += fztmpj;
  }

  a_f(i, 0) += fxtmpi;
  a_f(i, 1) += fytmpi;
  a_f(i, 2) += fztmpi;
}

template <class DeviceType>
template <int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION void
PairUF3Kokkos<DeviceType>::operator()(TagPairUF3ComputeFullA<NEIGHFLAG, EVFLAG>,
                                      const int &ii) const
{
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG, EVFLAG>(TagPairUF3ComputeFullA<NEIGHFLAG, EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION void
PairUF3Kokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j, const F_FLOAT &epair,
                                    const F_FLOAT &fpair, const F_FLOAT &delx, const F_FLOAT &dely,
                                    const F_FLOAT &delz) const
{

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for CUDA,
  // and neither for Serial

  auto a_eatom = escatter.access();
  auto a_vatom = vscatter.access();
  auto a_cvatom = cvscatter.access();

  if (eflag_atom) {
    const E_FLOAT epairhalf = 0.5 * epair;
    a_eatom[i] += epairhalf;
    a_eatom[j] += epairhalf;
  }

  if (vflag_either) {
    const E_FLOAT v0 = delx * delx * fpair;
    const E_FLOAT v1 = dely * dely * fpair;
    const E_FLOAT v2 = delz * delz * fpair;
    const E_FLOAT v3 = delx * dely * fpair;
    const E_FLOAT v4 = delx * delz * fpair;
    const E_FLOAT v5 = dely * delz * fpair;

    if (vflag_global) {
      ev.v[0] += v0;
      ev.v[1] += v1;
      ev.v[2] += v2;
      ev.v[3] += v3;
      ev.v[4] += v4;
      ev.v[5] += v5;
    }

    if (vflag_atom) {
      a_vatom(i, 0) += 0.5 * v0;
      a_vatom(i, 1) += 0.5 * v1;
      a_vatom(i, 2) += 0.5 * v2;
      a_vatom(i, 3) += 0.5 * v3;
      a_vatom(i, 4) += 0.5 * v4;
      a_vatom(i, 5) += 0.5 * v5;

      a_vatom(j, 0) += 0.5 * v0;
      a_vatom(j, 1) += 0.5 * v1;
      a_vatom(j, 2) += 0.5 * v2;
      a_vatom(j, 3) += 0.5 * v3;
      a_vatom(j, 4) += 0.5 * v4;
      a_vatom(j, 5) += 0.5 * v5;
    }

    if (cvflag_atom) {
      a_cvatom(i, 0) += 0.5 * v0;
      a_cvatom(i, 1) += 0.5 * v1;
      a_cvatom(i, 2) += 0.5 * v2;
      a_cvatom(i, 3) += 0.5 * v3;
      a_cvatom(i, 4) += 0.5 * v4;
      a_cvatom(i, 5) += 0.5 * v5;
      a_cvatom(i, 6) += 0.5 * v3;
      a_cvatom(i, 7) += 0.5 * v4;
      a_cvatom(i, 8) += 0.5 * v5;
      a_cvatom(j, 0) += 0.5 * v0;
      a_cvatom(j, 1) += 0.5 * v1;
      a_cvatom(j, 2) += 0.5 * v2;
      a_cvatom(j, 3) += 0.5 * v3;
      a_cvatom(j, 4) += 0.5 * v4;
      a_cvatom(j, 5) += 0.5 * v5;
      a_cvatom(j, 6) += 0.5 * v3;
      a_cvatom(j, 7) += 0.5 * v4;
      a_cvatom(j, 8) += 0.5 * v5;
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW and hbond potentials, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

template <class DeviceType>
template <int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION void
PairUF3Kokkos<DeviceType>::ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
                                     const F_FLOAT &evdwl, const F_FLOAT &ecoul, F_FLOAT *fj,
                                     F_FLOAT *fk, F_FLOAT *drji, F_FLOAT *drki) const
{
  F_FLOAT epairthird, v[6];

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for CUDA,
  // and neither for Serial

  auto a_eatom = escatter.access();
  auto a_vatom = vscatter.access();

  if (eflag_atom) {
    epairthird = THIRD * (evdwl + ecoul);
    a_eatom[i] += epairthird;
    a_eatom[j] += epairthird;
    a_eatom[k] += epairthird;
  }

  if (vflag_either) {
    v[0] = drji[0] * fj[0] + drki[0] * fk[0];
    v[1] = drji[1] * fj[1] + drki[1] * fk[1];
    v[2] = drji[2] * fj[2] + drki[2] * fk[2];
    v[3] = drji[0] * fj[1] + drki[0] * fk[1];
    v[4] = drji[0] * fj[2] + drki[0] * fk[2];
    v[5] = drji[1] * fj[2] + drki[1] * fk[2];

    if (vflag_global) {
      ev.v[0] += v[0];
      ev.v[1] += v[1];
      ev.v[2] += v[2];
      ev.v[3] += v[3];
      ev.v[4] += v[4];
      ev.v[5] += v[5];
    }

    if (vflag_atom) {
      a_vatom(i, 0) += THIRD * v[0];
      a_vatom(i, 1) += THIRD * v[1];
      a_vatom(i, 2) += THIRD * v[2];
      a_vatom(i, 3) += THIRD * v[3];
      a_vatom(i, 4) += THIRD * v[4];
      a_vatom(i, 5) += THIRD * v[5];

      a_vatom(j, 0) += THIRD * v[0];
      a_vatom(j, 1) += THIRD * v[1];
      a_vatom(j, 2) += THIRD * v[2];
      a_vatom(j, 3) += THIRD * v[3];
      a_vatom(j, 4) += THIRD * v[4];
      a_vatom(j, 5) += THIRD * v[5];

      a_vatom(k, 0) += THIRD * v[0];
      a_vatom(k, 1) += THIRD * v[1];
      a_vatom(k, 2) += THIRD * v[2];
      a_vatom(k, 3) += THIRD * v[3];
      a_vatom(k, 4) += THIRD * v[4];
      a_vatom(k, 5) += THIRD * v[5];
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW and hbond potentials, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

template <class DeviceType>
template <typename T, typename V>
void PairUF3Kokkos<DeviceType>::copy_2d(V &d, T **h, int m, int n)
{
  Kokkos::View<T **> tmp("pair::tmp", m, n); //Create tmp view(array) on
  //device memory

  //auto h_view = Kokkos::create_mirror_view(tmp);
  auto h_view = Kokkos::create_mirror(tmp); //Create a mirror of the device
  //view(array) tmp, as deep_copy is only possible for mirror views

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      h_view(i, j) = h[i][j]; //fill mirror
    }
    //views with data from normal array 'h' which always lives on host memory
  }

  Kokkos::deep_copy(tmp, h_view); //Deepcopy data from h_view(host) to tmp(device)

  d = tmp;
}

template <class DeviceType>
template <typename T, typename V>
void PairUF3Kokkos<DeviceType>::copy_3d(V &d, T ***h, int m, int n, int o)
{
  Kokkos::View<T ***> tmp("pair::tmp", m, n, o); //Create tmp view(array) on
  //device memory

  //auto h_view = Kokkos::create_mirror_view(tmp); //create_mirror always copies
  //the data. create_mirror_view only copies data if the host cannot access the
  //data
  auto h_view = Kokkos::create_mirror(tmp); //Create a mirror of the device
  //view(array) tmp, as deep_copy is only possible for mirror views
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < o; k++) { h_view(i, j, k) = h[i][j][k]; } //fill mirror
      //views with data from normal array 'h' which always lives on host memory
    }
  }

  Kokkos::deep_copy(tmp, h_view); //Deepcopy data from h_view(host) to tmp(device)

  d = tmp;
}

template <class DeviceType>
std::vector<F_FLOAT> PairUF3Kokkos<DeviceType>::get_constants(double *knots, double coefficient)
{

  std::vector<F_FLOAT> constants(16);

  constants[0] = coefficient *
      (-cube(knots[0]) /
       (-cube(knots[0]) + square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
        square(knots[0]) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  constants[1] = coefficient *
      (3 * square(knots[0]) /
       (-cube(knots[0]) + square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
        square(knots[0]) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  constants[2] = coefficient *
      (-3 * knots[0] /
       (-cube(knots[0]) + square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
        square(knots[0]) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  constants[3] = coefficient *
      (1 /
       (-cube(knots[0]) + square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
        square(knots[0]) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  constants[4] = coefficient *
      (square(knots[1]) * knots[4] /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) +
       square(knots[0]) * knots[2] /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) +
       knots[0] * knots[1] * knots[3] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])));
  constants[5] = coefficient *
      (-square(knots[1]) /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) -
       2 * knots[1] * knots[4] /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) -
       square(knots[0]) /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) -
       2 * knots[0] * knots[2] /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) -
       knots[0] * knots[1] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])) -
       knots[0] * knots[3] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])) -
       knots[1] * knots[3] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])));
  constants[6] = coefficient *
      (2 * knots[1] /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) +
       knots[4] /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) +
       2 * knots[0] /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) +
       knots[2] /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) +
       knots[0] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])) +
       knots[1] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])) +
       knots[3] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])));
  constants[7] = coefficient *
      (-1 /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) -
       1 /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) -
       1 /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])));
  constants[8] = coefficient *
      (-knots[0] * square(knots[3]) /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) -
       knots[1] * knots[3] * knots[4] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) -
       knots[2] * square(knots[4]) /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])));
  constants[9] = coefficient *
      (2 * knots[0] * knots[3] /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) +
       square(knots[3]) /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) +
       knots[1] * knots[3] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) +
       knots[1] * knots[4] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) +
       knots[3] * knots[4] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) +
       2 * knots[2] * knots[4] /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])) +
       square(knots[4]) /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])));
  constants[10] = coefficient *
      (-knots[0] /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) -
       2 * knots[3] /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) -
       knots[1] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) -
       knots[3] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) -
       knots[4] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) -
       knots[2] /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])) -
       2 * knots[4] /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])));
  constants[11] = coefficient *
      (1 /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) +
       1 /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) +
       1 /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])));
  constants[12] = coefficient *
      (cube(knots[4]) /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * square(knots[4]) +
        knots[2] * knots[3] * knots[4] - knots[2] * square(knots[4]) - knots[3] * square(knots[4]) +
        cube(knots[4])));
  constants[13] = coefficient *
      (-3 * square(knots[4]) /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * square(knots[4]) +
        knots[2] * knots[3] * knots[4] - knots[2] * square(knots[4]) - knots[3] * square(knots[4]) +
        cube(knots[4])));
  constants[14] = coefficient *
      (3 * knots[4] /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * square(knots[4]) +
        knots[2] * knots[3] * knots[4] - knots[2] * square(knots[4]) - knots[3] * square(knots[4]) +
        cube(knots[4])));
  constants[15] = coefficient *
      (-1 /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * square(knots[4]) +
        knots[2] * knots[3] * knots[4] - knots[2] * square(knots[4]) - knots[3] * square(knots[4]) +
        cube(knots[4])));

  return constants;
}

template <class DeviceType>
std::vector<F_FLOAT> PairUF3Kokkos<DeviceType>::get_dnconstants(double *knots, double coefficient)
{
  std::vector<F_FLOAT> constants(9);

  constants[0] = coefficient *
      (square(knots[0]) /
       (square(knots[0]) - knots[0] * knots[1] - knots[0] * knots[2] + knots[1] * knots[2]));
  constants[1] = coefficient *
      (-2 * knots[0] /
       (square(knots[0]) - knots[0] * knots[1] - knots[0] * knots[2] + knots[1] * knots[2]));
  constants[2] = coefficient *
      (1 / (square(knots[0]) - knots[0] * knots[1] - knots[0] * knots[2] + knots[1] * knots[2]));
  constants[3] = coefficient *
      (-knots[1] * knots[3] /
           (square(knots[1]) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) -
       knots[0] * knots[2] /
           (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + square(knots[2])));
  constants[4] = coefficient *
      (knots[1] /
           (square(knots[1]) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) +
       knots[3] /
           (square(knots[1]) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) +
       knots[0] /
           (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + square(knots[2])) +
       knots[2] /
           (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + square(knots[2])));
  constants[5] = coefficient *
      (-1 / (square(knots[1]) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) -
       1 / (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + square(knots[2])));
  constants[6] = coefficient *
      (square(knots[3]) /
       (knots[1] * knots[2] - knots[1] * knots[3] - knots[2] * knots[3] + square(knots[3])));
  constants[7] = coefficient *
      (-2 * knots[3] /
       (knots[1] * knots[2] - knots[1] * knots[3] - knots[2] * knots[3] + square(knots[3])));
  constants[8] = coefficient *
      (1 / (knots[1] * knots[2] - knots[1] * knots[3] - knots[2] * knots[3] + square(knots[3])));

  return constants;
}

template <class DeviceType>
double PairUF3Kokkos<DeviceType>::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                                         double /*factor_coul*/, double factor_lj, double &fforce)
{
  double value = 0.0;
  double r = sqrt(rsq);
  int interaction_id = map2b(itype, jtype);
  int start_index = 3;
  while (r > d_n2b_knot(interaction_id, start_index + 1)) start_index++;

  if (r < d_cutsq(itype, jtype)) {
    F_FLOAT r_values[4];
    r_values[0] = 1;
    r_values[1] = r;
    r_values[2] = r_values[1] * r_values[1];
    r_values[3] = r_values[2] * r_values[1];

    // Calculate energy
    value = constants_2b(interaction_id, start_index, 0);
    value += r_values[1] * constants_2b(interaction_id, start_index, 1);
    value += r_values[2] * constants_2b(interaction_id, start_index, 2);
    value += r_values[3] * constants_2b(interaction_id, start_index, 3);
    value += constants_2b(interaction_id, start_index - 1, 4);
    value += r_values[1] * constants_2b(interaction_id, start_index - 1, 5);
    value += r_values[2] * constants_2b(interaction_id, start_index - 1, 6);
    value += r_values[3] * constants_2b(interaction_id, start_index - 1, 7);
    value += constants_2b(interaction_id, start_index - 2, 8);
    value += r_values[1] * constants_2b(interaction_id, start_index - 2, 9);
    value += r_values[2] * constants_2b(interaction_id, start_index - 2, 10);
    value += r_values[3] * constants_2b(interaction_id, start_index - 2, 11);
    value += constants_2b(interaction_id, start_index - 3, 12);
    value += r_values[1] * constants_2b(interaction_id, start_index - 3, 13);
    value += r_values[2] * constants_2b(interaction_id, start_index - 3, 14);
    value += r_values[3] * constants_2b(interaction_id, start_index - 3, 15);

    // Calculate force
    fforce = dnconstants_2b(interaction_id, start_index - 1, 0);
    fforce += r_values[1] * dnconstants_2b(interaction_id, start_index - 1, 1);
    fforce += r_values[2] * dnconstants_2b(interaction_id, start_index - 1, 2);
    fforce += dnconstants_2b(interaction_id, start_index - 2, 3);
    fforce += r_values[1] * dnconstants_2b(interaction_id, start_index - 2, 4);
    fforce += r_values[2] * dnconstants_2b(interaction_id, start_index - 2, 5);
    fforce += dnconstants_2b(interaction_id, start_index - 3, 6);
    fforce += r_values[1] * dnconstants_2b(interaction_id, start_index - 3, 7);
    fforce += r_values[2] * dnconstants_2b(interaction_id, start_index - 3, 8);
  }

  return factor_lj * value;
}

namespace LAMMPS_NS {
template class PairUF3Kokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairUF3Kokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS
