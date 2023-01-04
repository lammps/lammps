// unit tests for accessing neighbor lists in a LAMMPS instance through the Fortran wrapper

#include "lammps.h"
#include "library.h"
#include "force.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "info.h"
//#include <cstdint>
//#include <cstdlib>
#include <mpi.h>
#include <string>

#include "gtest/gtest.h"

// prototypes for Fortran reverse wrapper functions
extern "C" {
void *f_lammps_with_args();
void f_lammps_close();
void f_lammps_setup_neigh_tests();
int f_lammps_pair_neighlist_test();
int f_lammps_fix_neighlist_test();
int f_lammps_compute_neighlist_test();
int f_lammps_neighlist_num_elements(int);
}

namespace LAMMPS_NS {

class LAMMPS_neighbors : public ::testing::Test {
protected:
    LAMMPS_NS::LAMMPS *lmp;
    LAMMPS_neighbors()           = default;
    ~LAMMPS_neighbors() override = default;

    void SetUp() override {
        ::testing::internal::CaptureStdout();
        lmp                = (LAMMPS_NS::LAMMPS *)f_lammps_with_args();
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 8).c_str(), "LAMMPS (");
    }
    void TearDown() override {
        ::testing::internal::CaptureStdout();
        f_lammps_close();
        std::string output = ::testing::internal::GetCapturedStdout();
        EXPECT_STREQ(output.substr(0, 16).c_str(), "Total wall time:");
        lmp = nullptr;
    }
};

TEST_F(LAMMPS_neighbors, pair)
{
    f_lammps_setup_neigh_tests();
    int pair_neighlist = f_lammps_pair_neighlist_test();
    Pair *pair = lmp->force->pair_match("lj/cut",1,0);
    int index = -2;
    if (pair != nullptr) {
        for (int i = 0; i < lmp->neighbor->nlist; i++) {
            NeighList *list = lmp->neighbor->lists[i];
            if ((list->requestor_type == NeighList::PAIR)
                    and (pair == list->requestor)
                    and (list->id == 0)) {
                index = i;
                break;
            }
        }
    }
    EXPECT_EQ(index, pair_neighlist);
};

TEST_F(LAMMPS_neighbors, fix)
{
    if (not Info::has_package("REPLICA")) GTEST_SKIP();
    f_lammps_setup_neigh_tests();
    auto fix = lmp->modify->get_fix_by_id("f");
    EXPECT_NE(fix, nullptr);
    int ilist = -2;
    for (int i = 0; i < lmp->neighbor->nlist; i++) {
      NeighList *list = lmp->neighbor->lists[i];
      if ( (list->requestor_type == NeighList::FIX)
              and (fix == list->requestor) and (list->id == 0) ) {
          ilist = i;
          break;
      }
    }
    EXPECT_EQ(ilist, f_lammps_fix_neighlist_test());
};

TEST_F(LAMMPS_neighbors, compute)
{
    f_lammps_setup_neigh_tests();
    auto compute = lmp->modify->get_compute_by_id("c");
    EXPECT_NE(compute,nullptr);
    int ilist = -2;
    for (int i=0; i < lmp->neighbor->nlist; i++) {
        NeighList *list = lmp->neighbor->lists[i];
        if ( (list->requestor_type == NeighList::COMPUTE)
                and (compute == list->requestor) and (list->id == 0) ) {
            ilist = i;
            break;
        }
    }
    EXPECT_EQ(ilist, f_lammps_compute_neighlist_test());
};

TEST_F(LAMMPS_neighbors, numelements)
{
    f_lammps_setup_neigh_tests();
    int num_neigh = 0;
    int pair_id = f_lammps_pair_neighlist_test();
    num_neigh = f_lammps_neighlist_num_elements(pair_id);
    EXPECT_EQ(num_neigh, lammps_neighlist_num_elements(lmp, pair_id));
    if (Info::has_package("REPLICA")) {
      int fix_id = f_lammps_fix_neighlist_test();
      num_neigh = f_lammps_neighlist_num_elements(fix_id);
      EXPECT_EQ(num_neigh, lammps_neighlist_num_elements(lmp, fix_id));
    }
    int compute_id = f_lammps_compute_neighlist_test();
    num_neigh = f_lammps_neighlist_num_elements(compute_id);
    EXPECT_EQ(num_neigh, lammps_neighlist_num_elements(lmp, compute_id));
};

} // LAMMPS_NS
