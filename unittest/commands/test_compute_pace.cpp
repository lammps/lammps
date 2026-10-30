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

// unit tests for compute pace and compute pace/atom (ML-PACE package).
// the output-style reference test compute-pace_dgrad.yaml pins the
// bikflag=1 dgradflag=1 global array to stored numbers; the tests here
// check that the other global array layouts and the per-atom array are
// consistent with it, that the virial rows are the strain derivatives of
// the summed descriptors, and that the computes behave across MPI ranks,
// multiple runs and a changing number of atoms

#include "../testing/core.h"
#include "../testing/test_mpi_main.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "fmt/format.h"
#include "group.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <string>
#include <vector>

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

namespace LAMMPS_NS {

static const std::string yace = STRINGIFY(TEST_INPUT_FOLDER) "/pace_hno.yace";

class ComputePACETest : public LAMMPSTest {
protected:
    // the small periodic H/N/O cell of the output-style reference test,
    // with pair_style lj/cut for the reference energy, forces and virial
    void hno_cell()
    {
        HIDE_OUTPUT([&] {
            command("variable input_dir index " STRINGIFY(TEST_INPUT_FOLDER));
            command("include " STRINGIFY(TEST_INPUT_FOLDER) "/in.pace-hno");
            command("pair_style lj/cut 5.7");
            command("pair_coeff * * 0.1 1.0");
        });
    }

    // a Lennard-Jones fcc lattice filling the region "atoms" of the box
    void lj_lattice(const std::string &box, const std::string &atoms)
    {
        HIDE_OUTPUT([&] {
            command("units lj");
            command("atom_style atomic");
            command("atom_modify map array sort 0 0.0");
            command("lattice fcc 0.8442");
            command("region box block " + box);
            command("region atoms block " + atoms);
            command("create_box 1 box");
            command("create_atoms 1 region atoms");
            command("mass 1 1.0");
            command("displace_atoms all random 0.05 0.05 0.05 12345 units box");
            command("pair_style lj/cut 5.7");
            command("pair_coeff 1 1 1.0 1.0");
        });
    }

    // the TEST_FAILURE macro matches the screen output, which only MPI rank 0
    // prints, so the error message is matched in the exception instead
    void expect_error(const std::string &errmsg, const std::string &cmd)
    {
        std::string mesg;
        BEGIN_HIDE_OUTPUT();
        try {
            command(cmd);
        } catch (LAMMPSException &e) {
            mesg = e.what();
        }
        END_HIDE_OUTPUT();
        EXPECT_THAT(mesg, ContainsRegex(errmsg));
    }

    // the per-atom descriptors of the local atoms must equal the per-atom
    // rows of the global array (bikflag = 1) in the column block of their type
    void check_peratom(Compute *peratom, Compute *global)
    {
        const int nvalues = peratom->size_peratom_cols;
        const auto *tag   = lmp->atom->tag;
        const auto *type  = lmp->atom->type;
        for (int i = 0; i < lmp->atom->nlocal; ++i)
            for (int k = 0; k < nvalues; ++k)
                EXPECT_DOUBLE_EQ(peratom->array_atom[i][k],
                                 global->array[tag[i] - 1][(type[i] - 1) * nvalues + k])
                    << "atom " << tag[i] << " column " << k + 1;
    }
};

// the four views of the same descriptors and gradients must agree:
//   pace/atom   per-atom array of B_i
//   pace 1 0    per-atom rows of B_i, rows of sum_i -dB_i/dr_j, virial rows
//   pace 0 0    one row of sum_i B_i instead of the per-atom rows
//   pace 1 1    per-atom rows of B_i and one row per pair (i,j) of -dB_i/dr_j
// the columns of pace 1 0 and pace 0 0 are blocks of nvalues per atom type,
// the last column holds the reference energy, forces and virial of the pair
// style, which pace 1 1 stores in the first three columns and a final row

TEST_F(ComputePACETest, layouts)
{
    hno_cell();
    HIDE_OUTPUT([&] {
        command("compute atom all pace/atom " + yace);
        command("compute b1 all pace " + yace + " 1 0");
        command("compute b0 all pace " + yace + " 0 0");
        command("compute d1 all pace " + yace + " 1 1");
        command("compute vir all pressure NULL virial");
        command("run 0 post no");
    });
    auto *peratom = lmp->modify->get_compute_by_id("atom");
    auto *b1      = lmp->modify->get_compute_by_id("b1");
    auto *b0      = lmp->modify->get_compute_by_id("b0");
    auto *d1      = lmp->modify->get_compute_by_id("d1");
    auto *vir     = lmp->modify->get_compute_by_id("vir");
    auto *pe      = lmp->modify->get_compute_by_id("thermo_pe");
    ASSERT_NE(peratom, nullptr);
    ASSERT_NE(b1, nullptr);
    ASSERT_NE(b0, nullptr);
    ASSERT_NE(d1, nullptr);
    ASSERT_NE(vir, nullptr);
    ASSERT_NE(pe, nullptr);
    peratom->compute_peratom();
    b1->compute_array();
    b0->compute_array();
    d1->compute_array();
    vir->compute_vector();
    const double energy = pe->compute_scalar();

    const int natoms  = lmp->atom->natoms;
    const int ntypes  = lmp->atom->ntypes;
    const int nvalues = peratom->size_peratom_cols;
    const int ncols   = nvalues * ntypes;
    ASSERT_EQ(natoms, 6);
    ASSERT_EQ(ntypes, 3);
    ASSERT_EQ(b1->size_array_rows, natoms + 3 * natoms + 6);
    ASSERT_EQ(b1->size_array_cols, ncols + 1);
    ASSERT_EQ(b0->size_array_rows, 1 + 3 * natoms + 6);
    ASSERT_EQ(b0->size_array_cols, ncols + 1);
    ASSERT_EQ(d1->size_array_rows, natoms + 3 * natoms * natoms + 1);
    ASSERT_EQ(d1->size_array_cols, nvalues + 3);

    // atom types by atom ID, collected from all ranks

    std::vector<int> mytype(natoms + 1, 0), type(natoms + 1, 0);
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        mytype[lmp->atom->tag[i]] = lmp->atom->type[i];
    MPI_Allreduce(mytype.data(), type.data(), natoms + 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // rows of the three global arrays

    auto b1_grad = [&](int t, int d) {
        return b1->array[3 * (t - 1) + natoms + d];
    };
    auto b0_grad = [&](int t, int d) {
        return b0->array[3 * (t - 1) + 1 + d];
    };
    auto d1_grad = [&](int i, int j, int d) {
        return d1->array[natoms + 3 * natoms * (j - 1) + 3 * (i - 1) + d];
    };
    auto b1_virial = [&](int v) {
        return b1->array[4 * natoms + v];
    };
    auto b0_virial = [&](int v) {
        return b0->array[1 + 3 * natoms + v];
    };
    const int last = ncols;

    // descriptors: pace/atom == per-atom rows of pace 1 0 == pace 1 1,
    // the other type blocks of pace 1 0 are exactly zero

    check_peratom(peratom, b1);
    double maxval = 0.0;
    for (int t = 1; t <= natoms; ++t) {
        for (int k = 0; k < nvalues; ++k) {
            const double value = b1->array[t - 1][(type[t] - 1) * nvalues + k];
            EXPECT_DOUBLE_EQ(d1->array[t - 1][3 + k], value) << "atom " << t << " column " << k + 1;
            maxval = std::max(maxval, std::fabs(value));
        }
        for (int c = 0; c < ncols; ++c)
            if (c / nvalues != type[t] - 1) EXPECT_EQ(b1->array[t - 1][c], 0.0);
    }
    EXPECT_GT(maxval, 0.1);

    // summed descriptors: the row of pace 0 0 is the sum of the rows of pace 1 0

    for (int c = 0; c < ncols; ++c) {
        double sum = 0.0;
        for (int t = 1; t <= natoms; ++t)
            sum += b1->array[t - 1][c];
        EXPECT_NEAR(b0->array[0][c], sum, 1.0e-12 * (1.0 + std::fabs(sum))) << "column " << c + 1;
    }

    // gradients: the rows of pace 1 0 and pace 0 0 are identical, and the
    // sum over the atoms i of each type of the rows of pace 1 1 (which carry
    // the indices i, j, d in their first three columns)

    double maxgrad = 0.0;
    for (int j = 1; j <= natoms; ++j) {
        for (int d = 0; d < 3; ++d) {
            for (int c = 0; c < ncols; ++c) {
                const double value = b1_grad(j, d)[c];
                EXPECT_DOUBLE_EQ(b0_grad(j, d)[c], value) << "atom " << j << " column " << c + 1;
                double sum = 0.0;
                for (int i = 1; i <= natoms; ++i)
                    if (type[i] - 1 == c / nvalues) sum += d1_grad(i, j, d)[3 + c % nvalues];
                EXPECT_NEAR(value, sum, 1.0e-12 * (1.0 + std::fabs(sum)))
                    << "atom " << j << " column " << c + 1;
                maxgrad = std::max(maxgrad, std::fabs(value));
            }
            for (int i = 1; i <= natoms; ++i) {
                EXPECT_EQ(d1_grad(i, j, d)[0], i - 1);
                EXPECT_EQ(d1_grad(i, j, d)[1], j - 1);
                EXPECT_EQ(d1_grad(i, j, d)[2], d);
            }
        }
    }
    EXPECT_GT(maxgrad, 0.1);

    // translational invariance: the gradients of each descriptor sum to zero
    // over the atoms j, both per pair (pace 1 1) and per type (pace 1 0)

    for (int d = 0; d < 3; ++d) {
        for (int i = 1; i <= natoms; ++i) {
            for (int k = 0; k < nvalues; ++k) {
                double sum = 0.0;
                for (int j = 1; j <= natoms; ++j)
                    sum += d1_grad(i, j, d)[3 + k];
                EXPECT_NEAR(sum, 0.0, 1.0e-10) << "atom " << i << " column " << k + 1;
            }
        }
        for (int c = 0; c < ncols; ++c) {
            double sum = 0.0;
            for (int j = 1; j <= natoms; ++j)
                sum += b1_grad(j, d)[c];
            EXPECT_NEAR(sum, 0.0, 1.0e-10) << "column " << c + 1;
        }
    }

    // virial rows: identical between pace 1 0 and pace 0 0, and non-zero

    double maxvirial = 0.0;
    for (int v = 0; v < 6; ++v) {
        for (int c = 0; c < ncols; ++c) {
            EXPECT_DOUBLE_EQ(b0_virial(v)[c], b1_virial(v)[c])
                << "row " << v + 1 << " column " << c + 1;
            maxvirial = std::max(maxvirial, std::fabs(b1_virial(v)[c]));
        }
    }
    EXPECT_GT(maxvirial, 0.1);

    // reference energy, forces and virial of the pair style: identical in
    // the three layouts and equal to the thermo_pe, atom forces and
    // compute pressure values (the virial rows are in the order
    // xx, yy, zz, yz, xz, xy)

    EXPECT_NE(energy, 0.0);
    EXPECT_DOUBLE_EQ(b1->array[0][last], energy);
    EXPECT_DOUBLE_EQ(b0->array[0][last], energy);
    EXPECT_DOUBLE_EQ(d1->array[natoms + 3 * natoms * natoms][0], energy);
    for (int t = 2; t <= natoms; ++t)
        EXPECT_EQ(b1->array[t - 1][last], 0.0);

    double maxforce = 0.0;
    for (int j = 1; j <= natoms; ++j) {
        for (int d = 0; d < 3; ++d) {
            const double force = b1_grad(j, d)[last];
            EXPECT_DOUBLE_EQ(b0_grad(j, d)[last], force) << "atom " << j;
            EXPECT_DOUBLE_EQ(d1->array[j - 1][d], force) << "atom " << j;
            maxforce = std::max(maxforce, std::fabs(force));
        }
    }
    EXPECT_GT(maxforce, 0.1);
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        for (int d = 0; d < 3; ++d)
            EXPECT_DOUBLE_EQ(b1_grad(lmp->atom->tag[i], d)[last], lmp->atom->f[i][d])
                << "atom " << lmp->atom->tag[i];

    const int voigt[6] = {0, 1, 2, 5, 4, 3};
    for (int v = 0; v < 6; ++v) {
        EXPECT_DOUBLE_EQ(b1_virial(v)[last], vir->vector[voigt[v]]) << "row " << v + 1;
        EXPECT_DOUBLE_EQ(b0_virial(v)[last], vir->vector[voigt[v]]) << "row " << v + 1;
    }
}

// the virial rows of the global array hold -sum_j r_j^b dB/dr_j^a for the
// summed descriptors B of each type, i.e. minus the derivative of the
// summed descriptors with respect to the strain r^a += eps * r^b of the
// whole periodic cell. check them against central finite differences

TEST_F(ComputePACETest, virial_finite_difference)
{
    hno_cell();
    HIDE_OUTPUT([&] {
        command("compute sum all pace " + yace + " 0 0");
        command("change_box all triclinic");
        command("run 0 post no");
    });
    auto *sum = lmp->modify->get_compute_by_id("sum");
    ASSERT_NE(sum, nullptr);
    sum->compute_array();

    const int natoms = lmp->atom->natoms;
    const int ncols  = sum->size_array_cols - 1;
    std::vector<std::vector<double>> virial(6, std::vector<double>(ncols));
    for (int v = 0; v < 6; ++v)
        for (int c = 0; c < ncols; ++c)
            virial[v][c] = sum->array[1 + 3 * natoms + v][c];

    // summed descriptors after a change_box command

    auto summed = [&](const std::string &deform) {
        HIDE_OUTPUT([&] {
            command(deform);
            command("run 0 post no");
        });
        sum->compute_array();
        return std::vector<double>(sum->array[0], sum->array[0] + ncols);
    };

    // strain eps of the box about its center (xx, yy, zz) and tilt factors
    // eps * box length (yz, xz, xy). the descriptors are rotationally
    // invariant, so their virial tensor is symmetric and the tilt xy, which
    // shifts x in proportion to y, probes the xy row, which is accumulated
    // from the y derivative times x

    const double eps = 1.0e-5;
    double lo[3], hi[3];
    for (int d = 0; d < 3; ++d) {
        lo[d] = lmp->domain->boxlo[d];
        hi[d] = lmp->domain->boxhi[d];
    }
    const char *dims[3]  = {"x", "y", "z"};
    const char *tilts[3] = {"yz", "xz", "xy"};
    const int tiltdim[3] = {2, 2, 1};

    for (int v = 0; v < 6; ++v) {
        std::string plus, minus, reset;
        if (v < 3) {
            const double center = 0.5 * (lo[v] + hi[v]);
            const double half   = 0.5 * (hi[v] - lo[v]);
            auto box            = [&](double scale) {
                return fmt::format("change_box all {} final {:.17g} {:.17g} remap", dims[v],
                                   center - scale * half, center + scale * half);
            };
            plus  = box(1.0 + eps);
            minus = box(1.0 - eps);
            reset = box(1.0);
        } else {
            const double length = hi[tiltdim[v - 3]] - lo[tiltdim[v - 3]];
            auto box            = [&](double tilt) {
                return fmt::format("change_box all {} final {:.17g} remap", tilts[v - 3],
                                   tilt * length);
            };
            plus  = box(eps);
            minus = box(-eps);
            reset = box(0.0);
        }
        const auto bplus  = summed(plus);
        const auto bminus = summed(minus);
        HIDE_OUTPUT([&] {
            command(reset);
        });
        for (int c = 0; c < ncols; ++c) {
            const double deriv = (bplus[c] - bminus[c]) / (2.0 * eps);
            EXPECT_NEAR(-deriv, virial[v][c], 1.0e-7 * (1.0 + std::fabs(virial[v][c])))
                << "row " << v + 1 << " column " << c + 1;
        }
    }
}

// compute pace restricted to a group: the descriptor rows of atoms outside
// the group are zero, those of atoms in the group equal the unrestricted
// values, and the dB/dr rows sum the per-pair gradients of pace 1 1 over the
// group atoms only, since the neighbors are not filtered by the group.  the
// reference energy and forces are those of the whole system.  (the virial
// rows cannot be rebuilt from the per-pair gradients of a periodic cell.)

TEST_F(ComputePACETest, group)
{
    hno_cell();
    HIDE_OUTPUT([&] {
        command("group sub id 1 3 5");
        command("compute full all pace " + yace + " 1 1");
        command("compute part sub pace " + yace + " 1 0");
        command("run 0 post no");
    });
    auto *full = lmp->modify->get_compute_by_id("full");
    auto *part = lmp->modify->get_compute_by_id("part");
    ASSERT_NE(full, nullptr);
    ASSERT_NE(part, nullptr);
    full->compute_array();
    part->compute_array();

    const int natoms  = lmp->atom->natoms;
    const int ntypes  = lmp->atom->ntypes;
    const int nvalues = full->size_array_cols - 3;
    const int ncols   = nvalues * ntypes;
    ASSERT_EQ(part->size_array_cols, ncols + 1);

    // atom types and group membership by atom ID, collected from all ranks

    const int groupbit = lmp->group->get_bitmask_by_id(FLERR, "sub", "test");
    std::vector<int> mytype(natoms + 1, 0), type(natoms + 1, 0);
    std::vector<int> mygroup(natoms + 1, 0), ingroup(natoms + 1, 0);
    for (int i = 0; i < lmp->atom->nlocal; ++i) {
        mytype[lmp->atom->tag[i]]  = lmp->atom->type[i];
        mygroup[lmp->atom->tag[i]] = (lmp->atom->mask[i] & groupbit) ? 1 : 0;
    }
    MPI_Allreduce(mytype.data(), type.data(), natoms + 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(mygroup.data(), ingroup.data(), natoms + 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    int nsub = 0;
    for (int t = 1; t <= natoms; ++t)
        nsub += ingroup[t];
    ASSERT_EQ(nsub, 3);

    auto part_grad = [&](int j, int d) {
        return part->array[3 * (j - 1) + natoms + d];
    };
    auto full_grad = [&](int i, int j, int d) {
        return full->array[natoms + 3 * natoms * (j - 1) + 3 * (i - 1) + d];
    };

    for (int t = 1; t <= natoms; ++t) {
        for (int c = 0; c < ncols; ++c) {
            const double expected = (ingroup[t] && (c / nvalues == type[t] - 1))
                                        ? full->array[t - 1][3 + c % nvalues]
                                        : 0.0;
            EXPECT_DOUBLE_EQ(part->array[t - 1][c], expected)
                << "atom " << t << " column " << c + 1;
        }
    }

    double maxgrad = 0.0;
    for (int j = 1; j <= natoms; ++j) {
        for (int d = 0; d < 3; ++d) {
            for (int c = 0; c < ncols; ++c) {
                double sum = 0.0;
                for (int i = 1; i <= natoms; ++i)
                    if (ingroup[i] && (type[i] - 1 == c / nvalues))
                        sum += full_grad(i, j, d)[3 + c % nvalues];
                EXPECT_NEAR(part_grad(j, d)[c], sum, 1.0e-12 * (1.0 + std::fabs(sum)))
                    << "atom " << j << " column " << c + 1;
                maxgrad = std::max(maxgrad, std::fabs(sum));
            }
            EXPECT_DOUBLE_EQ(part_grad(j, d)[ncols], full->array[j - 1][d]) << "atom " << j;
        }
    }
    EXPECT_GT(maxgrad, 0.1);
    EXPECT_DOUBLE_EQ(part->array[0][ncols], full->array[natoms + 3 * natoms * natoms][0]);
}

// compute pace and pace/atom on a system in which all atoms sit in one
// corner of the box, away from the periodic boundaries so that the random
// displacement cannot wrap them around, so that with multiple MPI ranks
// all ranks but one own no atoms

TEST_F(ComputePACETest, empty_ranks)
{
    lj_lattice("0 4 0 4 0 4", "0.5 1.5 0.5 1.5 0.5 1.5");
    HIDE_OUTPUT([&] {
        command("compute glob all pace " + yace + " 1 0");
        command("compute test all pace/atom " + yace);
        command("run 0 post no");
    });
    auto *glob = lmp->modify->get_compute_by_id("glob");
    auto *test = lmp->modify->get_compute_by_id("test");
    ASSERT_NE(glob, nullptr);
    ASSERT_NE(test, nullptr);
    glob->compute_array();
    test->compute_peratom();
    check_peratom(test, glob);

    int nlocal = lmp->atom->nlocal;
    int empty  = (nlocal == 0) ? 1 : 0;
    int maxlocal, allempty, nprocs;
    MPI_Allreduce(&nlocal, &maxlocal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&empty, &allempty, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    EXPECT_GT(maxlocal, 0);
    EXPECT_EQ(allempty, nprocs - 1);
}

// compute pace and pace/atom across several run commands, each of which
// re-enters init(): a uniform translation of all atoms (across periodic
// boundaries and, with multiple MPI ranks, between ranks) must leave the
// global array unchanged, a few steps of dynamics must change it, and the
// values the computes hold after a run, invoked by the thermo output on its
// last step, must equal a fresh evaluation.  the per-atom descriptors must
// always equal the global descriptor rows.

TEST_F(ComputePACETest, multiple_runs)
{
    lj_lattice("0 4 0 4 0 4", "0 4 0 4 0 4");
    HIDE_OUTPUT([&] {
        command("fix 1 all nve");
        command("timestep 0.001");
        command("compute glob all pace " + yace + " 1 0");
        command("compute test all pace/atom " + yace);
        command("variable b1 equal C_test[1][1]");
        command("thermo_style custom step pe c_glob[1][1] v_b1");
        command("thermo 1");
    });
    auto *glob = lmp->modify->get_compute_by_id("glob");
    auto *test = lmp->modify->get_compute_by_id("test");
    ASSERT_NE(glob, nullptr);
    ASSERT_NE(test, nullptr);
    const int nrows = glob->size_array_rows;
    const int ncols = glob->size_array_cols;
    const int nvals = test->size_peratom_cols;

    // copies of the global array and of the per-atom array of the local atoms
    auto global_copy = [&]() {
        std::vector<double> copy;
        for (int r = 0; r < nrows; ++r)
            for (int c = 0; c < ncols; ++c)
                copy.push_back(glob->array[r][c]);
        return copy;
    };
    auto peratom_copy = [&]() {
        std::vector<double> copy;
        for (int i = 0; i < lmp->atom->nlocal; ++i)
            for (int m = 0; m < nvals; ++m)
                copy.push_back(test->array_atom[i][m]);
        return copy;
    };

    // evaluate both computes and check the per-atom data against the
    // descriptor rows of the global array; return a copy of the global array
    auto evaluate = [&]() {
        glob->compute_array();
        test->compute_peratom();
        check_peratom(test, glob);
        return global_copy();
    };
    auto maxdiff = [](const std::vector<double> &a, const std::vector<double> &b) {
        double diff = 0.0;
        for (std::size_t n = 0; n < a.size(); ++n)
            diff = std::max(diff, std::fabs(a[n] - b[n]) / std::max(1.0, std::fabs(a[n])));
        return diff;
    };

    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    const auto first = evaluate();

    HIDE_OUTPUT([&] {
        command("displace_atoms all move 0.37 -0.21 0.13 units box");
        command("run 0 post no");
    });
    const auto translated = evaluate();
    EXPECT_LT(maxdiff(first, translated), 1.0e-8);

    HIDE_OUTPUT([&] {
        command("run 5 post no");
    });
    EXPECT_EQ(glob->invoked_array, 5);
    EXPECT_EQ(test->invoked_peratom, 5);
    const auto during_global  = global_copy();
    const auto during_peratom = peratom_copy();
    const auto moved          = evaluate();
    EXPECT_EQ(maxdiff(during_global, moved), 0.0);
    EXPECT_EQ(maxdiff(during_peratom, peratom_copy()), 0.0);
    EXPECT_GT(maxdiff(first, moved), 1.0e-4);

    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    const auto repeated = evaluate();
    EXPECT_LT(maxdiff(moved, repeated), 1.0e-10);
}

// the global array of compute pace is sized for the atom count at definition
// and indexed by atom ID, so atoms created during a run must be an error
// instead of writing past the end of the array; compute pace/atom follows
// the local atom count and is unaffected

TEST_F(ComputePACETest, atoms_created_during_run)
{
    lj_lattice("0 4 0 4 0 6", "0 4 0 4 0 4");
    HIDE_OUTPUT([&] {
        command("region top block 0 4 0 4 4.8 5.8");
        command("fix 1 all nve");
        command("fix 2 all deposit 100 1 1 12345 region top near 0.5 attempt 100");
        command("timestep 0.001");
        command("compute test all pace/atom " + yace);
        command("variable b1 equal C_test[1][1]");
        command("thermo_style custom step atoms v_b1");
        command("thermo 1");
        command("run 2 post no");
    });
    const bigint natoms = lmp->atom->natoms;
    EXPECT_EQ(natoms, 290);
    auto *test = lmp->modify->get_compute_by_id("test");
    ASSERT_NE(test, nullptr);
    EXPECT_EQ(test->invoked_peratom, 2);

    HIDE_OUTPUT([&] {
        command("compute glob all pace " + yace + " 1 0");
        command("thermo_style custom step atoms c_glob[1][1] v_b1");
    });

    expect_error("Compute pace was defined for 290 atoms but there are now 291", "run 1 post no");
}

// the ACE evaluator rejects very small interatomic distances with a C++
// exception, which must be reported as a LAMMPS error.  it is raised with
// error->one() by the rank that owns the atoms, so this runs on one rank only

TEST_F(ComputePACETest, overlapping_atoms)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs > 1) GTEST_SKIP();
    HIDE_OUTPUT([&] {
        command("units metal");
        command("atom_style atomic");
        command("region box block 0 10 0 10 0 10");
        command("create_box 3 box");
        command("create_atoms 1 single 5.0 5.0 5.0");
        command("create_atoms 2 single 5.0 5.0 5.0000000001");
        command("mass * 1.0");
        command("pair_style zero 5.7");
        command("pair_coeff * *");
        command("compute glob all pace " + yace + " 1 0");
        command("thermo_style custom step c_glob[1][1]");
    });
    expect_error("ERROR on proc 0: Encountered very small distance", "run 0 post no");
}

} // namespace LAMMPS_NS
