// unit tests for the management of fixes and computes in the Modify class

#include "lammps.h"

#include "atom.h"
#include "atom_vec.h"
#include "compute.h"
#include "domain.h"
#include "fix.h"
#include "group.h"
#include "info.h"
#include "library.h"
#include "modify.h"
#include "output.h"
#include "platform.h"
#include "region.h"
#include "utils.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>

using ::testing::ContainsRegex;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using namespace FixConst;

class ModifyTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ModifyTest";
        LAMMPSTest::SetUp();
    }

    void InitSystem() override
    {
        HIDE_OUTPUT([&] {
            command("units lj");
            command("atom_style atomic");
            command("lattice fcc 0.8442");
            command("region box block 0 2 0 2 0 2");
            command("create_box 1 box");
            command("create_atoms 1 box");
            command("mass 1 1.0");
            command("velocity all create 1.0 12345 loop geom");
            command("pair_style lj/cut 2.5");
            command("pair_coeff 1 1 1.0 1.0 2.5");
            command("neighbor 0.3 bin");
            command("neigh_modify delay 0 every 1 check yes");
        });
    }

    // ID of the fix or compute at list position idx, empty string if none

    std::string fix_id(int idx)
    {
        auto *f = lmp->modify->get_fix_by_index(idx);
        return f ? f->id : "";
    }

    std::string compute_id(int idx)
    {
        auto *c = lmp->modify->get_compute_by_index(idx);
        return c ? c->id : "";
    }

    // ID at position idx as reported through the C library interface

    std::string library_id(const std::string &category, int idx)
    {
        char buf[64];
        if (lammps_id_name(lmp, category.c_str(), idx, buf, 64) == 0) return "";
        return buf;
    }

    // number of occurrences of needle in text

    static int count(const std::string &text, const std::string &needle)
    {
        int n = 0;
        for (auto pos = text.find(needle); pos != std::string::npos;
             pos      = text.find(needle, pos + 1))
            ++n;
        return n;
    }

    // order of marker lines in captured output as a string of letters

    static std::string marker_order(const std::string &text)
    {
        std::string order;
        for (const auto &line : utils::split_lines(text)) {
            if (line.find("AAAA") != std::string::npos) order += 'A';
            if (line.find("BBBB") != std::string::npos) order += 'B';
        }
        return order;
    }
};

TEST_F(ModifyTest, FixListOrderAndIndex)
{
    HIDE_OUTPUT([&] {
        command("fix f1 all nve");
        command("fix f2 all setforce 0.0 0.0 0.0");
        command("fix f3 all print 10 \"xxx\" screen no");
        command("fix f4 all addforce 1.0 0.0 0.0");
    });
    auto *modify      = lmp->modify;
    const auto &fixes = modify->get_fix_list();
    ASSERT_EQ(fixes.size(), 4);
    EXPECT_EQ(fix_id(0), "f1");
    EXPECT_EQ(fix_id(1), "f2");
    EXPECT_EQ(fix_id(2), "f3");
    EXPECT_EQ(fix_id(3), "f4");
    EXPECT_EQ(fix_id(4), "");
    EXPECT_EQ(fix_id(-1), "");
    EXPECT_EQ(lammps_id_count(lmp, "fix"), 4);
    EXPECT_EQ(library_id("fix", 0), "f1");
    EXPECT_EQ(library_id("fix", 3), "f4");
    EXPECT_EQ(library_id("fix", 4), "");

    // pointers from the list and from the lookup functions agree
    EXPECT_EQ(modify->get_fix_by_id("f3"), fixes[2]);
    EXPECT_EQ(modify->get_fix_by_id("f5"), nullptr);
    EXPECT_EQ(modify->get_fix_by_id(""), nullptr);
    EXPECT_EQ(modify->get_fix_by_style("^setforce").size(), 1);
    EXPECT_EQ(modify->get_fix_by_style("force").size(), 2);

    // masks reflect what the fix styles asked for
    EXPECT_TRUE(modify->get_fix_mask(fixes[0]) & INITIAL_INTEGRATE);
    EXPECT_TRUE(modify->get_fix_mask(fixes[0]) & FINAL_INTEGRATE);
    EXPECT_FALSE(modify->get_fix_mask(fixes[0]) & POST_FORCE);
    EXPECT_TRUE(modify->get_fix_mask(fixes[1]) & POST_FORCE);
    EXPECT_TRUE(modify->get_fix_mask(fixes[2]) & END_OF_STEP);
    EXPECT_FALSE(modify->get_fix_mask(fixes[2]) & POST_FORCE);
    EXPECT_TRUE(modify->get_fix_mask(fixes[3]) & POST_FORCE);
    EXPECT_EQ(modify->get_fix_mask(nullptr), 0);

    // counters are computed by init(), i.e. during a run
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    EXPECT_TRUE(modify->has_fixes(INITIAL_INTEGRATE));
    EXPECT_TRUE(modify->has_fixes(FINAL_INTEGRATE));
    EXPECT_TRUE(modify->has_fixes(POST_FORCE_ANY));
    EXPECT_TRUE(modify->has_fixes(END_OF_STEP));
    EXPECT_FALSE(modify->has_fixes(PRE_FORCE));
    EXPECT_FALSE(modify->has_fixes(ENERGY_GLOBAL));
    EXPECT_TRUE(modify->has_fixes(PRE_FORCE | POST_FORCE));
    EXPECT_FALSE(modify->has_fixes(PRE_FORCE | PRE_REVERSE | POST_FORCE_GROUP));
    EXPECT_FALSE(modify->has_fixes(0));
}

TEST_F(ModifyTest, PostForceInvocationOrder)
{
    // fix setforce zeroes forces, fix addforce adds to them: the result
    // depends on which runs first and fixes run in the order they are defined

    HIDE_OUTPUT([&] {
        command("fix f1 all nve");
        command("fix sf all setforce 0.0 0.0 0.0");
        command("fix af all addforce 1.0 0.0 0.0");
        command("run 0 post no");
    });
    EXPECT_DOUBLE_EQ(lmp->atom->f[0][0], 1.0);
    EXPECT_DOUBLE_EQ(lmp->atom->f[0][1], 0.0);

    HIDE_OUTPUT([&] {
        command("unfix sf");
        command("unfix af");
        command("fix af all addforce 1.0 0.0 0.0");
        command("fix sf all setforce 0.0 0.0 0.0");
        command("run 0 post no");
    });
    EXPECT_DOUBLE_EQ(lmp->atom->f[0][0], 0.0);

    // same through the minimizer path (min_setup and min_post_force)

    HIDE_OUTPUT([&] {
        command("minimize 0.0 0.0 1 1");
    });
    EXPECT_DOUBLE_EQ(lmp->atom->f[0][0], 0.0);

    HIDE_OUTPUT([&] {
        command("unfix sf");
        command("unfix af");
        command("fix sf all setforce 0.0 0.0 0.0");
        command("fix af all addforce 1.0 0.0 0.0");
        command("minimize 0.0 0.0 1 1");
    });
    EXPECT_DOUBLE_EQ(lmp->atom->f[0][0], 1.0);
}

TEST_F(ModifyTest, EndOfStepOrderAndNevery)
{
    // fix print prints in setup() and then every nevery steps
    HIDE_OUTPUT([&] {
        command("fix f1 all nve");
        command("fix pa all print 1 \"AAAA\"");
        command("fix pb all print 2 \"BBBB\"");
    });
    auto output = CAPTURE_OUTPUT([&] {
        command("run 4 post no");
    });
    EXPECT_EQ(count(output, "AAAA"), 5);
    EXPECT_EQ(count(output, "BBBB"), 3);
    EXPECT_EQ(marker_order(output), "ABAABAAB");

    // reversed definition order changes the invocation order
    HIDE_OUTPUT([&] {
        command("unfix pa");
        command("unfix pb");
        command("fix pb all print 2 \"BBBB\"");
        command("fix pa all print 1 \"AAAA\"");
    });
    output = CAPTURE_OUTPUT([&] {
        command("run 4 post no");
    });
    EXPECT_EQ(count(output, "AAAA"), 5);
    EXPECT_EQ(count(output, "BBBB"), 3);
    EXPECT_EQ(marker_order(output), "BAABAABA");
}

TEST_F(ModifyTest, UnfixCompactsAndKeepsOrder)
{
    HIDE_OUTPUT([&] {
        command("fix f1 all nve");
        command("fix sf all setforce 0.0 0.0 0.0");
        command("fix pr all print 10 \"xxx\" screen no");
        command("fix af all addforce 1.0 0.0 0.0");
        command("run 0 post no");
    });
    auto *modify = lmp->modify;
    EXPECT_TRUE(modify->has_fixes(END_OF_STEP));
    EXPECT_TRUE(modify->has_fixes(POST_FORCE_ANY));

    HIDE_OUTPUT([&] {
        command("unfix pr");
    });
    EXPECT_EQ(lammps_id_count(lmp, "fix"), 3);
    EXPECT_EQ(fix_id(0), "f1");
    EXPECT_EQ(fix_id(1), "sf");
    EXPECT_EQ(fix_id(2), "af");
    EXPECT_EQ(fix_id(3), "");
    EXPECT_EQ(modify->get_fix_by_id("pr"), nullptr);
    EXPECT_EQ(modify->get_fix_list().size(), 3);

    // order of the remaining fixes is unchanged and counters are recomputed
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    EXPECT_DOUBLE_EQ(lmp->atom->f[0][0], 1.0);
    EXPECT_FALSE(modify->has_fixes(END_OF_STEP));
    EXPECT_TRUE(modify->has_fixes(POST_FORCE_ANY));

    // remove first and last fix
    HIDE_OUTPUT([&] {
        command("unfix f1");
        command("unfix af");
    });
    EXPECT_EQ(lammps_id_count(lmp, "fix"), 1);
    EXPECT_EQ(fix_id(0), "sf");
    EXPECT_EQ(fix_id(1), "");

    // unfix of unknown ID is an error, list is unchanged
    TEST_FAILURE(".*Could not find fix ID xx to delete.*", command("unfix xx"););
    EXPECT_EQ(lammps_id_count(lmp, "fix"), 1);
    EXPECT_EQ(fix_id(0), "sf");
}

TEST_F(ModifyTest, ReplaceSameIdSameStyle)
{
    HIDE_OUTPUT([&] {
        command("fix f1 all nve");
        command("fix sf all setforce 0.0 0.0 0.0");
        command("fix af all addforce 1.0 0.0 0.0");
    });
    auto *modify = lmp->modify;

    // redefinition with the same ID and style replaces the fix in place
    HIDE_OUTPUT([&] {
        command("fix sf all setforce 0.0 1.0 0.0");
    });
    EXPECT_EQ(lammps_id_count(lmp, "fix"), 3);
    EXPECT_EQ(fix_id(1), "sf");
    EXPECT_TRUE(modify->get_fix_mask(modify->get_fix_by_id("sf")) & POST_FORCE);
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    EXPECT_DOUBLE_EQ(lmp->atom->f[0][0], 1.0);
    EXPECT_DOUBLE_EQ(lmp->atom->f[0][1], 1.0);

    // a different style for the same ID is an error, list is unchanged
    TEST_FAILURE(".*Replacing a fix, but new style != old style.*",
                 command("fix sf all addforce 1.0 0.0 0.0"););
    EXPECT_EQ(lammps_id_count(lmp, "fix"), 3);
    EXPECT_EQ(fix_id(1), "sf");
    EXPECT_EQ(utils::strip_style_suffix(modify->get_fix_by_index(1)->style, lmp), "setforce");
}

TEST_F(ModifyTest, FailedReplaceLeavesListValid)
{
#if defined(LAMMPS_SKIP_DEATH_TESTS)
    GTEST_SKIP() << "the failing command is not executed without death tests";
#endif
    HIDE_OUTPUT([&] {
        command("fix f1 all nve");
        command("fix sf all setforce 0.0 0.0 0.0");
        command("fix af all addforce 1.0 0.0 0.0");
    });
    auto *modify = lmp->modify;

    // the old fix is deleted before the new one is created, so when the
    // constructor of the new fix fails, the ID is gone and the list is
    // compacted; it must not contain an empty slot afterwards
    TEST_FAILURE(".*fix setforce.*", command("fix sf all setforce 0.0 0.0"););
    EXPECT_EQ(modify->get_fix_by_id("sf"), nullptr);
    EXPECT_EQ(lammps_id_count(lmp, "fix"), 2);
    EXPECT_EQ(fix_id(0), "f1");
    EXPECT_EQ(fix_id(1), "af");
    for (const auto &ifix : modify->get_fix_list())
        EXPECT_NE(ifix, nullptr);

    auto output = CAPTURE_OUTPUT([&] {
        command("info fixes");
    });
    EXPECT_THAT(output, ContainsRegex("Fix\\[ *1\\]: +af,"));
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
#if defined(LMP_KOKKOS_SINGLE_DOUBLE) || defined(LMP_KOKKOS_SINGLE_SINGLE)
    EXPECT_FLOAT_EQ(lmp->atom->f[0][0], 1.0);
#else
    EXPECT_DOUBLE_EQ(lmp->atom->f[0][0], 1.0);
#endif
}

TEST_F(ModifyTest, ReplaceFixKeepsSlot)
{
    if (!Info::has_package("GRANULAR")) GTEST_SKIP();

    // pair style granular creates a placeholder fix when the pair style is
    // defined and replaces it with its NEIGH_HISTORY fix during init().
    // the replacement must keep the position of the placeholder in the list.

    HIDE_OUTPUT([&] {
        command("clear");
        command("units lj");
        command("atom_style sphere");
        command("lattice fcc 0.8442");
        command("region box block 0 2 0 2 0 2");
        command("create_box 1 box");
        command("create_atoms 1 box");
        command("velocity all create 1.0 12345 loop geom");
        command("comm_modify vel yes");
        command("fix f1 all nve/sphere");
        command("pair_style granular");
        command("pair_coeff * * hooke 1000.0 50.0 tangential linear_history 500.0 1.0 0.4 "
                "damping mass_velocity");
        command("fix f2 all setforce 0.0 0.0 0.0");
    });
    auto *modify = lmp->modify;
    ASSERT_EQ(lammps_id_count(lmp, "fix"), 3);
    EXPECT_EQ(fix_id(0), "f1");
    EXPECT_EQ(std::string(modify->get_fix_by_index(1)->style), "DUMMY");
    EXPECT_EQ(fix_id(2), "f2");

    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    ASSERT_EQ(lammps_id_count(lmp, "fix"), 3);
    EXPECT_EQ(fix_id(0), "f1");
    EXPECT_EQ(std::string(modify->get_fix_by_index(1)->style), "NEIGH_HISTORY");
    EXPECT_THAT(fix_id(1), ContainsRegex("^NEIGH_HISTORY_GRANULAR"));
    EXPECT_EQ(fix_id(2), "f2");
    EXPECT_TRUE(modify->get_fix_mask(modify->get_fix_by_index(1)) & PRE_EXCHANGE);

    // changing the pair style deletes the history fix and closes the gap
    HIDE_OUTPUT([&] {
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0 2.5");
    });
    ASSERT_EQ(lammps_id_count(lmp, "fix"), 2);
    EXPECT_EQ(fix_id(0), "f1");
    EXPECT_EQ(fix_id(1), "f2");
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    EXPECT_DOUBLE_EQ(lmp->atom->f[0][0], 0.0);
}

TEST_F(ModifyTest, CallbackListsFollowFixOrder)
{
    // fixes with per-atom arrays register callbacks with the Atom class in
    // their constructor, i.e. before Modify has stored them in its list.
    // the callback lists must be in the order of that list, since per-atom
    // data is packed and unpacked by walking them, also after a fix was
    // re-defined or deleted.

    auto *atom   = lmp->atom;
    auto *modify = lmp->modify;
    HIDE_OUTPUT([&] {
        command("fix a all property/atom i_a ghost yes");
        command("fix b all nve");
        command("fix c all property/atom d_c");
    });
    Fix *fa = modify->get_fix_by_id("a");
    Fix *fc = modify->get_fix_by_id("c");
    ASSERT_NE(fa, nullptr);
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(atom->extra_grow, (std::vector<Fix *>{fa, fc}));
    EXPECT_EQ(atom->extra_restart, (std::vector<Fix *>{fa, fc}));
    EXPECT_EQ(atom->extra_border, (std::vector<Fix *>{fa}));

    // re-defining a fix keeps its position in both lists

    HIDE_OUTPUT([&] {
        command("fix a all property/atom i_a2 ghost yes");
    });
    fa = modify->get_fix_by_id("a");
    ASSERT_NE(fa, nullptr);
    EXPECT_EQ(fix_id(0), "a");
    EXPECT_EQ(atom->extra_grow, (std::vector<Fix *>{fa, fc}));
    EXPECT_EQ(atom->extra_restart, (std::vector<Fix *>{fa, fc}));
    EXPECT_EQ(atom->extra_border, (std::vector<Fix *>{fa}));

    // deleting a fix removes its callbacks without renumbering

    HIDE_OUTPUT([&] {
        command("unfix a");
    });
    EXPECT_EQ(atom->extra_grow, (std::vector<Fix *>{fc}));
    EXPECT_EQ(atom->extra_restart, (std::vector<Fix *>{fc}));
    EXPECT_TRUE(atom->extra_border.empty());

    HIDE_OUTPUT([&] {
        command("unfix c");
        command("run 0 post no");
    });
    EXPECT_TRUE(atom->extra_grow.empty());
    EXPECT_TRUE(atom->extra_restart.empty());
    EXPECT_TRUE(atom->extra_border.empty());
}

TEST_F(ModifyTest, CallbackListsAfterReplaceFix)
{
    if (!Info::has_package("GRANULAR")) GTEST_SKIP();

    // the NEIGH_HISTORY fix that replaces the placeholder of pair style
    // granular registers callbacks from its constructor. it must be sorted
    // into the position of the placeholder, between the two other fixes.

    HIDE_OUTPUT([&] {
        command("clear");
        command("units lj");
        command("atom_style sphere");
        command("lattice fcc 0.8442");
        command("region box block 0 2 0 2 0 2");
        command("create_box 1 box");
        command("create_atoms 1 box");
        command("velocity all create 1.0 12345 loop geom");
        command("comm_modify vel yes");
        command("fix f1 all property/atom i_f1");
        command("pair_style granular");
        command("pair_coeff * * hooke 1000.0 50.0 tangential linear_history 500.0 1.0 0.4 "
                "damping mass_velocity");
        command("fix f2 all property/atom d_f2");
        command("run 0 post no");
    });
    auto *atom   = lmp->atom;
    auto *modify = lmp->modify;
    ASSERT_EQ(lammps_id_count(lmp, "fix"), 3);
    Fix *f1 = modify->get_fix_by_index(0);
    Fix *fh = modify->get_fix_by_index(1);
    Fix *f2 = modify->get_fix_by_index(2);
    EXPECT_EQ(std::string(f1->id), "f1");
    EXPECT_EQ(std::string(fh->style), "NEIGH_HISTORY");
    EXPECT_EQ(std::string(f2->id), "f2");
    EXPECT_EQ(atom->extra_grow, (std::vector<Fix *>{f1, fh, f2}));
    EXPECT_EQ(atom->extra_restart, (std::vector<Fix *>{f1, fh, f2}));

    // the history fix is deleted with the pair style

    HIDE_OUTPUT([&] {
        command("pair_style lj/cut 2.5");
        command("pair_coeff * * 1.0 1.0 2.5");
        command("run 0 post no");
    });
    EXPECT_EQ(atom->extra_grow, (std::vector<Fix *>{f1, f2}));
    EXPECT_EQ(atom->extra_restart, (std::vector<Fix *>{f1, f2}));
}

TEST_F(ModifyTest, FailedFixConstructorRemovesCallbacks)
{
#if defined(LAMMPS_SKIP_DEATH_TESTS)
    GTEST_SKIP() << "the failing command is not executed without death tests";
#endif
    if (!Info::has_package("RIGID")) GTEST_SKIP();

    auto *atom   = lmp->atom;
    auto *modify = lmp->modify;
    HIDE_OUTPUT([&] {
        command("fix a all property/atom i_a");
    });
    Fix *fa = modify->get_fix_by_id("a");
    ASSERT_EQ(atom->extra_grow, (std::vector<Fix *>{fa}));

    // fix rigid registers its callbacks before it parses the body style, so
    // its constructor fails after the registration. the destructor of the
    // Fix base class must remove the callbacks of the half-built fix again.

    TEST_FAILURE(".*Unknown fix rigid mode.*", command("fix r all rigid bogus"););
    EXPECT_EQ(modify->get_fix_by_id("r"), nullptr);
    EXPECT_EQ(atom->extra_grow, (std::vector<Fix *>{fa}));
    EXPECT_EQ(atom->extra_restart, (std::vector<Fix *>{fa}));

    // growing the per-atom arrays must not reach the deleted fix

    const int nmax = atom->nmax;
    atom->avec->grow(0);
    EXPECT_GT(atom->nmax, nmax);
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
}

TEST_F(ModifyTest, DynamicMaskChange)
{
    if (!Info::has_package("RIGID")) GTEST_SKIP();

    HIDE_OUTPUT([&] {
        command("fix r all rigid single");
        command("run 0 post no");
    });
    auto *modify = lmp->modify;
    auto *rigid  = modify->get_fix_by_id("r");
    ASSERT_NE(rigid, nullptr);
    EXPECT_FALSE(modify->get_fix_mask(rigid) & POST_FORCE);
    EXPECT_FALSE(modify->has_fixes(POST_FORCE_ANY));

    // fix_modify changes the mask, counters follow at the next init()
    HIDE_OUTPUT([&] {
        command("fix_modify r bodyforces early");
    });
    EXPECT_TRUE(modify->get_fix_mask(rigid) & POST_FORCE);
    EXPECT_FALSE(modify->has_fixes(POST_FORCE_ANY));
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    EXPECT_TRUE(modify->has_fixes(POST_FORCE_ANY));

    HIDE_OUTPUT([&] {
        command("fix_modify r bodyforces late");
        command("run 0 post no");
    });
    EXPECT_FALSE(modify->get_fix_mask(rigid) & POST_FORCE);
    EXPECT_FALSE(modify->has_fixes(POST_FORCE_ANY));

    // set_fix_mask()/clear_fix_mask() only touch the requested bits
    EXPECT_TRUE(modify->get_fix_mask(rigid) & INITIAL_INTEGRATE);
    modify->clear_fix_mask(rigid, INITIAL_INTEGRATE);
    EXPECT_FALSE(modify->get_fix_mask(rigid) & INITIAL_INTEGRATE);
    EXPECT_TRUE(modify->get_fix_mask(rigid) & FINAL_INTEGRATE);
    modify->set_fix_mask(rigid, INITIAL_INTEGRATE);
    EXPECT_TRUE(modify->get_fix_mask(rigid) & INITIAL_INTEGRATE);
    modify->set_fix_mask(nullptr, POST_FORCE);
    modify->clear_fix_mask(nullptr, POST_FORCE);
}

TEST_F(ModifyTest, GroupFixRunsFirstInPostForce)
{
    // the fix GROUP updating a dynamic group is defined AFTER the fix using the
    // group, but Modify calls fix GROUP first in post_force() so that the group
    // is current when fix setforce is applied

    HIDE_OUTPUT([&] {
        command("region r block INF 1.0 INF INF INF INF");
        command("group dyn region r");
        command("fix f1 all nve");
        command("fix sf dyn setforce 0.0 0.0 0.0");
        command("group dyn dynamic all region r every 1");
        command("run 100 post no");
    });
    auto *modify = lmp->modify;
    EXPECT_EQ(fix_id(0), "f1");
    EXPECT_EQ(fix_id(1), "sf");
    EXPECT_EQ(std::string(modify->get_fix_by_index(2)->style), "GROUP");
    // fix GROUP is counted in n_post_force_any, although its mask is empty
    EXPECT_TRUE(modify->has_fixes(POST_FORCE_ANY));

    auto *region = lmp->domain->get_region_by_id("r");
    ASSERT_NE(region, nullptr);
    region->prematch();
    const int groupbit = lmp->group->bitmask[lmp->group->find("dyn")];
    const int nlocal   = lmp->atom->nlocal;
    double **x         = lmp->atom->x;
    double **f         = lmp->atom->f;
    int *mask          = lmp->atom->mask;
    int ninside = 0, noutside = 0;
    for (int i = 0; i < nlocal; ++i) {
        const bool inside = region->match(x[i][0], x[i][1], x[i][2]);
        EXPECT_EQ(inside, (mask[i] & groupbit) != 0);
        if (inside) {
            ++ninside;
            EXPECT_DOUBLE_EQ(f[i][0], 0.0);
            EXPECT_DOUBLE_EQ(f[i][1], 0.0);
            EXPECT_DOUBLE_EQ(f[i][2], 0.0);
        } else {
            ++noutside;
        }
    }
    EXPECT_GT(ninside, 0);
    EXPECT_GT(noutside, 0);
}

TEST_F(ModifyTest, RunWithoutSetupAfterChangeIsError)
{
    // "run ... pre no" skips the setup in which Modify takes the snapshot of
    // the fixes to invoke, so it is an error when fixes, computes, or dumps
    // were added, deleted, or modified since the last setup

    auto *modify = lmp->modify;
    auto *output = lmp->output;
    HIDE_OUTPUT([&] {
        command("fix f1 all nve");
        command("run 2 post no");
    });
    EXPECT_FALSE(modify->changed_since_init());
    HIDE_OUTPUT([&] {
        command("fix pa all print 1 \"AAAA\"");
    });
    EXPECT_TRUE(modify->changed_since_init());
    TEST_FAILURE(".*Cannot use run with 'pre no'.*", command("run 2 pre no post no"););

    // after a run with setup the fix is invoked, and a run without setup is
    // allowed when nothing changed in between
    auto out = CAPTURE_OUTPUT([&] {
        command("run 2 post no");
    });
    EXPECT_EQ(count(out, "AAAA"), 3);
    EXPECT_FALSE(modify->changed_since_init());
    out = CAPTURE_OUTPUT([&] {
        command("run 2 pre no post no");
    });
    EXPECT_EQ(count(out, "AAAA"), 2);
    EXPECT_FALSE(modify->changed_since_init());

    HIDE_OUTPUT([&] {
        command("fix_modify f1 dynamic/dof yes");
    });
    EXPECT_TRUE(modify->changed_since_init());
    TEST_FAILURE(".*Cannot use run with 'pre no'.*", command("run 2 pre no post no"););
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    EXPECT_FALSE(modify->changed_since_init());

    HIDE_OUTPUT([&] {
        command("unfix pa");
    });
    EXPECT_TRUE(modify->changed_since_init());
    TEST_FAILURE(".*Cannot use run with 'pre no'.*", command("run 2 pre no post no"););
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    EXPECT_FALSE(modify->changed_since_init());

    HIDE_OUTPUT([&] {
        command("compute ke all ke");
    });
    EXPECT_TRUE(modify->changed_since_init());
    TEST_FAILURE(".*Cannot use run with 'pre no'.*", command("run 2 pre no post no"););
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    EXPECT_FALSE(modify->changed_since_init());
    HIDE_OUTPUT([&] {
        command("uncompute ke");
    });
    EXPECT_TRUE(modify->changed_since_init());
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });

    // same for dumps
    EXPECT_FALSE(output->changed_since_init());
    HIDE_OUTPUT([&] {
        command("dump d1 all atom 100 test_modify_class.dump");
    });
    EXPECT_TRUE(output->changed_since_init());
    EXPECT_FALSE(modify->changed_since_init());
    TEST_FAILURE(".*Cannot use run with 'pre no'.*", command("run 2 pre no post no"););
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    EXPECT_FALSE(output->changed_since_init());
    HIDE_OUTPUT([&] {
        command("dump_modify d1 scale no");
    });
    EXPECT_TRUE(output->changed_since_init());
    TEST_FAILURE(".*Cannot use run with 'pre no'.*", command("run 2 pre no post no"););
    HIDE_OUTPUT([&] {
        command("run 0 post no");
        command("undump d1");
    });
    EXPECT_TRUE(output->changed_since_init());
    TEST_FAILURE(".*Cannot use run with 'pre no'.*", command("run 2 pre no post no"););
    HIDE_OUTPUT([&] {
        command("run 0 post no");
    });
    EXPECT_FALSE(output->changed_since_init());
    platform::unlink("test_modify_class.dump");
}

TEST_F(ModifyTest, ComputeListAndUncompute)
{
    // thermo_temp, thermo_press, thermo_pe are always the first three computes
    auto *modify = lmp->modify;
    EXPECT_EQ(lammps_id_count(lmp, "compute"), 3);
    EXPECT_EQ(compute_id(0), "thermo_temp");
    EXPECT_EQ(compute_id(1), "thermo_press");
    EXPECT_EQ(compute_id(2), "thermo_pe");

    HIDE_OUTPUT([&] {
        command("compute c1 all temp");
        command("compute c2 all pe");
        command("compute c3 all ke");
    });
    EXPECT_EQ(lammps_id_count(lmp, "compute"), 6);
    EXPECT_EQ(compute_id(3), "c1");
    EXPECT_EQ(library_id("compute", 5), "c3");
    EXPECT_EQ(library_id("compute", 6), "");
    EXPECT_EQ(modify->get_compute_by_id("c2"), modify->get_compute_by_index(4));
    EXPECT_EQ(modify->get_compute_by_style("^temp").size(), 2);

    // reuse of a compute ID is an error
    TEST_FAILURE(".*Reuse of compute ID 'c1'.*", command("compute c1 all ke"););
    EXPECT_EQ(lammps_id_count(lmp, "compute"), 6);

    HIDE_OUTPUT([&] {
        command("uncompute c2");
    });
    EXPECT_EQ(lammps_id_count(lmp, "compute"), 5);
    EXPECT_EQ(compute_id(3), "c1");
    EXPECT_EQ(compute_id(4), "c3");
    EXPECT_EQ(compute_id(5), "");
    EXPECT_EQ(modify->get_compute_by_id("c2"), nullptr);
    for (const auto &icompute : modify->get_compute_list())
        EXPECT_NE(icompute, nullptr);

    HIDE_OUTPUT([&] {
        command("fix f1 all nve");
        command("thermo_style custom step c_c1 c_c3");
        command("thermo 1");
    });
    auto output = CAPTURE_OUTPUT([&] {
        command("run 2 post no");
    });
    EXPECT_THAT(output, ContainsRegex("Step +c_c1 +c_c3"));
    TEST_FAILURE(".*Could not find compute ID c2 to delete.*", command("uncompute c2"););
}

TEST_F(ModifyTest, DeleteBetweenRuns)
{
    HIDE_OUTPUT([&] {
        command("fix f1 all nve");
        command("fix pa all print 1 \"AAAA\"");
    });
    auto output = CAPTURE_OUTPUT([&] {
        command("run 2 post no");
    });
    EXPECT_EQ(count(output, "AAAA"), 3);

    HIDE_OUTPUT([&] {
        command("unfix pa");
    });
    output = CAPTURE_OUTPUT([&] {
        command("run 2 post no");
    });
    EXPECT_EQ(count(output, "AAAA"), 0);
    EXPECT_FALSE(lmp->modify->has_fixes(END_OF_STEP));

    // fix_modify energy yes adds the fix energy to the potential energy
    HIDE_OUTPUT([&] {
        command("fix sp all spring tether 10.0 0.0 0.0 0.0 0.0");
        command("run 0 post no");
    });
    auto *modify = lmp->modify;
    auto *pe     = modify->get_compute_by_id("thermo_pe");
    EXPECT_FALSE(modify->has_fixes(ENERGY_GLOBAL));
    const double pe_without = pe->compute_scalar();
    HIDE_OUTPUT([&] {
        command("fix_modify sp energy yes");
        command("run 0 post no");
    });
    EXPECT_TRUE(modify->has_fixes(ENERGY_GLOBAL));
    const double spring = modify->get_fix_by_id("sp")->compute_scalar();
    EXPECT_GT(spring, 0.0);
    EXPECT_NEAR(pe->compute_scalar() - pe_without, spring, 1.0e-10);
    EXPECT_NEAR(modify->energy_global(), spring, 1.0e-10);
}

} // namespace LAMMPS_NS

using namespace LAMMPS_NS;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (platform::mpi_vendor() == "Open MPI" && !Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. "
                     "Death tests will be skipped\n";

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();

    // finalize the KOKKOS package explicitly: otherwise Kokkos is torn down by
    // static destructors at program exit, leading to segfaults in some cases

    lammps_kokkos_finalize();

    MPI_Finalize();
    return rv;
}
