/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// unit tests for the thermo, thermo_style, and thermo_modify commands

#include "lammps.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "info.h"
#include "input.h"
#include "output.h"
#include "thermo.h"
#include "utils.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstring>
#include <string>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::ContainsRegex;
using ::testing::HasSubstr;
using ::testing::Not;

// gtest's ContainsRegex() uses different regular expression engines on different
// platforms (POSIX ERE vs. its own limited fallback on Windows) whose common
// subset is too small.  Content checks with patterns therefore use the bundled
// (and thus platform-independent) LAMMPS regex implementation instead.
#define ASSERT_MATCH(text, pattern) \
    ASSERT_TRUE(utils::strmatch(text, pattern)) << "no match for pattern: " << (pattern)

// small LJ system; only core functionality is used, so the test runs with
// any package selection.  LJ units => thermo_modify norm defaults to yes.

class ThermoTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "ThermoTest";
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
            command("velocity all create 1.0 4928459");
            command("pair_style lj/cut 2.5");
            command("pair_coeff 1 1 1.0 1.0");
            command("neighbor 0.3 bin");
            command("fix nve all nve");
        });
    }

    std::string run0()
    {
        return CAPTURE_OUTPUT([&] {
            command("run 0 post no");
        });
    }

    std::string run(int nsteps)
    {
        return CAPTURE_OUTPUT([&] {
            command(fmt::format("run {} post no", nsteps));
        });
    }

    Thermo *thermo() { return lmp->output->thermo; }
};

TEST_F(ThermoTest, Styles)
{
    // default style "one"
    auto output = run0();
    ASSERT_MATCH(output, "Step +Temp +E_pair +E_mol +TotEng +Press *\n");
    ASSERT_MATCH(output, "\n +0 +1 +-[0-9.]+ +0 +-[0-9.]+ +-?[0-9.]+ *\n");

    HIDE_OUTPUT([&] {
        command("thermo_style multi");
    });
    output = run0();
    ASSERT_MATCH(output, "-+ Step +0 -+ CPU = +[0-9.]+ \\(sec\\) -+");
    ASSERT_MATCH(output, "TotEng += +-[0-9.]+ +KinEng += +[0-9.]+ +Temp += +1.0000");
    ASSERT_MATCH(output, "PotEng += +-[0-9.]+ +E_bond += +0.0000 +E_angle += +0.0000");
    ASSERT_MATCH(output, "E_dihed += +0.0000 +E_impro += +0.0000 +E_vdwl += +-[0-9.]+");
    ASSERT_MATCH(output, "E_coul += +0.0000 +E_long += +0.0000 +Press += +-?[0-9.]+");

    HIDE_OUTPUT([&] {
        command("thermo_style yaml");
    });
    output = run(2);
    ASSERT_THAT(
        output,
        HasSubstr("---\nkeywords: ['Step', 'Temp', 'KinEng', 'PotEng', 'E_bond', 'E_angle', "
                  "'E_dihed', 'E_impro', 'E_vdwl', 'E_coul', 'E_long', 'Press', ]\ndata:\n"));
    ASSERT_MATCH(output, "\n  - \\[0, 1, [0-9.]+, -[0-9.]+, 0, 0, 0, 0, -[0-9.]+, 0, 0, -?[0-9.]+, \\]\n");
    ASSERT_MATCH(output, "\n  - \\[2, [0-9.]+, .*\\]\n\\.\\.\\.\n");

    // line keyword switches the format, not the content
    HIDE_OUTPUT([&] {
        command("thermo_style one");
        command("thermo_modify line multi");
    });
    output = run0();
    ASSERT_MATCH(output, "-+ Step +2 -+ CPU = +[0-9.]+ \\(sec\\) -+");
    ASSERT_MATCH(output, "Step += +2 +Temp += +[0-9.]+ +E_pair += +-[0-9.]+");
    ASSERT_MATCH(output, "E_mol += +0.0000 +TotEng += +-[0-9.]+ +Press += +-?[0-9.]+");
    ASSERT_THAT(output, Not(HasSubstr("E_bond   =")));

    HIDE_OUTPUT([&] {
        command("thermo_modify line yaml");
    });
    output = run0();
    ASSERT_THAT(output, HasSubstr("---\nkeywords: ['Step', 'Temp', 'E_pair', 'E_mol', 'TotEng', "
                                  "'Press', ]\ndata:\n  - [2, "));
    ASSERT_THAT(output, HasSubstr("]\n...\n"));

    HIDE_OUTPUT([&] {
        command("thermo_modify line one");
    });
    output = run0();
    ASSERT_MATCH(output, "Step +Temp +E_pair +E_mol +TotEng +Press *\n +2 ");

    TEST_FAILURE(".*ERROR: Unknown thermo_modify line argument: xxx.*",
                 command("thermo_modify line xxx"););
    TEST_FAILURE(".*ERROR: Unknown thermo style xxx.*", command("thermo_style xxx"););
    // the previous thermo settings survive a failed thermo_style command
    output = run0();
    ASSERT_MATCH(output, "Step +Temp +E_pair +E_mol +TotEng +Press *\n +2 ");
    TEST_FAILURE(".*ERROR: Illegal thermo_style command.*", command("thermo_style"););
    TEST_FAILURE(".*ERROR: Illegal thermo_modify command.*", command("thermo_modify"););
    TEST_FAILURE(".*ERROR: Unknown thermo_modify keyword: xxx.*",
                 command("thermo_modify xxx yes"););
}

TEST_F(ThermoTest, Custom)
{
    HIDE_OUTPUT([&] {
        command("compute ke all ke");
        command("compute rdc all reduce max vx vy vz");
        command("fix ave all ave/time 1 1 1 c_ke c_rdc[2] mode scalar");
        command("fix ave1 all ave/time 1 1 1 c_ke mode scalar");
        command("variable eq equal 2.0*temp");
        command("variable vec vector [1.5,2.5,3.5]");
        command(
            "thermo_style custom step elapsed elaplong dt time cpu tpcpu spcpu cpuremain part "
            "timeremain "
            "atoms temp press pe ke etotal enthalpy evdwl ecoul epair ebond eangle edihed eimp "
            "emol elong "
            "etail vol density lx ly lz xlo xhi ylo yhi zlo zhi xy xz yz xlat ylat zlat "
            "bonds angles dihedrals impropers pxx pyy pzz pxy pxz pyz fmax fnorm nbuild ndanger "
            "cella cellb cellc cellalpha cellbeta cellgamma ecouple econserve "
            "c_ke c_rdc[2] f_ave1 f_ave[1] f_ave[2] v_eq v_vec[2]");
        command("thermo_modify format float %.8g");
    });
    auto output = run(5);
    ASSERT_MATCH(output, "Step +Elapsed +Elaplong +Dt +Time +CPU +T/CPU +S/CPU +CPULeft +Part "
                      "+TimeoutLeft +Atoms"
                      " +Temp +Press +PotEng +KinEng +TotEng +Enthalpy +E_vdwl +E_coul +E_pair "
                      "+E_bond +E_angle +E_dihed"
                      " +E_impro +E_mol +E_long +E_tail +Volume +Density +Lx +Ly +Lz +Xlo +Xhi +Ylo"
                      " +Yhi +Zlo +Zhi +Xy +Xz +Yz +Xlat +Ylat +Zlat +Bonds +Angles +Diheds"
                      " +Impros +Pxx +Pyy +Pzz +Pxy +Pxz +Pyz +Fmax +Fnorm +Nbuild +Ndanger +Cella"
                      " +Cellb +Cellc +CellAlpha +CellBeta +CellGamma +Ecouple +Econserve +c_ke "
                      "+c_rdc\\[2\\] +f_ave1 +f_ave\\[1\\] +f_ave\\[2\\]"
                      " +v_eq +v_vec\\[2\\]");

    // consistency checks through the variable interface of the thermo keywords
    ASSERT_NEAR(get_variable_value("eq"), 2.0 * lmp->input->variable->compute_equal("temp"), 1e-12);
    HIDE_OUTPUT([&] {
        command("variable natoms equal atoms");
        command("variable step equal step");
        command("variable elapsed equal elapsed");
        command("variable time equal time");
        command("variable dt equal dt");
        command("variable vol equal vol");
        command("variable lx equal lx");
        command("variable ly equal ly");
        command("variable lz equal lz");
        command("variable xy equal xy");
        command("variable cella equal cella");
        command("variable cellalpha equal cellalpha");
        command("variable dens equal density");
        command("variable etot equal etotal");
        command("variable pe equal pe");
        command("variable ke equal ke");
        command("variable cke equal c_ke");
        command("variable evdwl equal evdwl");
        command("variable epair equal epair");
        command("variable press equal press");
        command("variable ptrace equal (pxx+pyy+pzz)/3.0");
        command("variable nbonds equal bonds");
        command("variable fave equal f_ave1");
        command("variable fave1 equal f_ave[1]");
        command("variable fave2 equal f_ave[2]");
        command("variable crdc2 equal c_rdc[2]");
        command("variable vvec2 equal v_vec[2]");
        command("variable ecouple equal ecouple");
        command("variable econserve equal econserve");
    });
    ASSERT_EQ(get_variable_value("natoms"), 32.0);
    ASSERT_EQ(get_variable_value("step"), 5.0);
    ASSERT_EQ(get_variable_value("elapsed"), 5.0);
    ASSERT_EQ(get_variable_value("dt"), 0.005);
    ASSERT_NEAR(get_variable_value("time"), 0.025, 1e-14);
    ASSERT_NEAR(get_variable_value("lx"), 2.0 * pow(4.0 / 0.8442, 1.0 / 3.0), 1e-12);
    ASSERT_NEAR(get_variable_value("vol"),
                get_variable_value("lx") * get_variable_value("ly") * get_variable_value("lz"),
                1e-12);
    ASSERT_EQ(get_variable_value("xy"), 0.0);
    ASSERT_NEAR(get_variable_value("cella"), get_variable_value("lx"), 1e-12);
    ASSERT_NEAR(get_variable_value("cellalpha"), 90.0, 1e-12);
    ASSERT_NEAR(get_variable_value("dens"), 0.8442, 1e-12);
    ASSERT_NEAR(get_variable_value("etot"), get_variable_value("pe") + get_variable_value("ke"),
                1e-12);
    // with norm yes (lj units) the ke keyword is per atom, while a compute reference in a
    // variable always returns the unnormalized value
    ASSERT_NEAR(get_variable_value("ke"), get_variable_value("cke") / 32.0, 1e-12);
    ASSERT_NEAR(get_variable_value("evdwl"), get_variable_value("epair"), 1e-12);
    ASSERT_NEAR(get_variable_value("press"), get_variable_value("ptrace"), 1e-10);
    ASSERT_EQ(get_variable_value("nbonds"), 0.0);
    ASSERT_NEAR(get_variable_value("fave"), get_variable_value("cke"), 1e-12);
    ASSERT_NEAR(get_variable_value("fave1"), get_variable_value("cke"), 1e-12);
    ASSERT_NEAR(get_variable_value("fave2"), get_variable_value("crdc2"), 1e-12);
    ASSERT_EQ(get_variable_value("vvec2"), 2.5);
    ASSERT_EQ(get_variable_value("ecouple"), 0.0);
    ASSERT_NEAR(get_variable_value("econserve"), get_variable_value("etot"), 1e-12);

    // errors for invalid custom keywords and references
    TEST_FAILURE(".*ERROR: Unknown keyword 'xxx' in thermo_style custom command.*",
                 command("thermo_style custom step xxx"););
    TEST_FAILURE(".*ERROR: Could not find thermo custom compute ID: xxx.*",
                 command("thermo_style custom step c_xxx"););
    TEST_FAILURE(".*ERROR: Could not find thermo custom fix ID: xxx.*",
                 command("thermo_style custom step f_xxx"););
    TEST_FAILURE(".*ERROR: Could not find thermo custom variable name: xxx.*",
                 command("thermo_style custom step v_xxx"););
    TEST_FAILURE(".*ERROR: Thermo custom compute ke does not compute a vector.*",
                 command("thermo_style custom step c_ke[1]"););
    TEST_FAILURE(".*ERROR: Thermo custom compute rdc does not compute a scalar.*",
                 command("thermo_style custom step c_rdc"););
    TEST_FAILURE(".*ERROR: Thermo custom compute rdc vector is accessed out-of-range.*",
                 command("thermo_style custom step c_rdc[4]"););
    TEST_FAILURE(".*ERROR: Thermo custom fix ave vector is accessed out-of-range.*",
                 command("thermo_style custom step f_ave[3]"););
    TEST_FAILURE(".*ERROR: Thermo custom variable vec is not an equal-style variable.*",
                 command("thermo_style custom step v_vec"););
    TEST_FAILURE(".*ERROR: Thermo custom variable eq is not a vector-style variable.*",
                 command("thermo_style custom step v_eq[1]"););
    HIDE_OUTPUT([&] {
        command("variable xxx equal xxx");
    });
    TEST_FAILURE(".*ERROR: Variable xxx: Invalid thermo keyword .xxx. in variable formula.*",
                 get_variable_value("xxx"););
}

TEST_F(ThermoTest, Modify)
{
    auto *th = thermo();
    ASSERT_EQ(th->modified, 0);
    ASSERT_EQ(th->lostflag, Thermo::ERROR);
    ASSERT_EQ(th->lostbond, Thermo::ERROR);

    // norm: lj units default to yes, pe keyword is per-atom then
    HIDE_OUTPUT([&] {
        command("variable pe equal pe");
        command("run 0 post no");
    });
    ASSERT_EQ(th->normflag, 1);
    double pe_norm = get_variable_value("pe");
    HIDE_OUTPUT([&] {
        command("thermo_modify norm no");
        command("run 0 post no");
    });
    ASSERT_EQ(th->normflag, 0);
    ASSERT_EQ(th->modified, 1);
    ASSERT_NEAR(get_variable_value("pe"), 32.0 * pe_norm, 1e-10);
    HIDE_OUTPUT([&] {
        command("thermo_modify norm yes");
        command("run 0 post no");
    });
    ASSERT_EQ(th->normflag, 1);
    ASSERT_NEAR(get_variable_value("pe"), pe_norm, 1e-12);
    TEST_FAILURE(".*ERROR: Illegal thermo_modify norm command: missing argument.*",
                 command("thermo_modify norm"););

    // lost, lost/bond, warn
    HIDE_OUTPUT([&] {
        command("thermo_modify lost ignore lost/bond warn");
    });
    ASSERT_EQ(th->lostflag, Thermo::IGNORE);
    ASSERT_EQ(th->lostbond, Thermo::WARN);
    HIDE_OUTPUT([&] {
        command("thermo_modify lost warn lost/bond ignore");
    });
    ASSERT_EQ(th->lostflag, Thermo::WARN);
    ASSERT_EQ(th->lostbond, Thermo::IGNORE);
    HIDE_OUTPUT([&] {
        command("thermo_modify lost error lost/bond error");
    });
    ASSERT_EQ(th->lostflag, Thermo::ERROR);
    ASSERT_EQ(th->lostbond, Thermo::ERROR);
    TEST_FAILURE(".*ERROR: Unknown thermo_modify lost argument: xxx.*",
                 command("thermo_modify lost xxx"););
    TEST_FAILURE(".*ERROR: Unknown thermo_modify lost/bond argument: xxx.*",
                 command("thermo_modify lost/bond xxx"););

    HIDE_OUTPUT([&] {
        command("thermo_modify warn ignore");
    });
    ASSERT_EQ(lmp->error->get_maxwarn(), -1);
    HIDE_OUTPUT([&] {
        command("thermo_modify warn always");
    });
    ASSERT_EQ(lmp->error->get_maxwarn(), 0);
    HIDE_OUTPUT([&] {
        command("thermo_modify warn 5");
    });
    ASSERT_EQ(lmp->error->get_maxwarn(), 5);
    HIDE_OUTPUT([&] {
        command("thermo_modify warn reset");
    });
    ASSERT_EQ(lmp->error->get_maxwarn(), 5);
    ASSERT_EQ(lmp->error->get_numwarn(), 0);
    HIDE_OUTPUT([&] {
        command("thermo_modify warn default");
    });
    ASSERT_EQ(lmp->error->get_maxwarn(), 100);
    TEST_FAILURE(".*ERROR: Expected integer parameter instead of 'xxx'.*",
                 command("thermo_modify warn xxx"););

    // flush only affects buffering; both settings must work
    auto output = CAPTURE_OUTPUT([&] {
        command("thermo_modify flush yes");
        command("run 0 post no");
        command("thermo_modify flush no");
        command("run 0 post no");
    });
    ASSERT_MATCH(output, "Step +Temp");
    TEST_FAILURE(".*ERROR: Expected boolean parameter instead of 'xxx'.*",
                 command("thermo_modify flush xxx"););
    TEST_FAILURE(".*ERROR: Illegal thermo_modify flush command: missing argument.*",
                 command("thermo_modify flush"););
}

TEST_F(ThermoTest, Format)
{
    HIDE_OUTPUT([&] {
        command("thermo_style custom step atoms temp pe");
    });
    auto output = run0();
    ASSERT_MATCH(output, "\n +0 +32 +1 +-[0-9.]+ *\n");

    // format line: one format per column
    HIDE_OUTPUT([&] {
        command("thermo_modify format line \"%8d %6d T=%.3f PE=%.4e\"");
    });
    output = run0();
    ASSERT_MATCH(output, "\n +0 +32 T=1.000 PE=-[0-9]\\.[0-9][0-9][0-9][0-9]e[-+][0-9][0-9] *\n");

    // format int / float override types, format N overrides one column
    HIDE_OUTPUT([&] {
        command("thermo_modify format none");
        command("thermo_modify format int %04d format float %.2f");
    });
    output = run0();
    ASSERT_MATCH(output, "\n0000 +0032 +1.00 +-[0-9]\\.[0-9][0-9] *\n");
    HIDE_OUTPUT([&] {
        command("thermo_modify format 3 %10.5f");
        command("thermo_modify format -1 %.1e");
    });
    output = run0();
    ASSERT_MATCH(output, "\n0000 +0032 +1.00000 -[0-9]\\.[0-9]e[-+][0-9][0-9] *\n");

    // format none resets everything
    HIDE_OUTPUT([&] {
        command("thermo_modify format none");
    });
    output = run0();
    ASSERT_MATCH(output, "\n +0 +32 +1 +-[0-9.]+ *\n");

    // integer format for a bigint column gets the correct conversion specifier
    HIDE_OUTPUT([&] {
        command("thermo_modify format int %3d");
    });
    output = run0();
    ASSERT_MATCH(output, "\n +0 +32 +1 +-[0-9.]+ *\n");

    TEST_FAILURE(".*ERROR: Illegal thermo_modify format command: missing argument.*",
                 command("thermo_modify format"););
    TEST_FAILURE(".*ERROR: Illegal thermo_modify format command: missing argument.*",
                 command("thermo_modify format line"););
    TEST_FAILURE(".*ERROR: Invalid thermo_modify format argument: xxx.*",
                 command("thermo_modify format xxx %d"););
    TEST_FAILURE(".*ERROR: Invalid thermo_modify format argument: 10.*",
                 command("thermo_modify format 10 %d"););
    TEST_FAILURE(".*ERROR: Invalid thermo_modify format argument: -10.*",
                 command("thermo_modify format -10 %d"););
    TEST_FAILURE(".*ERROR: Thermo_modify int format does not contain a d conversion character.*",
                 command("thermo_modify format int %f"););
}

TEST_F(ThermoTest, Colname)
{
    HIDE_OUTPUT([&] {
        command("unfix nve");
        command("fix nvt all nvt temp 1.0 1.0 0.1");
        command("compute rdc all reduce max vx vy");
        command("thermo_style custom step temp pe f_nvt f_nvt[1] c_rdc[2]");
    });
    auto output = run0();
    ASSERT_MATCH(output, "Step +Temp +PotEng +f_nvt +f_nvt\\[1\\] +c_rdc\\[2\\] *\n");

    // by index, by negative index, by keyword
    HIDE_OUTPUT([&] {
        command("thermo_modify colname 1 Timestep colname -1 MaxVy colname temp Temperature");
    });
    output = run0();
    ASSERT_MATCH(output, "Timestep +Temperature +PotEng +f_nvt +f_nvt\\[1\\] +MaxVy *\n");

    // auto names from fix nvt and compute reduce
    HIDE_OUTPUT([&] {
        command("thermo_modify colname auto");
    });
    output = run0();
    ASSERT_MATCH(output, "Timestep +Temperature +PotEng +f_nvt:ecouple "
                                      "+f_nvt:eta\\[1\\] +c_rdc:max\\(vy\\) *\n");

    // default restores all
    HIDE_OUTPUT([&] {
        command("thermo_modify colname default");
    });
    output = run0();
    ASSERT_MATCH(output, "Step +Temp +PotEng +f_nvt +f_nvt\\[1\\] +c_rdc\\[2\\] *\n");

    // custom names are also used in multi and yaml format
    HIDE_OUTPUT([&] {
        command("thermo_modify colname pe Energy line multi");
    });
    output = run0();
    ASSERT_MATCH(output, "Energy += +-[0-9.]+");
    HIDE_OUTPUT([&] {
        command("thermo_modify line yaml");
    });
    output = run0();
    ASSERT_THAT(
        output,
        HasSubstr("keywords: ['Step', 'Temp', 'Energy', 'f_nvt', 'f_nvt[1]', 'c_rdc[2]', ]"));

    TEST_FAILURE(".*ERROR: Illegal thermo_modify colname command: missing argument.*",
                 command("thermo_modify colname"););
    TEST_FAILURE(".*ERROR: Illegal thermo_modify colname command: missing argument.*",
                 command("thermo_modify colname 1"););
    TEST_FAILURE(".*ERROR: Invalid thermo_modify colname argument: 7.*",
                 command("thermo_modify colname 7 xxx"););
    TEST_FAILURE(".*ERROR: Invalid thermo_modify colname argument: -7.*",
                 command("thermo_modify colname -7 xxx"););
    TEST_FAILURE(".*ERROR: Invalid thermo_modify colname argument: xxx.*",
                 command("thermo_modify colname xxx yyy"););
}

TEST_F(ThermoTest, TempPress)
{
    HIDE_OUTPUT([&] {
        command("group half id 1:16");
        command("compute alltemp all temp");
        command("compute halftemp half temp");
        command("compute ke all ke");
        command("compute mypress all pressure alltemp");
        command("variable temp equal temp");
        command("variable press equal press");
        command("run 0 post no");
    });
    double temp  = get_variable_value("temp");
    double press = get_variable_value("press");

    // replacing the temperature compute with an equivalent one changes nothing
    auto output = CAPTURE_OUTPUT([&] {
        command("thermo_modify temp alltemp");
        command("run 0 post no");
    });
    ASSERT_THAT(output, Not(HasSubstr("WARNING")));
    ASSERT_NEAR(get_variable_value("temp"), temp, 1e-12);
    ASSERT_NEAR(get_variable_value("press"), press, 1e-12);

    // a temperature compute for a subset of atoms is flagged and changes the pressure
    output = CAPTURE_OUTPUT([&] {
        command("thermo_modify temp halftemp");
        command("run 0 post no");
    });
    ASSERT_THAT(output, HasSubstr("WARNING: Temperature for thermo pressure is not for group all"));
    ASSERT_NE(get_variable_value("temp"), temp);
    ASSERT_NE(get_variable_value("press"), press);

    output = CAPTURE_OUTPUT([&] {
        command("thermo_modify temp alltemp press mypress");
        command("run 0 post no");
    });
    ASSERT_THAT(output, Not(HasSubstr("WARNING")));
    ASSERT_NEAR(get_variable_value("press"), press, 1e-12);

    TEST_FAILURE(".*ERROR: Could not find thermo_modify temperature compute xxx.*",
                 command("thermo_modify temp xxx"););
    TEST_FAILURE(".*ERROR: Thermo_modify compute ke does not compute temperature.*",
                 command("thermo_modify temp ke"););
    TEST_FAILURE(".*ERROR: Could not find thermo_modify pressure compute xxx.*",
                 command("thermo_modify press xxx"););
    TEST_FAILURE(".*ERROR: Thermo_modify compute ke does not compute pressure.*",
                 command("thermo_modify press ke"););
    TEST_FAILURE(".*ERROR: Illegal thermo_modify temp command: missing argument.*",
                 command("thermo_modify temp"););
    TEST_FAILURE(".*ERROR: Illegal thermo_modify press command: missing argument.*",
                 command("thermo_modify press"););

    // thermo styles without temperature or pressure output reject the keywords
    HIDE_OUTPUT([&] {
        command("thermo_style custom step pe");
    });
    TEST_FAILURE(".*ERROR: Thermo style does not use temp.*",
                 command("thermo_modify temp alltemp"););
    TEST_FAILURE(".*ERROR: Thermo style does not use press.*",
                 command("thermo_modify press mypress"););

    // restricted triclinic box: general triclinic output is rejected
    HIDE_OUTPUT([&] {
        command("thermo_style custom step temp press");
        command("change_box all triclinic");
    });
    TEST_FAILURE(".*ERROR: Thermo_modify triclinic/general cannot be used if simulation box is not "
                 "general triclinic.*",
                 command("thermo_modify triclinic/general yes"););
    HIDE_OUTPUT([&] {
        command("thermo_modify triclinic/general no");
    });
}

TEST_F(ThermoTest, LostAtoms)
{
    // replace the system with two atoms, one of which leaves the box through a fixed boundary
    HIDE_OUTPUT([&] {
        command("clear");
        command("units lj");
        command("atom_style atomic");
        command("boundary f p p");
        command("region box block 0 10 0 10 0 10");
        command("create_box 1 box");
        command("create_atoms 1 single 9.9 5.0 5.0");
        command("create_atoms 1 single 5.0 5.0 5.0");
        command("mass 1 1.0");
        command("velocity all set 1.0 0.0 0.0");
        command("pair_style zero 1.0");
        command("pair_coeff * *");
        command("fix nve all nve");
        command("thermo 10");
    });
    ASSERT_EQ(lmp->atom->natoms, 2);

    // default: error
    TEST_FAILURE(".*ERROR: Lost atoms: original 2 current 1.*", command("run 50 post no"););

    // warn: one warning, atom count updated
    HIDE_OUTPUT([&] {
        command("clear");
        command("units lj");
        command("atom_style atomic");
        command("boundary f p p");
        command("region box block 0 10 0 10 0 10");
        command("create_box 1 box");
        command("create_atoms 1 single 9.9 5.0 5.0");
        command("create_atoms 1 single 5.0 5.0 5.0");
        command("mass 1 1.0");
        command("velocity all set 1.0 0.0 0.0");
        command("pair_style zero 1.0");
        command("pair_coeff * *");
        command("fix nve all nve");
        command("thermo 10");
        command("thermo_modify lost warn");
    });
    auto output = run(50);
    ASSERT_THAT(output, HasSubstr("WARNING: Lost atoms: original 2 current 1"));
    ASSERT_EQ(lmp->atom->natoms, 1);
    // the warning is printed only once
    output = run(50);
    ASSERT_THAT(output, Not(HasSubstr("WARNING: Lost atoms")));

    // ignore: no message, atom count updated
    HIDE_OUTPUT([&] {
        command("clear");
        command("units lj");
        command("atom_style atomic");
        command("boundary f p p");
        command("region box block 0 10 0 10 0 10");
        command("create_box 1 box");
        command("create_atoms 1 single 9.9 5.0 5.0");
        command("create_atoms 1 single 5.0 5.0 5.0");
        command("mass 1 1.0");
        command("velocity all set 1.0 0.0 0.0");
        command("pair_style zero 1.0");
        command("pair_coeff * *");
        command("fix nve all nve");
        command("thermo 10");
        command("thermo_modify lost ignore");
    });
    output = run(50);
    ASSERT_THAT(output, Not(HasSubstr("Lost atoms")));
    ASSERT_EQ(lmp->atom->natoms, 1);
}

TEST_F(ThermoTest, TriclinicGeneral)
{
    // general triclinic box rotated by about 37 degrees around z relative to its restricted form
    HIDE_OUTPUT([&] {
        command("clear");
        command("units lj");
        command("atom_style atomic");
        command("lattice custom 1.0 a1 0.8 0.6 0.0 a2 -0.6 0.8 0.0 a3 0.0 0.3 1.0 basis 0.0 0.0 "
                "0.0 triclinic/general");
        command("create_box 1 NULL 0 3 0 3 0 3");
        command("create_atoms 1 box");
        command("mass 1 1.0");
        command("velocity all create 1.0 12345");
        command("pair_style lj/cut 2.5");
        command("pair_coeff 1 1 1.0 1.0");
        command("thermo_style custom step press pxx pyy pzz pxy pxz pyz avecx avecy avecz bvecx "
                "bvecy bvecz cvecx cvecy cvecz");
        command("thermo_modify format float %.10g");
        command("variable press equal press");
        command("variable pxx equal pxx");
        command("variable pyy equal pyy");
        command("variable pzz equal pzz");
        command("variable pxy equal pxy");
        command("variable pxz equal pxz");
        command("variable pyz equal pyz");
        command("variable avecx equal avecx");
        command("variable avecy equal avecy");
        command("variable bvecx equal bvecx");
        command("variable cvecy equal cvecy");
        command("variable cvecz equal cvecz");
    });
    ASSERT_TRUE(lmp->domain->triclinic_general);

    // restricted triclinic frame (default)
    auto output = run0();
    double p[6], q[6];
    const char *names[] = {"pxx", "pyy", "pzz", "pxy", "pxz", "pyz"};
    for (int i = 0; i < 6; ++i)
        p[i] = get_variable_value(names[i]);
    double press = get_variable_value("press");
    ASSERT_NEAR(press, (p[0] + p[1] + p[2]) / 3.0, 1e-10);
    ASSERT_NEAR(get_variable_value("avecx"), lmp->domain->xprd, 1e-12);
    ASSERT_EQ(get_variable_value("avecy"), 0.0);
    ASSERT_NEAR(get_variable_value("bvecx"), lmp->domain->xy, 1e-12);
    ASSERT_NEAR(get_variable_value("cvecz"), lmp->domain->zprd, 1e-12);
    // the printed columns agree with the variable values
    ASSERT_THAT(output, HasSubstr(fmt::format(" {:.10g} {:.10g} {:.10g} ", p[0], p[1], p[2])));

    // general triclinic frame: tensor rotated into the frame the box was created in
    HIDE_OUTPUT([&] {
        command("thermo_modify triclinic/general yes");
    });
    output = run0();
    for (int i = 0; i < 6; ++i)
        q[i] = get_variable_value(names[i]);
    ASSERT_THAT(output, HasSubstr(fmt::format(" {:.10g} {:.10g} {:.10g} ", q[0], q[1], q[2])));
    ASSERT_THAT(output, HasSubstr(fmt::format(" {:.10g} {:.10g} {:.10g} ", q[3], q[4], q[5])));
    // rotation about z: different in-plane components, same invariants, same zz component
    ASSERT_GT(fabs(q[0] - p[0]), 0.1);
    ASSERT_GT(fabs(q[3] - p[3]), 0.1);
    ASSERT_NEAR(get_variable_value("press"), press, 1e-10);
    ASSERT_NEAR(q[0] + q[1] + q[2], p[0] + p[1] + p[2], 1e-9);
    ASSERT_NEAR(q[2], p[2], 1e-9);
    ASSERT_NEAR(q[4] * q[4] + q[5] * q[5], p[4] * p[4] + p[5] * p[5], 1e-8);
    double fp = 0.0, fq = 0.0;
    for (int i = 0; i < 3; ++i) {
        fp += p[i] * p[i] + 2.0 * p[i + 3] * p[i + 3];
        fq += q[i] * q[i] + 2.0 * q[i + 3] * q[i + 3];
    }
    ASSERT_NEAR(fq, fp, 1e-7);
    // box edge vectors are reported in the general frame
    ASSERT_NEAR(get_variable_value("avecx"), 0.8 * 3.0, 1e-12);
    ASSERT_NEAR(get_variable_value("avecy"), 0.6 * 3.0, 1e-12);
    ASSERT_NEAR(get_variable_value("bvecx"), -0.6 * 3.0, 1e-12);
    ASSERT_NEAR(get_variable_value("cvecy"), 0.3 * 3.0, 1e-12);
    ASSERT_NEAR(get_variable_value("cvecz"), 3.0, 1e-12);

    // and back
    HIDE_OUTPUT([&] {
        command("thermo_modify triclinic/general no");
    });
    run0();
    for (int i = 0; i < 6; ++i)
        ASSERT_NEAR(get_variable_value(names[i]), p[i], 1e-12);
    ASSERT_NEAR(get_variable_value("avecx"), lmp->domain->xprd, 1e-12);
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
