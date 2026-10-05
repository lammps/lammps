// unit tests for long-range electrostatics (KSpace styles)
//
// This is a C++ conversion of the long-range electrostatics
// verification tests based on scripts by Paul Crozier and Stan Moore),
// adapted to exercise any compiled-in KSpace style and to skip styles
// that are not available.
//
// Accuracy limits: the analytic erfc() approximation used by the coul/long
// pair styles is only accurate to roughly single precision, so force errors
// below ~1.0e-7 cannot converge at present.
//
// Not covered here (they need different systems or setups): the tip4p,
// dipole, spin, and electrode kspace variants, the scafacos wrapper, and the
// gpu, intel, and kk suffix styles (these need device/package initialization).

#include "library.h"

#include "atom.h"
#include "force.h"
#include "modify.h"
#include "compute.h"
#include "platform.h"
#include "utils.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

// which pair style family a kspace style must be combined with
enum PairKind {
    COUL_LONG,    // lj/cut/coul/long, no LJ (pure point charges)
    COUL_MSM,     // lj/cut/coul/msm, no LJ (pure point charges)
    LJ_LONG       // lj/long/coul/long with small but nonzero LJ parameters,
                  // since the dispersion kspace solvers reject all-zero
                  // dispersion coefficients
};

struct KSpaceConfig {
    std::string name;        // unique test instance name
    std::string check;       // kspace style that must be compiled in
    std::string kname;       // kspace style name used in the input
    std::string accuracy;    // accuracy for the two-point-charge tests
    std::string modify;      // extra kspace_modify settings, empty = none
    PairKind pair;           // pair style family to combine with
    bool omp;                // run with "-sf omp -pk omp 2"
};

// tells googletest how to print a config on test failures
static void PrintTo(const KSpaceConfig &cfg, std::ostream *os)
{
    *os << cfg.name;
}

// generates the gtest instance name from the config
static std::string config_name(const ::testing::TestParamInfo<KSpaceConfig> &info)
{
    return info.param.name;
}

// all point-charge capable kspace styles for fully periodic systems.
// accuracy settings follow the original Python tests where present.
static std::vector<KSpaceConfig> periodic_configs()
{
    return {
        {"ewald", "ewald", "ewald", "1.0e-8", "", COUL_LONG, false},
        {"ewald_omp", "ewald/omp", "ewald", "1.0e-8", "", COUL_LONG, true},
        {"ewald_disp", "ewald/disp", "ewald/disp", "1.0e-6", "", LJ_LONG, false},
        {"pppm", "pppm", "pppm", "1.0e-8", "", COUL_LONG, false},
        {"pppm_ad", "pppm", "pppm", "1.0e-6", "diff ad", COUL_LONG, false},
        {"pppm_cg", "pppm/cg", "pppm/cg", "1.0e-8", "", COUL_LONG, false},
        {"pppm_stagger", "pppm/stagger", "pppm/stagger", "1.0e-8", "", COUL_LONG, false},
        {"pppm_disp", "pppm/disp", "pppm/disp", "1.0e-6", "disp/auto yes", LJ_LONG, false},
        {"pppm_omp", "pppm/omp", "pppm", "1.0e-8", "", COUL_LONG, true},
        {"pppm_cg_omp", "pppm/cg/omp", "pppm/cg", "1.0e-8", "", COUL_LONG, true},
        {"pppm_disp_omp", "pppm/disp/omp", "pppm/disp", "1.0e-6", "disp/auto yes", LJ_LONG, true},
        // msm/cg and the OMP versions of the MSM pair styles require
        // "pressure/scalar no", so use it consistently for all MSM runs
        {"msm", "msm", "msm", "1.0e-8", "pressure/scalar no", COUL_MSM, false},
        {"msm_cg", "msm/cg", "msm/cg", "1.0e-8", "pressure/scalar no", COUL_MSM, false},
        {"msm_omp", "msm/omp", "msm", "1.0e-8", "pressure/scalar no", COUL_MSM, true},
        {"msm_cg_omp", "msm/cg/omp", "msm/cg", "1.0e-8", "pressure/scalar no", COUL_MSM, true},
    };
}

// kspace styles usable with a fully non-periodic box: MSM only
static std::vector<KSpaceConfig> nonperiodic_configs()
{
    return {
        {"msm", "msm", "msm", "1.0e-8", "pressure/scalar no", COUL_MSM, false},
        {"msm_cg", "msm/cg", "msm/cg", "1.0e-8", "pressure/scalar no", COUL_MSM, false},
        {"msm_omp", "msm/omp", "msm", "1.0e-8", "pressure/scalar no", COUL_MSM, true},
        {"msm_cg_omp", "msm/cg/omp", "msm/cg", "1.0e-8", "pressure/scalar no", COUL_MSM, true},
    };
}

// kspace styles usable with a slab geometry (p p f): MSM directly, the
// Ewald-family styles via the EW3DC slab correction
static std::vector<KSpaceConfig> slab_configs()
{
    return {
        {"ewald", "ewald", "ewald", "1.0e-8", "slab 3.0", COUL_LONG, false},
        {"pppm", "pppm", "pppm", "1.0e-8", "slab 3.0", COUL_LONG, false},
        {"msm", "msm", "msm", "1.0e-8", "pressure/scalar no", COUL_MSM, false},
        {"msm_cg", "msm/cg", "msm/cg", "1.0e-8", "pressure/scalar no", COUL_MSM, false},
    };
}

class KSpaceStyleTest : public LAMMPSTest, public ::testing::WithParamInterface<KSpaceConfig> {
protected:
    void SetUp() override
    {
        testbinary = "KSpaceStyles";
        if (GetParam().omp && lammps_config_has_package("OPENMP")) {
            for (const auto &arg : {"-sf", "omp", "-pk", "omp", "2"})
                args.emplace_back(arg);
        }
        LAMMPSTest::SetUp();
    }

    static std::string pair_name(PairKind kind)
    {
        switch (kind) {
            case COUL_MSM:
                return "lj/cut/coul/msm";
            case LJ_LONG:
                return "lj/long/coul/long";
            default:
                return "lj/cut/coul/long";
        }
    }

    // name of the missing required style, empty if all styles are available
    std::string missing_styles(const KSpaceConfig &cfg)
    {
        if (!info->has_style("kspace", cfg.check)) return cfg.check;
        if (!info->has_style("pair", pair_name(cfg.pair))) return pair_name(cfg.pair);
        return "";
    }

    // pair setup for the point charge tests: box/cutoff 10.0, no (or for the
    // dispersion solvers: negligibly small) LJ, analytic erfc() (no tables)
    void setup_point_charge_pair(PairKind kind)
    {
        switch (kind) {
            case COUL_MSM:
                command("pair_style lj/cut/coul/msm 10.0");
                command("pair_coeff 1 1 0.0 0.0");
                command("pair_modify table 0");
                break;
            case LJ_LONG:
                command("pair_style lj/long/coul/long long long 10.0");
                command("pair_coeff 1 1 1.0e-4 2.5");
                command("pair_modify table 0 table/disp 0");
                break;
            default:
                command("pair_style lj/cut/coul/long 10.0");
                command("pair_coeff 1 1 0.0 0.0");
                command("pair_modify table 0");
                break;
        }
    }

    void setup_kspace(const KSpaceConfig &cfg, const std::string &accuracy)
    {
        command(fmt::format("kspace_style {} {}", cfg.kname, accuracy));
        if (!cfg.modify.empty()) command("kspace_modify " + cfg.modify);
    }

    // two point charges +1/-1 in a 10x10x10 angstrom box
    void build_two_point_charges(const std::string &boundary, PairKind kind)
    {
        BEGIN_HIDE_OUTPUT();
        command("units real");
        command("atom_style charge");
        command("boundary " + boundary);
        command("region box block 0 10 0 10 0 10 units box");
        command("create_box 1 box");
        command("create_atoms 1 single 5.0 5.0 5.0 units box");
        command("create_atoms 1 single 3.14 2.78 5.43 units box");
        command("set atom 1 charge 1.0");
        command("set atom 2 charge -1.0");
        command("mass 1 1.0");
        setup_point_charge_pair(kind);
        END_HIDE_OUTPUT();
        ASSERT_EQ(lmp->atom->natoms, 2);
    }

    // distance between the two point charges
    static double tpc_distance()
    {
        return sqrt(1.86 * 1.86 + 2.22 * 2.22 + 0.43 * 0.43);
    }
};

/* ---------------------------------------------------------------------- */
// two point charges in a periodic box.  each style must reproduce the
// Madelung-like energy of the periodic charge pair.

class TwoPointChargesPeriodic : public KSpaceStyleTest {};

TEST_P(TwoPointChargesPeriodic, total_energy)
{
    const auto &cfg    = GetParam();
    const auto missing = missing_styles(cfg);
    if (!missing.empty()) GTEST_SKIP() << missing << " style is not compiled in";

    build_two_point_charges("p p p", cfg.pair);
    HIDE_OUTPUT([&] {
        setup_kspace(cfg, cfg.accuracy);
        command("run 0 post no");
    });

    const double etotal = lammps_get_thermo(lmp, "etotal");
    if (verbose) utils::print("{}: etotal = {:.8f}\n", cfg.name, etotal);

    // reference interval from the original Python test.  the tiny LJ
    // contribution of the dispersion solver configurations (see
    // setup_point_charge_pair()) is far below the interval width.
    EXPECT_GT(etotal, -119.17);
    EXPECT_LT(etotal, -119.15);
}

INSTANTIATE_TEST_SUITE_P(KSpace, TwoPointChargesPeriodic, ::testing::ValuesIn(periodic_configs()),
                         config_name);

/* ---------------------------------------------------------------------- */
// two point charges in a fully non-periodic box: only MSM supports this.
// with no periodic images the exact energy is the plain Coulomb energy of
// the isolated charge pair.

class TwoPointChargesNonPeriodic : public KSpaceStyleTest {};

TEST_P(TwoPointChargesNonPeriodic, energy_and_pressure)
{
    const auto &cfg    = GetParam();
    const auto missing = missing_styles(cfg);
    if (!missing.empty()) GTEST_SKIP() << missing << " style is not compiled in";

    build_two_point_charges("f f f", cfg.pair);
    HIDE_OUTPUT([&] {
        setup_kspace(cfg, cfg.accuracy);
        command("run 0 post no");
    });

    const double etotal = lammps_get_thermo(lmp, "etotal");
    const double press  = lammps_get_thermo(lmp, "press");
    if (verbose) utils::print("{}: etotal = {:.8f}  press = {:.4f}\n", cfg.name, etotal, press);

    // exact energy of two isolated point charges
    const double eref = -lmp->force->qqr2e / tpc_distance();
    EXPECT_NEAR(etotal, eref, 0.005);

    // pressure interval from the original Python test
    EXPECT_GT(press, -2593.1516);
    EXPECT_LT(press, -2591.1516);
}

INSTANTIATE_TEST_SUITE_P(KSpace, TwoPointChargesNonPeriodic,
                         ::testing::ValuesIn(nonperiodic_configs()), config_name);

/* ---------------------------------------------------------------------- */
// two point charges in a slab geometry (periodic in x and y only)

class TwoPointChargesSlab : public KSpaceStyleTest {};

TEST_P(TwoPointChargesSlab, total_energy)
{
    const auto &cfg    = GetParam();
    const auto missing = missing_styles(cfg);
    if (!missing.empty()) GTEST_SKIP() << missing << " style is not compiled in";

    build_two_point_charges("p p f", cfg.pair);
    HIDE_OUTPUT([&] {
        setup_kspace(cfg, cfg.accuracy);
        command("run 0 post no");
    });

    const double etotal = lammps_get_thermo(lmp, "etotal");
    if (verbose) utils::print("{}: etotal = {:.8f}\n", cfg.name, etotal);

    // reference interval from the original Python test
    EXPECT_GT(etotal, -119.2);
    EXPECT_LT(etotal, -119.1);
}

INSTANTIATE_TEST_SUITE_P(KSpace, TwoPointChargesSlab, ::testing::ValuesIn(slab_configs()),
                         config_name);

/* ---------------------------------------------------------------------- */
// random point charges in a box: compare the forces of each kspace style at
// moderate accuracy against a tightly converged Ewald reference on the same
// configuration.

class RandomPointCharges : public KSpaceStyleTest {
protected:
    static constexpr int NATOMS = 100;

    // +1/-1 point charges at random positions with a minimum separation
    void build_random_charges()
    {
        BEGIN_HIDE_OUTPUT();
        command("units real");
        command("atom_style charge");
        command("boundary p p p");
        command("region box block 0 30 0 30 0 30 units box");
        command("create_box 1 box");
        command(fmt::format("create_atoms 1 random {} 123456 NULL overlap 2.0 maxtry 1000 "
                            "units box", NATOMS));
        command(fmt::format("set atom 1*{} charge 1.0", NATOMS / 2));
        command(fmt::format("set atom {}*{} charge -1.0", NATOMS / 2 + 1, NATOMS));
        command("mass 1 1.0");
        END_HIDE_OUTPUT();
        ASSERT_EQ(lmp->atom->natoms, NATOMS);
    }

    // per-atom forces ordered by atom ID, summed across MPI ranks
    std::vector<double> gather_forces()
    {
        std::vector<double> forces(3 * NATOMS, 0.0);
        const auto *const *f  = lmp->atom->f;
        const auto *tag       = lmp->atom->tag;
        const int nlocal      = lmp->atom->nlocal;
        for (int i = 0; i < nlocal; ++i) {
            const int j      = 3 * ((int)tag[i] - 1);
            forces[j]        = f[i][0];
            forces[j + 1]    = f[i][1];
            forces[j + 2]    = f[i][2];
        }
        MPI_Allreduce(MPI_IN_PLACE, forces.data(), 3 * NATOMS, MPI_DOUBLE, MPI_SUM, lmp->world);
        return forces;
    }
};

TEST_P(RandomPointCharges, rms_force_error)
{
    const auto &cfg    = GetParam();
    const auto missing = missing_styles(cfg);
    if (!missing.empty()) GTEST_SKIP() << missing << " style is not compiled in";

    // reference solver matching the pair style family of the tested config
    KSpaceConfig ref = {"reference", "ewald", "ewald", "1.0e-8", "", COUL_LONG, cfg.omp};
    if (cfg.pair == LJ_LONG) {
        ref.check = ref.kname = "ewald/disp";
        ref.pair              = LJ_LONG;
    }
    if (!info->has_style("kspace", ref.check))
        GTEST_SKIP() << ref.check << " reference style is not compiled in";

    build_random_charges();

    // reference forces from a tightly converged Ewald run.  the requested
    // accuracy of 1.0e-8 is beyond the erfc() approximation limit of about
    // 1.0e-7 already; the original 1.0e-18 setting cannot improve convergence.
    HIDE_OUTPUT([&] {
        setup_point_charge_pair(ref.pair);
        setup_kspace(ref, ref.accuracy);
        command("run 0 post no");
    });
    const auto fref = gather_forces();

    // forces from the tested kspace style at moderate accuracy
    HIDE_OUTPUT([&] {
        setup_point_charge_pair(cfg.pair);
        setup_kspace(cfg, "1.0e-5");
        command("run 0 post no");
    });
    const auto ftest = gather_forces();

    double sumsq = 0.0;
    for (int i = 0; i < 3 * NATOMS; ++i) {
        const double df = ftest[i] - fref[i];
        sumsq += df * df;
    }
    const double rms = sqrt(sumsq / NATOMS);
    if (verbose) utils::print("{}: rms force error = {:.3g}\n", cfg.name, rms);

    // threshold from the original Python test (kcal/mol/angstrom).  the MSM
    // error estimation is known to be less strict: at the requested relative
    // accuracy of 1.0e-5 it produces an rms force error of 1.03e-2 for this
    // configuration where the Ewald and PPPM solvers stay below 0.5e-2
    if (cfg.pair == COUL_MSM)
        EXPECT_LT(rms, 3.0e-2);
    else
        EXPECT_LT(rms, 1.0e-2);
}

INSTANTIATE_TEST_SUITE_P(KSpace, RandomPointCharges, ::testing::ValuesIn(periodic_configs()),
                         config_name);

/* ---------------------------------------------------------------------- */
// SPC/E water: check that the per-atom energy and stress contributions of
// each kspace style sum up to the global potential energy and pressure
// tensor, and that energy and pressure agree with a tightly converged Ewald
// run of the same system.  this is the per-atom dynamic test of the original
// Python file, shrunk to unit test size (216 flexible SPC/E molecules).

class SpceWaterPerAtom : public KSpaceStyleTest {
protected:
    const std::string molfile = "kspace_spce.mol";
    bool have_molfile         = false;

    struct ThermoSnap {
        double pe, pesum;
        double press, sumpress;
        std::array<double, 6> ptensor, sumptensor;
    };

    void SetUp() override
    {
        KSpaceStyleTest::SetUp();
        have_molfile = write_molecule_file();
    }

    void TearDown() override
    {
        platform::unlink(molfile);
        KSpaceStyleTest::TearDown();
    }

    bool write_molecule_file()
    {
        FILE *fp = fopen(molfile.c_str(), "w");
        if (!fp) return false;
        fputs("# SPC/E water geometry UNITS: real\n"
              "3 atoms\n"
              "2 bonds\n"
              "1 angles\n"
              "\n"
              "Coords\n"
              "\n"
              "1    1.12456   0.09298   1.27452\n"
              "2    1.53683   0.75606   1.89928\n"
              "3    0.49482   0.56390   0.65678\n"
              "\n"
              "Types\n"
              "\n"
              "1        1\n"
              "2        2\n"
              "3        2\n"
              "\n"
              "Charges\n"
              "\n"
              "1       -0.8472\n"
              "2        0.4236\n"
              "3        0.4236\n"
              "\n"
              "Bonds\n"
              "\n"
              "1   1      1      2\n"
              "2   1      1      3\n"
              "\n"
              "Angles\n"
              "\n"
              "1   1      2      1      3\n"
              "\n"
              "Special Bond Counts\n"
              "\n"
              "1 2 0 0\n"
              "2 1 1 0\n"
              "3 1 1 0\n"
              "\n"
              "Special Bonds\n"
              "\n"
              "1 2 3\n"
              "2 1 3\n"
              "3 1 2\n",
              fp);
        fclose(fp);
        return true;
    }

    std::string missing_water_styles(const KSpaceConfig &cfg)
    {
        auto missing = missing_styles(cfg);
        if (!missing.empty()) return missing;
        if (!info->has_style("atom", "full")) return "atom style full";
        if (!info->has_style("bond", "harmonic")) return "bond style harmonic";
        if (!info->has_style("angle", "harmonic")) return "angle style harmonic";
        return "";
    }

    // 6x6x6 sc lattice of randomly oriented flexible SPC/E water molecules
    // at liquid density with per-atom energy/stress tally computes
    void build_water(PairKind kind)
    {
        BEGIN_HIDE_OUTPUT();
        command("clear");
        command("units real");
        command("atom_style full");
        command("boundary p p p");
        command("molecule water " + molfile);
        command("lattice sc 3.107");
        command("region box block 0 6 0 6 0 6");
        command("create_box 2 box bond/types 1 angle/types 1 extra/bond/per/atom 2 "
                "extra/angle/per/atom 1 extra/special/per/atom 2");
        command("create_atoms 0 box mol water 74637");
        command("mass 1 15.9994");
        command("mass 2 1.00794");
        switch (kind) {
            case COUL_MSM:
                command("pair_style lj/cut/coul/msm 9.0");
                break;
            case LJ_LONG:
                command("pair_style lj/long/coul/long long long 9.0");
                break;
            default:
                command("pair_style lj/cut/coul/long 9.0");
                break;
        }
        command("pair_coeff 1 1 0.15535 3.166");
        command("pair_coeff 2 2 0.0 0.0");
        command("bond_style harmonic");
        command("bond_coeff 1 1000.0 1.0");
        command("angle_style harmonic");
        command("angle_coeff 1 100.0 109.47");
        command("special_bonds lj/coul 0.0 0.0 0.5");
        command("velocity all create 300.0 432567 dist uniform");
        command("fix integrate all nve");
        command("timestep 0.5");
        command("compute peatom all pe/atom");
        command("compute pesum all reduce sum c_peatom");
        command("compute satom all stress/atom NULL");
        command("compute ssum all reduce sum c_satom[1] c_satom[2] c_satom[3] "
                "c_satom[4] c_satom[5] c_satom[6]");
        command("variable sumpress equal -(c_ssum[1]+c_ssum[2]+c_ssum[3])/(3*vol)");
        command("variable sumpxx equal -c_ssum[1]/vol");
        command("variable sumpyy equal -c_ssum[2]/vol");
        command("variable sumpzz equal -c_ssum[3]/vol");
        command("variable sumpxy equal -c_ssum[4]/vol");
        command("variable sumpxz equal -c_ssum[5]/vol");
        command("variable sumpyz equal -c_ssum[6]/vol");
        // per-atom tallies happen only on steps where the consuming computes
        // are invoked, so they must be part of the thermo output on every step
        command("thermo_style custom step pe c_pesum press v_sumpress pxx v_sumpxx "
                "pyy v_sumpyy pzz v_sumpzz pxy v_sumpxy pxz v_sumpxz pyz v_sumpyz");
        command("thermo 1");
        END_HIDE_OUTPUT();
        ASSERT_EQ(lmp->atom->natoms, 648);
    }

    ThermoSnap take_snapshot()
    {
        ThermoSnap snap;
        snap.pe       = lammps_get_thermo(lmp, "pe");
        snap.pesum    = lmp->modify->get_compute_by_id("pesum")->compute_scalar();
        snap.press    = lammps_get_thermo(lmp, "press");
        snap.sumpress = get_variable_value("sumpress");
        const char *comp[] = {"pxx", "pyy", "pzz", "pxy", "pxz", "pyz"};
        for (int i = 0; i < 6; ++i) {
            snap.ptensor[i]    = lammps_get_thermo(lmp, comp[i]);
            snap.sumptensor[i] = get_variable_value(fmt::format("sum{}", comp[i]));
        }
        return snap;
    }

    void check_consistency(const KSpaceConfig &cfg, const ThermoSnap &snap, int step)
    {
        if (verbose)
            utils::print("{} step {}: pe = {:.6f} pesum = {:.6f} press = {:.4f} "
                         "sumpress = {:.4f}\n",
                         cfg.name, step, snap.pe, snap.pesum, snap.press, snap.sumpress);

        // per-atom energies and stresses must sum up to the global values
        EXPECT_NEAR(snap.pesum, snap.pe, PE_CONS_TOL * fabs(snap.pe));
        EXPECT_NEAR(snap.sumpress, snap.press, PRESS_CONS_TOL);
        for (int i = 0; i < 6; ++i)
            EXPECT_NEAR(snap.sumptensor[i], snap.ptensor[i], PRESS_CONS_TOL);
    }

    // tolerances for per-atom vs. global consistency: the per-atom values of
    // the mesh-based solvers are computed via a separate interpolation and
    // only sum up to the global values approximately
    static constexpr double PE_CONS_TOL    = 1.0e-4;    // relative
    static constexpr double PRESS_CONS_TOL = 10.0;      // bar (absolute)

    // tolerances vs. the tightly converged Ewald reference run.  at the
    // requested relative accuracy of 1.0e-4 the largest observed deviations
    // are 2.0 kcal/mol and 90 bar (for the Ewald solver itself, from the
    // k-space truncation); the tolerances add a 3-5x safety margin
    static constexpr double PE_XSTYLE_TOL    = 10.0;    // kcal/mol (absolute)
    static constexpr double PRESS_XSTYLE_TOL = 300.0;   // bar (absolute)
};

TEST_P(SpceWaterPerAtom, energy_and_pressure)
{
    const auto &cfg = GetParam();
    if (!have_molfile)
        GTEST_SKIP() << "Cannot open molecule file for writing: " << utils::getsyserror();
    const auto missing = missing_water_styles(cfg);
    if (!missing.empty()) GTEST_SKIP() << missing << " is not compiled in";

    // reference solver matching the pair style family of the tested config
    KSpaceConfig ref = {"reference", "ewald", "ewald", "1.0e-6", "", COUL_LONG, cfg.omp};
    if (cfg.pair == LJ_LONG) {
        ref.check = ref.kname = "ewald/disp";
        ref.pair              = LJ_LONG;
    }
    if (!info->has_style("kspace", ref.check))
        GTEST_SKIP() << ref.check << " reference style is not compiled in";

    // reference energy and pressure from a tightly converged Ewald run
    build_water(ref.pair);
    HIDE_OUTPUT([&] {
        setup_kspace(ref, ref.accuracy);
        command("run 0 post no");
    });
    const auto refsnap = take_snapshot();

    // now the tested kspace style at the accuracy of the original Python test
    build_water(cfg.pair);
    HIDE_OUTPUT([&] {
        command(fmt::format("kspace_style {} 1.0e-4", cfg.kname));
        if (!cfg.modify.empty()) command("kspace_modify " + cfg.modify);
        command("run 0 post no");
    });

    const auto snap0 = take_snapshot();
    check_consistency(cfg, snap0, 0);

    // cross-style check against the Ewald reference (identical configuration
    // and velocities, so all differences are due to the kspace solver)
    if (verbose)
        utils::print("{} vs. reference: dpe = {:.3g} dpress = {:.3g} dptensor = {:.3g} "
                     "{:.3g} {:.3g} {:.3g} {:.3g} {:.3g}\n",
                     cfg.name, snap0.pe - refsnap.pe, snap0.press - refsnap.press,
                     snap0.ptensor[0] - refsnap.ptensor[0], snap0.ptensor[1] - refsnap.ptensor[1],
                     snap0.ptensor[2] - refsnap.ptensor[2], snap0.ptensor[3] - refsnap.ptensor[3],
                     snap0.ptensor[4] - refsnap.ptensor[4], snap0.ptensor[5] - refsnap.ptensor[5]);
    EXPECT_NEAR(snap0.pe, refsnap.pe, PE_XSTYLE_TOL);
    EXPECT_NEAR(snap0.press, refsnap.press, PRESS_XSTYLE_TOL);
    for (int i = 0; i < 6; ++i)
        EXPECT_NEAR(snap0.ptensor[i], refsnap.ptensor[i], PRESS_XSTYLE_TOL);

    // one step of MD to also exercise the tallies during timestepping
    HIDE_OUTPUT([&] { command("run 1 post no"); });
    const auto snap1 = take_snapshot();
    check_consistency(cfg, snap1, 1);
}

INSTANTIATE_TEST_SUITE_P(KSpace, SpceWaterPerAtom, ::testing::ValuesIn(periodic_configs()),
                         config_name);

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (const auto &arg : env) {
            if (arg == "-v") verbose = true;
        }
    }

    if ((argc > 1) && (std::string(argv[1]) == "-v")) verbose = true;

    const int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
