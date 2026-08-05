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

// Closed-form ("analytic") reference solutions for the BPM verification tests.
// Each model is selected by name from the YAML 'analytic_model' key and reads
// its parameters from the 'variables' block, so new functional forms can be
// exercised purely from YAML once implemented here.

#include "test_analytic_models.h"

#include "atom.h"
#include "lammps.h"

#include "gtest/gtest.h"

#include <cmath>
#include <map>
#include <string>

using namespace LAMMPS_NS;

// parse the variables block into name->double (non-numeric values are skipped)
static std::map<std::string, double> as_doubles(const TestConfig &cfg)
{
    std::map<std::string, double> vars;
    for (const auto &v : cfg.variables) {
        try {
            vars[v.first] = std::stod(v.second);
        } catch (std::exception &) {
            // ignore non-numeric variables
        }
    }
    return vars;
}

static double var_or(const std::map<std::string, double> &vars, const std::string &name,
                     double fallback)
{
    auto it = vars.find(name);
    return (it != vars.end()) ? it->second : fallback;
}

// local index of the atom with the given tag, or -1 if not present on this rank
static int find_local(LAMMPS *lmp, tagint id)
{
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        if (lmp->atom->tag[i] == id) return i;
    return -1;
}

// magnitude of the force on the atom with the given tag (0 if not on this rank)
static double force_mag(LAMMPS *lmp, tagint id)
{
    const int i = find_local(lmp, id);
    if (i < 0) return 0.0;
    double *f = lmp->atom->f[i];
    return std::sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
}

// current separation between atoms with tags id1 and id2 (single rank)
static double separation(LAMMPS *lmp, tagint id1, tagint id2)
{
    const int i = find_local(lmp, id1);
    const int j = find_local(lmp, id2);
    if ((i < 0) || (j < 0)) return 0.0;
    double **x = lmp->atom->x;
    const double dx = x[i][0] - x[j][0];
    const double dy = x[i][1] - x[j][1];
    const double dz = x[i][2] - x[j][2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static void expect_rel(double expected, double actual, double tol, const std::string &what)
{
    const double denom  = std::max(std::fabs(expected), 1.0e-300);
    const double relerr = std::fabs(actual - expected) / denom;
    EXPECT_LE(relerr, tol) << what << ": expected " << expected << " got " << actual;
}

// PMB single-bond force: for one bond between atoms 1 and 2 the equal-and-
// opposite force magnitude on each node is |f| = c * vfrac * |stretch|, with
// stretch = (r - r0)/r0 (the bond-direction normalization cancels the 1/r).
// The reference length r0 is the initial separation (variable r0); the current
// separation r is read live, so the law is checked at any run segment as the
// two nodes oscillate.  Variables provide c, vfrac and r0.
static void check_pmb_force(const TestConfig &cfg, LAMMPS *lmp)
{
    const auto vars    = as_doubles(cfg);
    const double c     = var_or(vars, "c", 0.0);
    const double vfrac = var_or(vars, "vfrac", 0.0);
    const double r0    = var_or(vars, "r0", 0.0);
    if (r0 <= 0.0) {
        ADD_FAILURE() << "pmb_force requires a positive r0 variable";
        return;
    }
    const double r        = separation(lmp, 1, 2);
    const double stretch  = (r - r0) / r0;
    const double expected = c * vfrac * std::fabs(stretch);

    expect_rel(expected, force_mag(lmp, 1), cfg.analytic_tol, "pmb_force atom 1");
    expect_rel(expected, force_mag(lmp, 2), cfg.analytic_tol, "pmb_force atom 2");
}

// value of a per-atom custom (d_) property for the atom with the given tag
static double custom_value(LAMMPS *lmp, const std::string &name, tagint id)
{
    int flag, cols;
    const int idx = lmp->atom->find_custom(name.c_str(), flag, cols);
    if ((idx < 0) || (flag != 1) || (cols != 0)) return 0.0;
    const int i = find_local(lmp, id);
    return (i < 0) ? 0.0 : lmp->atom->dvector[idx][i];
}

// LPS single-bond dilatation: with one bond the weighted volume is m = r0*vfrac,
// so theta = (3/m)*dr*vfrac = 3*dr/r0 = 3*stretch.  Reads the per-step dilatation
// from the d_theta property (which the YAML declares) and compares it with three
// times the live stretch; r0 is supplied as a variable.
static void check_lps_dilatation(const TestConfig &cfg, LAMMPS *lmp)
{
    const auto vars = as_doubles(cfg);
    const double r0 = var_or(vars, "r0", 0.0);
    if (r0 <= 0.0) {
        ADD_FAILURE() << "lps_dilatation requires a positive r0 variable";
        return;
    }
    const double stretch  = (separation(lmp, 1, 2) - r0) / r0;
    const double expected = 3.0 * stretch;
    expect_rel(expected, custom_value(lmp, "theta", 1), cfg.analytic_tol, "lps_dilatation atom 1");
    expect_rel(expected, custom_value(lmp, "theta", 2), cfg.analytic_tol, "lps_dilatation atom 2");
}

// LPS single-bond force: for one bond the weighted volume cancels the nodal
// volume, so the deviatoric (shear) term contributes nothing and the force
// magnitude reduces to the purely volumetric |f| = 18*K*|stretch|/r0 (the shear
// modulus G drops out for a single colinear bond -- there is no shear).  This
// complements check_lps_dilatation: that one validates the dilatation theta,
// this one validates the force assembled from theta and the weighted volume.
// Variables provide the bulk modulus kbulk and the reference length r0.
static void check_lps_force(const TestConfig &cfg, LAMMPS *lmp)
{
    const auto vars    = as_doubles(cfg);
    const double kbulk = var_or(vars, "kbulk", 0.0);
    const double r0    = var_or(vars, "r0", 0.0);
    if (r0 <= 0.0) {
        ADD_FAILURE() << "lps_force requires a positive r0 variable";
        return;
    }
    const double stretch  = (separation(lmp, 1, 2) - r0) / r0;
    const double expected = 18.0 * kbulk * std::fabs(stretch) / r0;

    expect_rel(expected, force_mag(lmp, 1), cfg.analytic_tol, "lps_force atom 1");
    expect_rel(expected, force_mag(lmp, 2), cfg.analytic_tol, "lps_force atom 2");
}

// bond_style bpm/zero: the defining property is that the bond exerts *no* force
// regardless of stretch, while still tracking and breaking bonds.  This checks
// both halves: the force on each bonded node is identically zero, and the live
// global bond count matches the expected value (variable nbonds_expect) so that
// the breaking machinery is exercised -- a run tuned to break the bond sets
// nbonds_expect 0, an intact run sets it to the initial bond count.
static void check_bpm_zero(const TestConfig &cfg, LAMMPS *lmp)
{
    const auto vars  = as_doubles(cfg);
    const double tol = cfg.analytic_tol;
    EXPECT_LE(force_mag(lmp, 1), tol) << "bpm_zero: force on atom 1 must be zero";
    EXPECT_LE(force_mag(lmp, 2), tol) << "bpm_zero: force on atom 2 must be zero";

    const double expect = var_or(vars, "nbonds_expect", -1.0);
    if (expect >= 0.0)
        EXPECT_EQ((double) lmp->atom->nbonds, expect) << "bpm_zero: unexpected surviving bond count";
}

// bond_style bpm/spring (smooth off, gamma 0, normalize off): a linear central
// spring, so the equal-and-opposite force magnitude is Hooke's law in the bond
// displacement, |f| = k*|r - r0|.  Variables provide the stiffness k and the
// reference length r0; the current separation r is read live.
static void check_spring_force(const TestConfig &cfg, LAMMPS *lmp)
{
    const auto vars = as_doubles(cfg);
    const double k  = var_or(vars, "k", 0.0);
    const double r0 = var_or(vars, "r0", 0.0);
    if (r0 <= 0.0) {
        ADD_FAILURE() << "spring_force requires a positive r0 variable";
        return;
    }
    const double expected = k * std::fabs(separation(lmp, 1, 2) - r0);
    expect_rel(expected, force_mag(lmp, 1), cfg.analytic_tol, "spring_force atom 1");
    expect_rel(expected, force_mag(lmp, 2), cfg.analytic_tol, "spring_force atom 2");
}

// bond_style bpm/spring/plastic (smooth off, gamma 0, normalize off): once the
// bond strain exceeds the plastic threshold eplastic, the plastic strain ep
// tracks the total strain so the elastic part saturates at eplastic, giving an
// elastic-perfectly-plastic force plateau |f| = k*eplastic*r0 independent of
// further stretch (the 1/r and the bond-direction projection cancel).  The deck
// must drive the bond past yield before the checked segment; variables provide
// k, eplastic and r0.
static void check_spring_plastic(const TestConfig &cfg, LAMMPS *lmp)
{
    const auto vars      = as_doubles(cfg);
    const double k       = var_or(vars, "k", 0.0);
    const double eplas   = var_or(vars, "eplastic", 0.0);
    const double r0      = var_or(vars, "r0", 0.0);
    const double expected = k * eplas * r0;
    expect_rel(expected, force_mag(lmp, 1), cfg.analytic_tol, "spring_plastic atom 1");
    expect_rel(expected, force_mag(lmp, 2), cfg.analytic_tol, "spring_plastic atom 2");
}

// bond_style bpm/peri, eps model -- cyclic plasticity on an octahedral cluster.
// A central node (tag 1) bonded to six axial neighbors is stretched along x past
// the deviatoric yield and unloaded back to the reference configuration (a fix
// move wiggle run over an integer half-period).  A single colinear bond has zero
// deviatoric extension (ed = dr - theta*r0/3 = 0) and can never yield -- which is
// why the two-body eps deck only ever exercises the elastic path; the central
// node's six non-colinear bonds carry a nonzero deviatoric force state, so the
// anisotropic stretch drives genuine plastic flow.
//
// The defining cyclic-plasticity signature is irreversibility: after a full
// load/unload the accumulated plastic multiplier lambda on the central node is
// strictly positive (an elastic cycle leaves it exactly zero), while lambda and
// every nodal force stay finite and bounded.  This is the regression guard for the
// EPS radial-return map (bond_bpm_peri.cpp): the corrected projection lands on the
// yield surface in one step and keeps lambda bounded, whereas the legacy
// edp += rkNew*deltalambda update (~2/r0 gain) overshoots and diverges -- which
// here would surface as a blown-up lambda or non-finite forces.  Variable
// lambda_min sets the positive floor that proves plasticity engaged.
static void check_eps_cyclic(const TestConfig &cfg, LAMMPS *lmp)
{
    const auto vars         = as_doubles(cfg);
    const double lambda_min = var_or(vars, "lambda_min", 1.0e-8);

    const double lambda_c = custom_value(lmp, "lambda", 1);    // central node (tag 1)
    EXPECT_GT(lambda_c, lambda_min)
        << "eps_cyclic: central node did not yield (lambda=" << lambda_c
        << ") -- the deviatoric plastic path was not exercised";
    EXPECT_TRUE(std::isfinite(lambda_c))
        << "eps_cyclic: plastic multiplier diverged (lambda=" << lambda_c
        << ") -- the EPS radial-return map is unstable";

    // a divergent return map blows the deviatoric force state up to non-finite
    // values; every node in the cluster must keep a finite force
    for (tagint id = 1; id <= 7; ++id) {
        const int i = find_local(lmp, id);
        if (i < 0) continue;
        const double *f = lmp->atom->f[i];
        EXPECT_TRUE(std::isfinite(f[0]) && std::isfinite(f[1]) && std::isfinite(f[2]))
            << "eps_cyclic: non-finite force on node " << id;
    }
}

void check_analytic_model(const TestConfig &cfg, LAMMPS *lmp, int segment)
{
    if (!cfg.analytic_enable) return;
    // a negative analytic_segment means "the last segment"
    const int target =
        (cfg.analytic_segment < 0) ? (int) cfg.run_segments.size() - 1 : cfg.analytic_segment;
    if (segment != target) return;

    if (cfg.analytic_model == "pmb_force") {
        check_pmb_force(cfg, lmp);
    } else if (cfg.analytic_model == "lps_dilatation") {
        check_lps_dilatation(cfg, lmp);
    } else if (cfg.analytic_model == "lps_force") {
        check_lps_force(cfg, lmp);
    } else if (cfg.analytic_model == "bpm_zero") {
        check_bpm_zero(cfg, lmp);
    } else if (cfg.analytic_model == "spring_force") {
        check_spring_force(cfg, lmp);
    } else if (cfg.analytic_model == "spring_plastic") {
        check_spring_plastic(cfg, lmp);
    } else if (cfg.analytic_model == "eps_cyclic") {
        check_eps_cyclic(cfg, lmp);
    } else {
        ADD_FAILURE() << "unknown analytic_model: " << cfg.analytic_model;
    }
}
