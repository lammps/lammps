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

// Closed-form ("analytic") reference solutions for the DEM verification
// tests.  Each model is selected by name from the YAML 'analytic_model' key
// and reads its parameters from the 'variables' block, so new functional
// forms can be exercised purely from YAML once implemented here.

#include "test_analytic_models.h"

#include "atom.h"
#include "lammps.h"
#include "update.h"

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
            // ignore non-numeric variables (e.g. NULL placeholders)
        }
    }
    return vars;
}

// local index of the atom with the given tag, or -1 if not present on this rank
static int find_local(LAMMPS *lmp, tagint id)
{
    for (int i = 0; i < lmp->atom->nlocal; ++i)
        if (lmp->atom->tag[i] == id) return i;
    return -1;
}

// relative-error comparison consistent with the regression epsilon checks
static void expect_rel(double expected, double actual, double tol, const std::string &what)
{
    const double denom  = std::max(std::fabs(expected), 1.0e-300);
    const double relerr = std::fabs(actual - expected) / denom;
    EXPECT_LE(relerr, tol) << what << ": expected " << expected << " got " << actual;
}

static double var_or(const std::map<std::string, double> &vars, const std::string &name,
                     double fallback)
{
    auto it = vars.find(name);
    return (it != vars.end()) ? it->second : fallback;
}

void check_analytic_model(const TestConfig &cfg, LAMMPS *lmp, int segment)
{
    if (!cfg.analytic_enable) return;

    int target = cfg.analytic_segment;
    if (target < 0) target = (int) cfg.run_segments.size() - 1;
    if (segment != target) return;

    const auto vars = as_doubles(cfg);
    const double t  = (double) lmp->update->ntimestep * lmp->update->dt;

    if (cfg.analytic_model == "freefall") {
        // ballistic motion of atom 1 before any contact: z = z0 - g t^2/2, vz = -g t.
        // velocity-Verlet integrates constant acceleration exactly, so this is tight.
        const double g  = var_or(vars, "grav", 0.0);
        const double z0 = var_or(vars, "z0", 0.0);
        const int i     = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "freefall: atom with tag 1 not found";
        expect_rel(z0 - 0.5 * g * t * t, lmp->atom->x[i][2], cfg.analytic_tol, "freefall z");
        expect_rel(-g * t, lmp->atom->v[i][2], cfg.analytic_tol, "freefall vz");
    } else {
        ADD_FAILURE() << "unknown analytic_model: '" << cfg.analytic_model << "'";
    }
}
