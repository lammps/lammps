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
    } else if (cfg.analytic_model == "bounce_height") {
        // hard-sphere limit: the apex (center) height after the k-th bounce is
        //   h_k = r + e^(2k) (h0 - r).
        // Evaluated at a free-flight segment via energy conservation,
        // apex = z + vz^2/(2g), so the segment need not end exactly at the apex.
        // The match is approximate (soft-sphere, finite stiffness) -> loose tol.
        const double g  = var_or(vars, "grav", 0.0);
        const double e  = var_or(vars, "restitution", var_or(vars, "en", 1.0));
        const double r  = var_or(vars, "radius", 0.0);
        const double h0 = var_or(vars, "h0", 0.0);
        const double k  = var_or(vars, "bounce_k", 1.0);
        const int i     = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "bounce_height: atom with tag 1 not found";
        const double z    = lmp->atom->x[i][2];
        const double vz   = lmp->atom->v[i][2];
        const double apex = z + vz * vz / (2.0 * g);
        expect_rel(r + std::pow(e, 2.0 * k) * (h0 - r), apex, cfg.analytic_tol, "bounce_height apex");
    } else if (cfg.analytic_model == "stack_energy") {
        // Two particles (tags 1 lower, 2 upper) stacked between a floor (ylo)
        // and ceiling (yhi).  For the elastic (e=1) linear-spring case the total
        // mechanical energy KE + gravitational PE + contact spring PE is
        // conserved; compare it to the initial value (particles start at rest).
        // Masses and radii are read from the live simulation so the comparison
        // does not depend on reproducing LAMMPS' mass formula.
        const double g    = var_or(vars, "grav", 0.0);
        const double kn   = var_or(vars, "knorm", 0.0);
        const double ylo  = var_or(vars, "ylo", 0.0);
        const double yhi  = var_or(vars, "yhi", 0.0);
        const double y1_0 = var_or(vars, "y1", 0.0);
        const double y2_0 = var_or(vars, "y2", 0.0);
        const int i1      = find_local(lmp, 1);
        const int i2      = find_local(lmp, 2);
        ASSERT_GE(i1, 0) << "stack_energy: atom with tag 1 not found";
        ASSERT_GE(i2, 0) << "stack_energy: atom with tag 2 not found";
        const double m1 = lmp->atom->rmass[i1];
        const double m2 = lmp->atom->rmass[i2];
        const double r1 = lmp->atom->radius[i1];
        const double r2 = lmp->atom->radius[i2];

        // linear-spring contact PE from the three possible overlaps:
        // lower particle vs floor, the pair, upper particle vs ceiling
        auto spring_pe = [&](double ya, double yb) {
            const double df  = std::max(0.0, r1 - (ya - ylo));
            const double dc  = std::max(0.0, (yb + r2) - yhi);
            const double dpp = std::max(0.0, (r1 + r2) - (yb - ya));
            return 0.5 * kn * (df * df + dc * dc + dpp * dpp);
        };

        const double e0 = (m1 * g * y1_0 + m2 * g * y2_0) + spring_pe(y1_0, y2_0);
        const double ya = lmp->atom->x[i1][1];
        const double yb = lmp->atom->x[i2][1];
        const double va = lmp->atom->v[i1][1];
        const double vb = lmp->atom->v[i2][1];
        const double ec = 0.5 * m1 * va * va + 0.5 * m2 * vb * vb + (m1 * g * ya + m2 * g * yb) +
                          spring_pe(ya, yb);
        expect_rel(e0, ec, cfg.analytic_tol, "stack_energy total energy");
    } else if (cfg.analytic_model == "slip_cessation") {
        // sphere (tag 1) launched along +x with no spin on a rough floor (normal +z):
        // kinetic friction decelerates u and spins it up about +y until the no-slip
        // condition u = omega_y r is reached.  Thereafter u = 5 u0/7 and omega_y = u/r.
        // Evaluate at a segment past the slip-cessation time t_s = 2 u0 / (7 mu g).
        const double u0 = var_or(vars, "u0", 0.0);
        const double r  = var_or(vars, "radius", 0.0);
        const int i     = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "slip_cessation: atom with tag 1 not found";
        const double u_final = 5.0 * u0 / 7.0;
        expect_rel(u_final, lmp->atom->v[i][0], cfg.analytic_tol, "slip_cessation vx");
        expect_rel(u_final / r, lmp->atom->omega[i][1], cfg.analytic_tol, "slip_cessation omega_y");
    } else if (cfg.analytic_model == "oblique_impact") {
        // grazing impact of sphere (tag 1) on a wall with normal +z, in the
        // gross-sliding regime (tangential velocity never reverses during
        // contact).  With approach velocity (vx_in, 0, -vz_in) and no spin:
        //   vz_out = en vz_in,  vx_out = vx_in - mu(1+en) vz_in,
        //   omega_y_out = (5/2) mu (1+en) vz_in / r.
        // Evaluate at a free-flight segment after the rebound.
        const double vx_in = var_or(vars, "vx_in", 0.0);
        const double vz_in = var_or(vars, "vz_in", 0.0);
        const double en    = var_or(vars, "en", 1.0);
        const double mu    = var_or(vars, "xmu", 0.0);
        const double r     = var_or(vars, "radius", 0.0);
        const int i        = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "oblique_impact: atom with tag 1 not found";
        const double dvt = mu * (1.0 + en) * vz_in;    // tangential velocity decrement
        expect_rel(en * vz_in, lmp->atom->v[i][2], cfg.analytic_tol, "oblique_impact vz_out");
        expect_rel(vx_in - dvt, lmp->atom->v[i][0], cfg.analytic_tol, "oblique_impact vx_out");
        expect_rel(2.5 * dvt / r, lmp->atom->omega[i][1], cfg.analytic_tol, "oblique_impact omega_y");
    } else {
        ADD_FAILURE() << "unknown analytic_model: '" << cfg.analytic_model << "'";
    }
}
