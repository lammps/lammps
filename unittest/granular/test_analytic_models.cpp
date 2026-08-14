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
#include "math_const.h"
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

    if (cfg.analytic_model == "slip_cessation") {
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
    } else if (cfg.analytic_model == "rolling_decay") {
        // sphere (tag 1) spinning about +y on a flat wall, damped only by the
        // rolling-resistance torque M = mu_r R N (N = m g).  In the gross-rolling
        // (Coulomb-capped) regime the spin decays linearly:
        //   omega_y(t) = omega0 - (5 mu_r g)/(2 R) t.
        const double g      = var_or(vars, "grav", 0.0);
        const double mur    = var_or(vars, "mur", 0.0);
        const double omega0 = var_or(vars, "omega0", 0.0);
        const int i         = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "rolling_decay: atom with tag 1 not found";
        const double r = lmp->atom->radius[i];
        expect_rel(omega0 - (5.0 * mur * g) / (2.0 * r) * t, lmp->atom->omega[i][1], cfg.analytic_tol,
                   "rolling_decay omega_y");
    } else if (cfg.analytic_model == "pulloff_dmt") {
        // DMT cohesive contact held at (near-)zero overlap: the normal force on
        // particle (tag 1) has magnitude equal to the pull-off force
        //   F_pull = 4 pi gamma R_eff
        // since the Hertzian part ~0 at delta~0.  gamma = variable 'coh',
        // R_eff = variable 'reff' (= R for a wall, R/2 for two equal spheres).
        const double coh  = var_or(vars, "coh", 0.0);
        const double reff = var_or(vars, "reff", 0.0);
        const int i       = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "pulloff_dmt: atom with tag 1 not found";
        const double *f    = lmp->atom->f[i];
        const double fmag  = std::sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
        expect_rel(4.0 * MathConst::MY_PI * coh * reff, fmag, cfg.analytic_tol, "pulloff_dmt force");
    } else if (cfg.analytic_model == "collision_restitution") {
        // head-on collision of two equal spheres (tags 1,2) moving along x at
        // +/- vin: the total x-momentum stays ~0 and the relative normal speed
        // reverses with the coefficient of restitution e:
        //   -(v1x - v2x)/(2 vin) = e.
        const double en  = var_or(vars, "en", 1.0);
        const double vin = var_or(vars, "vin", 0.0);
        const int i1     = find_local(lmp, 1);
        const int i2     = find_local(lmp, 2);
        ASSERT_GE(i1, 0) << "collision_restitution: atom with tag 1 not found";
        ASSERT_GE(i2, 0) << "collision_restitution: atom with tag 2 not found";
        const double m1 = lmp->atom->rmass[i1];
        const double m2 = lmp->atom->rmass[i2];
        const double v1 = lmp->atom->v[i1][0];
        const double v2 = lmp->atom->v[i2][0];
        EXPECT_LE(std::fabs(m1 * v1 + m2 * v2), cfg.analytic_tol * (m1 + m2) * vin)
            << "collision_restitution: total x-momentum not conserved";
        expect_rel(en, -(v1 - v2) / (2.0 * vin), cfg.analytic_tol, "collision_restitution e");
    } else if (cfg.analytic_model == "hertz_normal_impact") {
        // Elastic Hertzian normal impact at peak compression (Chung & Ooi 2011,
        // Tests 1 & 2).  This segment must be timed to land at peak compression,
        // where the relative normal velocity is ~0 and all relative kinetic
        // energy is stored as Hertzian elastic potential energy.  For a contact
        // force F = K delta^{3/2}, that PE is (2/5) K delta^{5/2} = (2/5) F delta,
        // so  (1/2) mu_red V_rela^2 = (2/5) P_max alpha_max.  The 2/5 factor is
        // the signature of the 3/2 power law (a linear spring would give 1/2);
        // the relation is independent of the stiffness convention.  mu_red is
        // mred_factor * m (1/2 for two equal spheres, 1 for a sphere on a rigid
        // wall) and V_rela is the relative approach speed (2*vin for two spheres,
        // vin for a wall) -- both supplied via the variables block.
        const double vrela = var_or(vars, "vrela", 0.0);
        const double mredf = var_or(vars, "mred_factor", 1.0);
        const double floor = var_or(vars, "floor", 0.0);
        const int i        = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "hertz_normal_impact: atom with tag 1 not found";
        const double m    = lmp->atom->rmass[i];
        const double mred = mredf * m;
        // peak overlap: two-sphere (tags 1,2) or sphere-on-wall (single atom)
        double alpha;
        if (lmp->atom->nlocal >= 2) {
            const int j = find_local(lmp, 2);
            ASSERT_GE(j, 0) << "hertz_normal_impact: atom with tag 2 not found";
            const double dx = lmp->atom->x[i][0] - lmp->atom->x[j][0];
            const double dy = lmp->atom->x[i][1] - lmp->atom->x[j][1];
            const double dz = lmp->atom->x[i][2] - lmp->atom->x[j][2];
            alpha = lmp->atom->radius[i] + lmp->atom->radius[j] -
                std::sqrt(dx * dx + dy * dy + dz * dz);
        } else {
            alpha = lmp->atom->radius[i] - (lmp->atom->x[i][2] - floor);
        }
        ASSERT_GT(alpha, 0.0) << "hertz_normal_impact: atoms not in contact (segment "
                                 "not timed at peak compression?)";
        const double fx   = lmp->atom->f[i][0];
        const double fy   = lmp->atom->f[i][1];
        const double fz   = lmp->atom->f[i][2];
        const double pmax = std::sqrt(fx * fx + fy * fy + fz * fz);
        const double ke   = 0.5 * mred * vrela * vrela;
        const double pe   = 0.4 * pmax * alpha;    // (2/5) P_max alpha_max
        expect_rel(ke, pe, cfg.analytic_tol, "hertz_normal_impact peak energy balance");
    } else if (cfg.analytic_model == "spin_impact") {
        // A sphere (tag 1) impacts a much heavier/rigid partner -- a wall or a
        // large dense sphere -- approaching normally in -z at speed vin while
        // spinning about +y with omega0 (Chung & Ooi 2011, Tests 6 and 8).  The
        // spin gives a contact-point tangential velocity ~r*omega0 along x; in
        // the gross-sliding regime the Coulomb friction impulse mu(1+en) vin
        // (per unit mass) acts throughout contact, so for the spinning sphere
        //   vz_out = en vin (rebound),
        //   vx_out = mu (1+en) vin,
        //   omega_y_out = omega0 - (5/2) mu (1+en) vin / r.
        // Evaluate at a free-flight segment after rebound (needs omega0 large
        // enough that the contact never sticks).
        const double vin = var_or(vars, "vin", 0.0);
        const double w0  = var_or(vars, "omega0", 0.0);
        const double en  = var_or(vars, "en", 1.0);
        const double mu  = var_or(vars, "xmu", 0.0);
        const double r   = var_or(vars, "radius", 0.0);
        const int i      = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "spin_impact: atom with tag 1 not found";
        const double dvt = mu * (1.0 + en) * vin;
        expect_rel(en * vin, lmp->atom->v[i][2], cfg.analytic_tol, "spin_impact vz_out");
        expect_rel(dvt, lmp->atom->v[i][0], cfg.analytic_tol, "spin_impact vx_out");
        expect_rel(w0 - 2.5 * dvt / r, lmp->atom->omega[i][1], cfg.analytic_tol,
                   "spin_impact omega_y");
    } else {
        ADD_FAILURE() << "unknown analytic_model: '" << cfg.analytic_model << "'";
    }
}
