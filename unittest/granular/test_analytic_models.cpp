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

    if (cfg.analytic_model == "freefall") {
        // ballistic motion of atom 1 before any contact: z = z0 - g t^2/2, vz = -g t.
        // velocity-Verlet integrates constant acceleration exactly, so this is tight.
        const double g  = var_or(vars, "grav", 0.0);
        const double z0 = var_or(vars, "z0", 0.0);
        const int i     = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "freefall: atom with tag 1 not found";
        expect_rel(z0 - 0.5 * g * t * t, lmp->atom->x[i][2], cfg.analytic_tol, "freefall z");
        expect_rel(-g * t, lmp->atom->v[i][2], cfg.analytic_tol, "freefall vz");
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
    } else if (cfg.analytic_model == "energy_dissipation") {
        // For a frictional collision with no external forces (no gravity) the
        // total mechanical (translational + rotational) kinetic energy of the
        // sphere (tag 1) must not increase.  This guards against the
        // grazing-impact energy-injection bug of the classic tangential model.
        // Initial state: velocity (vx_in, 0, -vz_in) with no spin; sphere moment
        // of inertia I = (2/5) m r^2.  analytic_tol is the (small) fractional
        // excess over the initial energy that is tolerated.
        const double vx_in = var_or(vars, "vx_in", 0.0);
        const double vz_in = var_or(vars, "vz_in", 0.0);
        const int i        = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "energy_dissipation: atom with tag 1 not found";
        const double m  = lmp->atom->rmass[i];
        const double r  = lmp->atom->radius[i];
        const double *v = lmp->atom->v[i];
        const double *w = lmp->atom->omega[i];
        const double e_init  = 0.5 * m * (vx_in * vx_in + vz_in * vz_in);
        const double ke_tr   = 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        const double ke_rot  = 0.5 * (0.4 * m * r * r) * (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
        const double e_final = ke_tr + ke_rot;
        EXPECT_LE(e_final, e_init * (1.0 + cfg.analytic_tol))
            << "energy_dissipation: final energy " << e_final << " exceeds initial " << e_init;
    } else if (cfg.analytic_model == "terminal_velocity_linear") {
        // particle (tag 1) falling under gravity g with linear (Stokes) drag
        // F = -gamma v reaches terminal speed v_term = m g / gamma.
        const double g     = var_or(vars, "grav", 0.0);
        const double gamma = var_or(vars, "gamma", 0.0);
        const int i        = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "terminal_velocity_linear: atom with tag 1 not found";
        const double m = lmp->atom->rmass[i];
        expect_rel(m * g / gamma, -lmp->atom->v[i][2], cfg.analytic_tol,
                   "terminal_velocity_linear");
    } else if (cfg.analytic_model == "terminal_velocity_schiller_naumann") {
        // particle (tag 1) falling under gravity g with Schiller-Naumann drag
        // (quiescent gas): terminal speed solves m g = 1/2 Cd rho_g pi r^2 v^2.
        const double g   = var_or(vars, "grav", 0.0);
        const double rho = var_or(vars, "rho_gas", 0.0);
        const double mu  = var_or(vars, "mu_gas", 0.0);
        const int i      = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "terminal_velocity_schiller_naumann: atom with tag 1 not found";
        const double m    = lmp->atom->rmass[i];
        const double r    = lmp->atom->radius[i];
        const double area = MathConst::MY_PI * r * r;
        const double mg   = m * g;
        auto drag         = [&](double v) {
            if (v <= 0.0) return 0.0;
            const double re = rho * v * (2.0 * r) / mu;
            const double cd = (24.0 / re) * (1.0 + 0.15 * std::pow(re, 0.687));
            return 0.5 * cd * rho * area * v * v;
        };
        // bracket then bisect for the terminal speed (drag is monotone in v)
        double vlo = 0.0, vhi = 1.0;
        int guard = 0;
        while ((drag(vhi) < mg) && (guard++ < 200)) vhi *= 2.0;
        // without a valid bracket the bisection below would return a meaningless
        // speed, so fail here instead.  this also catches missing or zero
        // rho_gas / mu_gas settings, for which drag() is not a number.
        ASSERT_GE(drag(vhi), mg) << "terminal_velocity_schiller_naumann: could not bracket the "
                                    "terminal speed; check the grav, rho_gas, and mu_gas settings";
        for (int it = 0; it < 100; ++it) {
            const double vm = 0.5 * (vlo + vhi);
            if (drag(vm) < mg)
                vlo = vm;
            else
                vhi = vm;
        }
        expect_rel(0.5 * (vlo + vhi), -lmp->atom->v[i][2], cfg.analytic_tol,
                   "terminal_velocity_schiller_naumann");
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
    } else if (cfg.analytic_model == "oblique_impact_pair") {
        // symmetric oblique impact of two equal spheres (tags 1,2) in the
        // gross-sliding regime (Chung & Ooi 2011, Test 5).  Sphere 1 starts at
        // (-sep,-y0,0) with velocity (+vn_in,+vt_in,0), sphere 2 mirrored, with
        // y0 = vt_in*(sep-radius)/vn_in so the line of centers is along x at
        // first touch.  The per-sphere impulse solution mirrors the wall case
        // (the reduced mass m/2 cancels between normal and tangential impulse):
        //   v1x' = -en vn_in,  v1y' = vt_in - mu(1+en) vn_in,
        //   omega1z' = omega2z' = -(5/2) mu(1+en) vn_in / r,
        // and sphere 2 keeps v2 = -v1 by symmetry.  Gross sliding throughout
        // requires vt_in > (7/2) mu (1+en) vn_in.  Evaluate at a free-flight
        // segment after the rebound.
        const double vn = var_or(vars, "vn_in", 0.0);
        const double vt = var_or(vars, "vt_in", 0.0);
        const double en = var_or(vars, "en", 1.0);
        const double mu = var_or(vars, "xmu", 0.0);
        const double r  = var_or(vars, "radius", 0.0);
        const int i1    = find_local(lmp, 1);
        const int i2    = find_local(lmp, 2);
        ASSERT_GE(i1, 0) << "oblique_impact_pair: atom with tag 1 not found";
        ASSERT_GE(i2, 0) << "oblique_impact_pair: atom with tag 2 not found";
        const double dvt = mu * (1.0 + en) * vn;    // tangential velocity decrement
        expect_rel(-en * vn, lmp->atom->v[i1][0], cfg.analytic_tol, "oblique_impact_pair v1x");
        expect_rel(vt - dvt, lmp->atom->v[i1][1], cfg.analytic_tol, "oblique_impact_pair v1y");
        expect_rel(en * vn, lmp->atom->v[i2][0], cfg.analytic_tol, "oblique_impact_pair v2x");
        expect_rel(-(vt - dvt), lmp->atom->v[i2][1], cfg.analytic_tol, "oblique_impact_pair v2y");
        expect_rel(-2.5 * dvt / r, lmp->atom->omega[i1][2], cfg.analytic_tol,
                   "oblique_impact_pair omega1z");
        expect_rel(-2.5 * dvt / r, lmp->atom->omega[i2][2], cfg.analytic_tol,
                   "oblique_impact_pair omega2z");
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
    } else if (cfg.analytic_model == "spin_no_friction") {
        // Two identical spheres (tags 1,2) collide head-on along z while spinning
        // about y with equal and opposite omega0, arranged so the relative
        // tangential velocity at the contact point is zero (Chung & Ooi 2011,
        // Test 7).  No tangential force should be generated: each sphere's spin
        // is preserved and it gains no tangential (x,y) velocity.  This guards
        // against a model spuriously creating friction from spin alone.
        const double w0  = var_or(vars, "omega0", 0.0);
        const double tol = cfg.analytic_tol;
        const int i      = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "spin_no_friction: atom with tag 1 not found";
        expect_rel(w0, lmp->atom->omega[i][1], tol, "spin_no_friction omega_y preserved");
        EXPECT_LE(std::fabs(lmp->atom->v[i][0]), tol) << "spin_no_friction: spurious vx";
        EXPECT_LE(std::fabs(lmp->atom->v[i][1]), tol) << "spin_no_friction: spurious vy";
    } else if (cfg.analytic_model == "hertz_peak") {
        // Per-quantity check of an elastic Hertzian impact at peak compression
        // (Chung & Ooi 2011, Tests 1 and 2, Eqs. 2 and 3).  With the contact
        // force written explicitly as F = kfac delta^{3/2} (kfac supplied via
        // the variables block, spelling out the stiffness convention of the
        // model under test), energy conservation gives the peak overlap and
        // peak force separately:
        //   alpha_max = (5 mu_red V_rela^2 / (4 kfac))^{2/5}
        //   P_max     = kfac alpha_max^{3/2}
        // The segment must be timed to land at peak compression, like
        // hertz_normal_impact (which checks only the combined energy balance).
        const double vrela = var_or(vars, "vrela", 0.0);
        const double mredf = var_or(vars, "mred_factor", 1.0);
        const double kfac  = var_or(vars, "kfac", 0.0);
        const double floor = var_or(vars, "floor", 0.0);
        const int i        = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "hertz_peak: atom with tag 1 not found";
        const double mred = mredf * lmp->atom->rmass[i];
        // measured peak overlap: two-sphere (tags 1,2) or sphere-on-wall
        double alpha;
        if (lmp->atom->nlocal >= 2) {
            const int j = find_local(lmp, 2);
            ASSERT_GE(j, 0) << "hertz_peak: atom with tag 2 not found";
            const double dx = lmp->atom->x[i][0] - lmp->atom->x[j][0];
            const double dy = lmp->atom->x[i][1] - lmp->atom->x[j][1];
            const double dz = lmp->atom->x[i][2] - lmp->atom->x[j][2];
            alpha = lmp->atom->radius[i] + lmp->atom->radius[j] -
                std::sqrt(dx * dx + dy * dy + dz * dz);
        } else {
            alpha = lmp->atom->radius[i] - (lmp->atom->x[i][2] - floor);
        }
        ASSERT_GT(alpha, 0.0) << "hertz_peak: atoms not in contact (segment "
                                 "not timed at peak compression?)";
        const double fx    = lmp->atom->f[i][0];
        const double fy    = lmp->atom->f[i][1];
        const double fz    = lmp->atom->f[i][2];
        const double pmax  = std::sqrt(fx * fx + fy * fy + fz * fz);
        const double apred = std::pow(5.0 * mred * vrela * vrela / (4.0 * kfac), 0.4);
        expect_rel(apred, alpha, cfg.analytic_tol, "hertz_peak alpha_max");
        expect_rel(kfac * std::pow(apred, 1.5), pmax, cfg.analytic_tol, "hertz_peak P_max");
    } else if (cfg.analytic_model == "slip_transient") {
        // sphere (tag 1) launched along +x on a rough floor: DURING the sliding
        // phase (t < t_s = 2 u0 / (7 mu g)) kinetic friction gives the exact
        // linear laws  u(t) = u0 - mu g t  and  omega_y(t) = (5/2) (mu g / r) t.
        // Together they also pin down the slip-cessation time t_s where the two
        // lines meet.  The segment boundary must lie strictly inside the
        // sliding phase.
        const double u0 = var_or(vars, "u0", 0.0);
        const double mu = var_or(vars, "xmu", 0.0);
        const double g  = var_or(vars, "grav", 0.0);
        const int i     = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "slip_transient: atom with tag 1 not found";
        const double r  = lmp->atom->radius[i];
        const double ts = 2.0 * u0 / (7.0 * mu * g);
        ASSERT_LT(t, ts) << "slip_transient: segment boundary is not inside the sliding phase";
        expect_rel(u0 - mu * g * t, lmp->atom->v[i][0], cfg.analytic_tol, "slip_transient u");
        expect_rel(2.5 * mu * g * t / r, lmp->atom->omega[i][1], cfg.analytic_tol,
                   "slip_transient omega_y");
    } else if (cfg.analytic_model == "incline_rolling") {
        // sphere (tag 1) released at rest on an inclined floor (gravity vector
        // (sin_t, 0, -cos_t) * g, wall normal +z) with rolling resistance mu_r
        // (rolling sds, Coulomb-capped).  Rolling without slipping down the
        // incline:  a = (5/7) g (sin(theta) - mu_r cos(theta)), v = a t,
        // omega_y = v / r.  For mu_r >= tan(theta) the sphere stays at rest
        // (the model then bounds |v| and |omega| by analytic_tol, absolute).
        const double g    = var_or(vars, "grav", 0.0);
        const double mur  = var_or(vars, "mur", 0.0);
        const double sint = var_or(vars, "sin_t", 0.0);
        const double cost = var_or(vars, "cos_t", 1.0);
        const int i       = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "incline_rolling: atom with tag 1 not found";
        const double r     = lmp->atom->radius[i];
        const double accel = (5.0 / 7.0) * g * (sint - mur * cost);
        if (accel > 0.0) {
            expect_rel(accel * t, lmp->atom->v[i][0], cfg.analytic_tol, "incline_rolling v");
            expect_rel(accel * t / r, lmp->atom->omega[i][1], cfg.analytic_tol,
                       "incline_rolling omega_y");
        } else {
            EXPECT_LE(std::fabs(lmp->atom->v[i][0]), cfg.analytic_tol)
                << "incline_rolling: sphere should stay at rest (v)";
            EXPECT_LE(std::fabs(lmp->atom->omega[i][1] * r), cfg.analytic_tol)
                << "incline_rolling: sphere should stay at rest (omega)";
        }
    } else if (cfg.analytic_model == "wall_restitution") {
        // sphere (tag 1) launched along +x bounces off a wall (plane or region)
        // and returns with vx_out = -e vx_in.  Verifies the restitution of
        // wall styles for which no other closed form applies (e.g. region
        // walls); evaluate at a free-flight segment after the rebound.
        const double en    = var_or(vars, "en", 1.0);
        const double vx_in = var_or(vars, "vx_in", 0.0);
        const int i        = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "wall_restitution: atom with tag 1 not found";
        expect_rel(-en * vx_in, lmp->atom->v[i][0], cfg.analytic_tol, "wall_restitution vx_out");
    } else if (cfg.analytic_model == "momentum_conservation") {
        // two particles (tags 1,2), no external forces: total linear momentum
        // and total angular momentum about the origin (orbital m r x v plus
        // per-particle angular momentum for ellipsoid-type particles) are
        // conserved.  The initial state is particle 1 at (x1,y1,z1) moving
        // with (vx,0,0) and particle 2 at rest without spin, so
        //   P0 = m1 (vx,0,0),   L0 = m1 (x1,y1,z1) x (vx,0,0).
        // Requires a non-periodic box (no minimum-image jumps in r x v).
        const double vx = var_or(vars, "vx", 0.0);
        const double x1 = var_or(vars, "x1", 0.0);
        const double y1 = var_or(vars, "y1", 0.0);
        const double z1 = var_or(vars, "z1", 0.0);
        const int i1    = find_local(lmp, 1);
        const int i2    = find_local(lmp, 2);
        ASSERT_GE(i1, 0) << "momentum_conservation: atom with tag 1 not found";
        ASSERT_GE(i2, 0) << "momentum_conservation: atom with tag 2 not found";
        const double m1 = lmp->atom->rmass[i1];
        const double m2 = lmp->atom->rmass[i2];
        double ptot[3], ltot[3];
        for (int k = 0; k < 3; ++k)
            ptot[k] = m1 * lmp->atom->v[i1][k] + m2 * lmp->atom->v[i2][k];
        for (int k = 0; k < 3; ++k) {
            const int ka = (k + 1) % 3, kb = (k + 2) % 3;
            ltot[k] = m1 * (lmp->atom->x[i1][ka] * lmp->atom->v[i1][kb] -
                            lmp->atom->x[i1][kb] * lmp->atom->v[i1][ka]) +
                      m2 * (lmp->atom->x[i2][ka] * lmp->atom->v[i2][kb] -
                            lmp->atom->x[i2][kb] * lmp->atom->v[i2][ka]);
            if (lmp->atom->angmom_flag)
                ltot[k] += lmp->atom->angmom[i1][k] + lmp->atom->angmom[i2][k];
        }
        const double p0[3] = {m1 * vx, 0.0, 0.0};
        const double l0[3] = {0.0, m1 * z1 * vx, -m1 * y1 * vx};
        (void) x1;    // enters L0 only via components that vanish for v = (vx,0,0)
        const double pscale = std::fabs(m1 * vx);
        const double lscale = std::fabs(m1 * vx) * std::sqrt(y1 * y1 + z1 * z1);
        for (int k = 0; k < 3; ++k) {
            EXPECT_LE(std::fabs(ptot[k] - p0[k]), cfg.analytic_tol * pscale)
                << "momentum_conservation: linear momentum component " << k;
            EXPECT_LE(std::fabs(ltot[k] - l0[k]), cfg.analytic_tol * lscale)
                << "momentum_conservation: angular momentum component " << k;
        }
    } else if (cfg.analytic_model == "twist_decay") {
        // sphere (tag 1) resting on a floor (normal +z, N = m g) spinning about
        // the contact normal with omega0, damped only by the Coulomb-capped
        // twisting torque of the sds model, M = mu_t N (mu_t carries length
        // units).  While the cap is active the spin decays linearly:
        //   omega_z(t) = omega0 - (5 mu_t g)/(2 r^2) t.
        const double g      = var_or(vars, "grav", 0.0);
        const double mut    = var_or(vars, "mut", 0.0);
        const double omega0 = var_or(vars, "omega0", 0.0);
        const int i         = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "twist_decay: atom with tag 1 not found";
        const double r = lmp->atom->radius[i];
        expect_rel(omega0 - (5.0 * mut * g) / (2.0 * r * r) * t, lmp->atom->omega[i][2],
                   cfg.analytic_tol, "twist_decay omega_z");
    } else if (cfg.analytic_model == "twist_decay_marshall") {
        // like twist_decay, but for the marshall twisting model whose Coulomb
        // cap is derived from the tangential friction and the contact radius
        // a = sqrt(delta R_eff) (Marshall 2009, eq 44): M = (2/3) mu_t a N.
        // The contact radius is measured from the live overlap (R_eff = r for
        // a wall):  omega_z(t) = omega0 - (5 mu_t a g)/(3 r^2) t.
        const double g      = var_or(vars, "grav", 0.0);
        const double mu     = var_or(vars, "xmu", 0.0);
        const double omega0 = var_or(vars, "omega0", 0.0);
        const int i         = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "twist_decay_marshall: atom with tag 1 not found";
        const double r     = lmp->atom->radius[i];
        const double delta = r - lmp->atom->x[i][2];
        ASSERT_GT(delta, 0.0) << "twist_decay_marshall: sphere not in contact with the floor";
        const double a = std::sqrt(delta * r);
        expect_rel(omega0 - (5.0 * mu * a * g) / (3.0 * r * r) * t, lmp->atom->omega[i][2],
                   cfg.analytic_tol, "twist_decay_marshall omega_z");
    } else if (cfg.analytic_model == "heat_equilibration") {
        // two spheres (tags 1,2) held in static contact (no integrator) with
        // initial temperatures t1_0 and t2_0, coupled by granular heat
        // conduction and integrated by fix heat/flow (constant specific heat
        // cp).  The conductance H is h pi a^2 (area model, variable htc_area)
        // or 2 k a (radius model, variable htc_radius) with the contact radius
        // a = sqrt(delta R_eff) measured from the live overlap.  Then
        //   T1 - T2 = (t1_0 - t2_0) exp(-t H (1/(cp m1) + 1/(cp m2)))
        // and the mass-weighted mean temperature stays constant.
        const double t1_0 = var_or(vars, "t1_0", 0.0);
        const double t2_0 = var_or(vars, "t2_0", 0.0);
        const double cp   = var_or(vars, "cp", 0.0);
        const double ha   = var_or(vars, "htc_area", 0.0);
        const double hr   = var_or(vars, "htc_radius", 0.0);
        const int i1      = find_local(lmp, 1);
        const int i2      = find_local(lmp, 2);
        ASSERT_GE(i1, 0) << "heat_equilibration: atom with tag 1 not found";
        ASSERT_GE(i2, 0) << "heat_equilibration: atom with tag 2 not found";
        ASSERT_TRUE(lmp->atom->temperature_flag) << "heat_equilibration: no temperature property";
        const double m1 = lmp->atom->rmass[i1];
        const double m2 = lmp->atom->rmass[i2];
        const double r1 = lmp->atom->radius[i1];
        const double r2 = lmp->atom->radius[i2];
        const double dx = lmp->atom->x[i1][0] - lmp->atom->x[i2][0];
        const double dy = lmp->atom->x[i1][1] - lmp->atom->x[i2][1];
        const double dz = lmp->atom->x[i1][2] - lmp->atom->x[i2][2];
        const double delta = r1 + r2 - std::sqrt(dx * dx + dy * dy + dz * dz);
        ASSERT_GT(delta, 0.0) << "heat_equilibration: spheres not in contact";
        const double reff = r1 * r2 / (r1 + r2);
        const double a    = std::sqrt(delta * reff);
        const double hcond = (ha > 0.0) ? ha * MathConst::MY_PI * a * a : 2.0 * hr * a;
        const double rate  = hcond * (1.0 / (cp * m1) + 1.0 / (cp * m2));
        const double T1    = lmp->atom->temperature[i1];
        const double T2    = lmp->atom->temperature[i2];
        expect_rel((t1_0 - t2_0) * std::exp(-rate * t), T1 - T2, cfg.analytic_tol,
                   "heat_equilibration temperature difference");
        expect_rel((m1 * t1_0 + m2 * t2_0) / (m1 + m2), (m1 * T1 + m2 * T2) / (m1 + m2),
                   cfg.analytic_tol, "heat_equilibration mean temperature");
    } else if (cfg.analytic_model == "pulloff_jkr") {
        // JKR cohesion: at zero overlap the (tensile) contact force is
        // (8/9) of the pull-off force F_po = 3 pi gamma R_eff (the LAMMPS
        // cohesion convention, see the pair granular documentation), i.e.
        // |F(delta=0)| = (8/3) pi gamma R_eff.  The reference places the pair
        // at a vanishing overlap so this validates the pull-off force through
        // the exact 8/9 relation without having to time a detachment event.
        const double coh  = var_or(vars, "coh", 0.0);
        const double reff = var_or(vars, "reff", 0.0);
        const int i       = find_local(lmp, 1);
        ASSERT_GE(i, 0) << "pulloff_jkr: atom with tag 1 not found";
        const double fx = lmp->atom->f[i][0];
        const double fy = lmp->atom->f[i][1];
        const double fz = lmp->atom->f[i][2];
        const double fmag = std::sqrt(fx * fx + fy * fy + fz * fz);
        expect_rel((8.0 / 3.0) * MathConst::MY_PI * coh * reff, fmag, cfg.analytic_tol,
                   "pulloff_jkr force");
    } else {
        ADD_FAILURE() << "unknown analytic_model: '" << cfg.analytic_model << "'";
    }
}
