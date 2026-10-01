#!/usr/bin/env python3
# Analytic audit for unittest/granular. Per checked quantity, prints: the closed-form
# equation, the value it gives from the YAML variables (computed here, independent of
# the harness), the harness's own expected, the simulated value, and relative errors.
# Forces the harness to print by running a temp copy with a tiny tolerance. Leaves
# nothing behind. Run from the build dir:  python3 dem_audit.py [tests_dir]

import sys, os, re, math, glob, subprocess, tempfile
PI = math.pi
HERE = os.path.dirname(os.path.abspath(__file__))
# tests dir: 1st arg, else ./tests next to this script, else the cwd's tests
TESTS = os.path.expanduser(sys.argv[1]) if len(sys.argv) > 1 else (
    os.path.join(HERE, "tests") if os.path.isdir(os.path.join(HERE, "tests"))
    else os.path.join(os.getcwd(), "tests"))

def find_build():
    # 2nd arg wins; else search likely locations for the test_dem_01 binary
    if len(sys.argv) > 2:
        return os.path.expanduser(sys.argv[2])
    if os.getenv("DEM_BUILD"):
        return os.path.expanduser(os.getenv("DEM_BUILD"))
    cands = [os.getcwd()]
    # walk up from here looking for a sibling 'build' dir, then scan it
    up = HERE
    for _ in range(6):
        up = os.path.dirname(up)
        if not up or up == "/":
            break
        cands.append(os.path.join(up, "build"))
    for base in cands:
        if not base or not os.path.isdir(base):
            continue
        if os.path.exists(os.path.join(base, "test_dem_01")):
            return base
        hit = glob.glob(os.path.join(base, "**", "test_dem_01"), recursive=True)
        if hit:
            return os.path.dirname(hit[0])
    return os.getcwd()

BUILD = find_build()

def variables(txt):
    m = re.search(r'^variables:.*?\n(.*?)\n[A-Za-z_]\w*:', txt, re.S | re.M)
    if not m: return {}
    d = {}
    for ln in m.group(1).splitlines():
        mm = re.match(r'\s+([A-Za-z_]\w*)\s+(\S+)', ln)
        if mm:
            try: d[mm.group(1)] = float(mm.group(2))
            except ValueError: pass
    return d

def radius(v): return v.get('radius', v.get('diam', 0.0) / 2.0)
def dim_of(txt):
    m = re.search(r'\n\s+dimension\s+(\d)', txt)
    return int(m.group(1)) if m else 3
def massof(v, txt):
    d, dm = v.get('dens', 0.0), v.get('diam', 0.0)
    if (dim_of(txt) == 2):
        return d * (PI / 4.0) * dm ** 2
    else:
        return d * (PI / 6.0) * dm ** 3

def seg_time(txt):
    try:
        segs = [int(x) for x in re.search(r'^run_segments:.*?\n\s+([\d ]+)', txt, re.S | re.M).group(1).split()]
        aseg = int(re.search(r'^analytic_segment:\s*(-?\d+)', txt, re.M).group(1))
        dt = float(re.search(r'\n\s+dt\s+(\S+)', txt).group(1))
        if aseg < 0: aseg = len(segs) - 1
        return sum(segs[:aseg + 1]) * dt
    except Exception:
        return None


def derived(model, v, txt):
    r = radius(v); m = massof(v, txt); out = []
    out.append(("radius r", r))
    out.append(("mass m = dens*(pi/4)*diam^2 [2D disc]" if dim_of(txt)==2 else "mass m = dens*(pi/6)*diam^3", m))
    if model in ('collision_restitution',):
        out.append(("m_red = m/2 (two equal spheres)", m/2))
    if model == 'hertz_normal_impact':
        out.append(("m_red = mred_factor*m", v.get('mred_factor',1.0)*m))
        out.append(("(1/2)*m_red*vrela^2", 0.5*v.get('mred_factor',1.0)*m*v.get('vrela',0.0)**2))
    if model == 'rolling_decay':
        t = seg_time(txt)
        if t is not None: out.append(("t = sum(run_segments[:aseg+1])*dt", t))
    if model == 'pulloff_dmt':
        out.append(("R_eff = reff", v.get('reff')))
    return out

def formulas(model, v, txt):
    r = radius(v); f = {}
    try:
        if model == 'collision_restitution':
            f['collision_restitution e'] = ("e = -(v1x - v2x)/(2*vin)", v['en'])
        elif model == 'oblique_impact':
            e, mu, vz, vx = v['en'], v['xmu'], v['vz_in'], v['vx_in']
            f['oblique_impact vz_out']  = ("vz_out = e*vz_in", e * vz)
            f['oblique_impact vx_out']  = ("vx_out = vx_in - mu*(1+e)*vz_in", vx - mu * (1 + e) * vz)
            f['oblique_impact omega_y'] = ("omega_y = (5/2)*mu*(1+e)*vz_in/r", 2.5 * mu * (1 + e) * vz / r)
        elif model == 'oblique_impact_pair':
            e, mu, vn, vt = v['en'], v['xmu'], v['vn_in'], v['vt_in']
            dvt = mu * (1 + e) * vn
            f['oblique_impact_pair v1x'] = ("v1x = -e*vn_in", -e * vn)
            f['oblique_impact_pair v1y'] = ("v1y = vt_in - mu*(1+e)*vn_in", vt - dvt)
            f['oblique_impact_pair v2x'] = ("v2x = +e*vn_in", e * vn)
            f['oblique_impact_pair v2y'] = ("v2y = -(vt_in - mu*(1+e)*vn_in)", -(vt - dvt))
            f['oblique_impact_pair omega1z'] = ("omega1z = -(5/2)*mu*(1+e)*vn_in/r", -2.5 * dvt / r)
            f['oblique_impact_pair omega2z'] = ("omega2z = -(5/2)*mu*(1+e)*vn_in/r", -2.5 * dvt / r)
        elif model == 'slip_cessation':
            u0 = v['u0']
            f['slip_cessation vx']      = ("vx = 5*u0/7", 5 * u0 / 7)
            f['slip_cessation omega_y'] = ("omega_y = (5*u0/7)/r", (5 * u0 / 7) / r)
        elif model == 'spin_impact':
            e, mu, vin, w0 = v['en'], v['xmu'], v['vin'], v['omega0']
            f['spin_impact vz_out']  = ("vz_out = e*vin", e * vin)
            f['spin_impact vx_out']  = ("vx_out = mu*(1+e)*vin", mu * (1 + e) * vin)
            f['spin_impact omega_y'] = ("omega_y = omega0 - (5/2)*mu*(1+e)*vin/r", w0 - 2.5 * mu * (1 + e) * vin / r)
        elif model == 'rolling_decay':
            t = seg_time(txt)
            if t is not None:
                f['rolling_decay omega_y'] = (f"omega_y(t) = omega0 - (5*mur*g)/(2r)*t  [t={t:.4g}s]",
                                              v['omega0'] - (5 * v['mur'] * v['grav']) / (2 * r) * t)
            else:
                f['rolling_decay omega_y'] = ("omega_y(t) = omega0 - (5*mur*g)/(2r)*t", None)
        elif model == 'pulloff_dmt':
            f['pulloff_dmt force'] = ("F = 4*pi*coh*reff", 4 * PI * v['coh'] * v['reff'])
        elif model == 'spin_no_friction':
            f['spin_no_friction omega_y preserved'] = ("omega_y = omega0 (no tangential force)", v['omega0'])
        elif model == 'energy_dissipation':
            m = massof(v, txt)
            f['energy_dissipation initial energy'] = (
                "E_init = (1/2)*m*(vx_in^2+vz_in^2)  [final E must not exceed this]",
                0.5 * m * (v['vx_in'] ** 2 + v['vz_in'] ** 2))
        elif model == 'terminal_velocity_linear':
            m = massof(v, txt)
            f['terminal_velocity_linear'] = ("v_term = m*g/gamma", m * v['grav'] / v['gamma'])
        elif model == 'terminal_velocity_schiller_naumann':
            m, r = massof(v, txt), radius(v)
            rho, mu, g = v['rho_gas'], v['mu_gas'], v['grav']
            def drag(u):
                re = rho * u * (2.0 * r) / mu
                cd = (24.0 / re) * (1.0 + 0.15 * re ** 0.687)
                return 0.5 * cd * rho * PI * r * r * u * u
            lo, hi = 1e-12, 1.0
            while drag(hi) < m * g:
                hi *= 2.0
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if drag(mid) < m * g:
                    lo = mid
                else:
                    hi = mid
            f['terminal_velocity_schiller_naumann'] = (
                "v_term solves m*g = (1/2)*Cd(Re)*rho_g*pi*r^2*v^2 (Schiller-Naumann)",
                0.5 * (lo + hi))
        elif model == 'hertz_normal_impact':
            f['hertz_normal_impact peak energy balance'] = (
                "(1/2)*m_red*vrela^2 = (2/5)*Pmax*alpha  [LHS here; RHS from sim force]",
                0.5 * v['mred_factor'] * massof(v, txt) * v['vrela'] ** 2)
        elif model == 'hertz_peak':
            mred = v['mred_factor'] * massof(v, txt)
            alpha = (5.0 * mred * v['vrela'] ** 2 / (4.0 * v['kfac'])) ** 0.4
            f['hertz_peak alpha_max'] = ("alpha_max = (5*m_red*vrela^2/(4*kfac))^(2/5)", alpha)
            f['hertz_peak P_max'] = ("P_max = kfac*alpha_max^(3/2)", v['kfac'] * alpha ** 1.5)
        elif model == 'slip_transient':
            tt = seg_time(txt)
            if tt is not None:
                f['slip_transient u'] = (f"u(t) = u0 - mu*g*t  [t={tt:.4g}s]",
                                         v['u0'] - v['xmu'] * v['grav'] * tt)
                f['slip_transient omega_y'] = ("omega_y(t) = (5/2)*mu*g*t/r",
                                               2.5 * v['xmu'] * v['grav'] * tt / radius(v))
        elif model == 'incline_rolling':
            tt = seg_time(txt)
            a = (5.0 / 7.0) * v['grav'] * (v['sin_t'] - v['mur'] * v['cos_t'])
            if a > 0 and tt is not None:
                f['incline_rolling v'] = (f"v(t) = (5/7)*g*(sin-mur*cos)*t  [t={tt:.4g}s]", a * tt)
                f['incline_rolling omega_y'] = ("omega_y = v/r", a * tt / radius(v))
            elif a <= 0:
                f['incline_rolling v'] = ("mur >= tan(theta): stays at rest, v = 0", None)
        elif model == 'wall_restitution':
            f['wall_restitution vx_out'] = ("vx_out = -e*vx_in", -v['en'] * v['vx_in'])
        elif model == 'pulloff_jkr':
            f['pulloff_jkr force'] = ("|F(delta=0)| = (8/3)*pi*coh*reff = (8/9)*F_pulloff",
                                      (8.0 / 3.0) * PI * v['coh'] * v['reff'])
        elif model == 'twist_decay':
            tt = seg_time(txt)
            if tt is not None:
                r = radius(v)
                f['twist_decay omega_z'] = (
                    f"omega_z(t) = omega0 - (5*mut*g)/(2*r^2)*t  [t={tt:.4g}s]",
                    v['omega0'] - 5.0 * v['mut'] * v['grav'] / (2.0 * r * r) * tt)
        elif model == 'twist_decay_marshall':
            tt = seg_time(txt)
            if tt is not None:
                r = radius(v)
                a = math.sqrt((r - v['z0']) * r)
                f['twist_decay_marshall omega_z'] = (
                    f"omega_z(t) = omega0 - (5*xmu*a*g)/(3*r^2)*t  [a={a:.4g}, t={tt:.4g}s]",
                    v['omega0'] - 5.0 * v['xmu'] * a * v['grav'] / (3.0 * r * r) * tt)
        elif model == 'heat_equilibration':
            tt = seg_time(txt)
            if tt is not None:
                r = radius(v)
                delta = v['diam'] - 2.0 * v['sep']
                a = math.sqrt(delta * 0.5 * r)    # R_eff = r/2 for an equal pair
                hcond = v.get('htc_area', 0.0) * PI * a * a if 'htc_area' in v \
                    else 2.0 * v['htc_radius'] * a
                m = massof(v, txt)
                rate = 2.0 * hcond / (v['cp'] * m)
                f['heat_equilibration temperature difference'] = (
                    f"T1-T2 = dT0*exp(-t*H*2/(cp*m))  [H={hcond:.4g}, t={tt:.4g}s]",
                    (v['t1_0'] - v['t2_0']) * math.exp(-rate * tt))
                f['heat_equilibration mean temperature'] = (
                    "mean T constant (equal masses)", 0.5 * (v['t1_0'] + v['t2_0']))
        elif model == 'freefall':
            tt = seg_time(txt)
            if tt is not None:
                f['freefall z'] = (f"z(t) = z0 - g*t^2/2  [t={tt:.4g}s]",
                                   v['z0'] - 0.5 * v['grav'] * tt * tt)
                f['freefall vz'] = ("vz(t) = -g*t", -v['grav'] * tt)
        elif model == 'stack_energy':
            r = v['diam'] / 2.0
            m1 = v['dens1'] * (PI / 6.0) * v['diam'] ** 3
            m2 = v['dens2'] * (PI / 6.0) * v['diam'] ** 3
            df = max(0.0, r - (v['y1'] - v['ylo']))
            dc = max(0.0, (v['y2'] + r) - v['yhi'])
            dpp = max(0.0, 2.0 * r - (v['y2'] - v['y1']))
            pe0 = 0.5 * v['knorm'] * (df * df + dc * dc + dpp * dpp)
            f['stack_energy total energy'] = (
                "E0 = m1*g*y1 + m2*g*y2 + (1/2)*kn*(d_floor^2+d_pair^2+d_ceil^2)",
                m1 * v['grav'] * v['y1'] + m2 * v['grav'] * v['y2'] + pe0)
    except KeyError:
        pass
    return f

def harness_numbers(yaml_path, driver):
    txt = open(yaml_path).read()
    t = tempfile.NamedTemporaryFile('w', suffix='.yaml', delete=False)

    # use negative tolerance such that it is never satisfied to force comparison output
    t.write(re.sub(r'analytic_tol:\s*\S+', 'analytic_tol: -1', txt)); t.close()
    try:
        r = subprocess.run([driver, t.name, '-v'], capture_output=True, text=True, cwd=BUILD, timeout=300)
        out = r.stdout + r.stderr
    except Exception:
        out = ""
    finally:
        os.unlink(t.name)
    seen = {}
    for what, E, G in re.findall(r'([A-Za-z_][\w ]*?): expected (\S+) got (\S+)', out):
        try: seen[what.strip()] = (float(E), float(G))
        except ValueError: pass
    return seen

def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)

def main():
    print(f"tests: {TESTS}\nbuild: {BUILD}")
    files = sorted(glob.glob(os.path.join(TESTS, "dem*.yaml")))
    skipped = []
    ntest = 0
    nfail = 0
    for y in files:
        txt = open(y).read()
        mm = re.search(r'^analytic_model:\s*(\S+)', txt, re.M)
        nn = re.search(r'dem(\d\d)', os.path.basename(y))
        if not mm or not nn:
            skipped.append(os.path.basename(y)); continue
        model = mm.group(1)
        driver = os.path.join(BUILD, f"test_dem_{nn.group(1)}")
        v = variables(txt)
        fx = formulas(model, v, txt)
        hv = harness_numbers(y, driver) if os.path.exists(driver) else {}
        print(f"\n### {os.path.basename(y)}   [{model}]")
        if v:
            print("  variables: " + ", ".join(f"{k}={val:g}" for k, val in v.items()))
            dv = derived(model, v, txt)
            if dv:
                print("  derived:   " + ",  ".join(f"{name} = {val:.6g}" for name, val in dv))
        if not fx and not hv:
            print("     (no variables/quantities parsed)"); continue
        for k in (fx or hv):
            eq, mine = fx.get(k, ("(see test_analytic_models.cpp)", None))
            E, G = hv.get(k, (None, None))
            label = k.split(None, 1)[-1] if ' ' in k else k
            print(f"  {label:<26} {eq}")
            if E is not None:
                line = f"       harness expected = {E:.6g}   simulation = {G:.6g}   relerr = {rel(G, E):.2e}"
                if mine is not None:
                    ntest += 1
                    if rel(mine, E) <= 1e-4:
                        line += "   (formula = harness, OK)"
                    else:
                        line += f"   (formula = {mine:.6g} != harness), fail!!"
                        nfail += 1
                print(line)
            elif mine is not None:
                print(f"       formula value = {mine:.6g}   (harness value needs the driver)")
    if skipped:
        print("\nskipped (no analytic_model / not a demNN file):", ", ".join(skipped))

    print(f"Failed {nfail} out of {ntest} tests\n")
    print()

if __name__ == "__main__":
    main()
