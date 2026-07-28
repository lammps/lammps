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
        elif model == 'hertz_normal_impact':
            f['hertz_normal_impact peak energy balance'] = (
                "(1/2)*m_red*vrela^2 = (2/5)*Pmax*alpha  [LHS here; RHS from sim force]",
                0.5 * v['mred_factor'] * massof(v, txt) * v['vrela'] ** 2)
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
