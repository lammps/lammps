#!/usr/bin/env python3
"""Accelerator coverage scan for LAMMPS pair styles.

Lists every base pair style with the accelerator variants it has (OPENMP, GPU,
KOKKOS, INTEL, OPT), how many force-style YAML references exercise it, how many
example and benchmark inputs use it, the size of its implementation, and a few
structural signals that predict porting effort.  The output is what the tables in
accelerator-porting-backlog.md are derived from; rerun it instead of trusting them.

Usage (from anywhere; the repository root is found relative to this file):
  accelerator-coverage-scan.py                 # full table, tab separated
  accelerator-coverage-scan.py --summary       # totals only
  accelerator-coverage-scan.py --missing omp   # styles without an /omp variant
  accelerator-coverage-scan.py --missing gpu --pkg EXTRA-PAIR --yaml
  accelerator-coverage-scan.py --root /path/to/lammps

Signals (comma separated in the last column):
  comm    forward_comm/reverse_comm in compute()   (per-atom intermediates)
  map     atom->map() lookups                       (topology or partner search)
  full    requests a full neighbor list             (many-body or two-sided)
  ctab    Coulomb table lookup (ncoultablebits)     (no tables on the GPU)
  nbody   ev_tally3/ev_tally4/v_tally*              (three-body or custom tally)
  bonus   ellipsoid/line/tri/body bonus data
  fix     creates or looks up a fix
  rng     random number generator
  torque  writes torque or reads dipoles
  rmass   per-atom mass (finite-size particles)
  single  implements single()

Only the Python standard library is used; safe to run with python3 -I.
"""

import argparse
import collections
import glob
import os
import re
import sys

SUFFIXES = ("omp", "gpu", "kk", "intel", "opt")
FLAGS = "OGKIP"
REG = re.compile(r"^\s*PairStyle\(\s*([^,\s]+)\s*,\s*([^)]+)\)", re.M)


def find_root(start):
    here = os.path.abspath(start)
    while here != os.path.dirname(here):
        if os.path.isdir(os.path.join(here, "src")) and os.path.isdir(
            os.path.join(here, "cmake")
        ):
            return here
        here = os.path.dirname(here)
    sys.exit("cannot find the LAMMPS repository root; use --root")


def read(path):
    with open(path, encoding="utf-8", errors="replace") as fp:
        return fp.read()


def collect_styles(root):
    """Return {base_name: header_path} and {base_name: {suffix: header_path}}."""
    styles = {}
    for path in glob.glob(os.path.join(root, "src", "**", "*.h"), recursive=True):
        if os.sep + "STUBS" + os.sep in path:
            continue
        for m in REG.finditer(read(path)):
            styles[m.group(1)] = path
    base, variants = {}, collections.defaultdict(dict)
    for name, path in styles.items():
        parts = name.split("/")
        if len(parts) > 2 and parts[-2] == "kk" and parts[-1] in ("device", "host"):
            continue
        if len(parts) > 1 and parts[-1] in SUFFIXES:
            variants["/".join(parts[:-1])][parts[-1]] = path
        else:
            base[name] = path
    return base, variants


def collect_yaml(root, base):
    """Count YAML references per base style; hybrid sub-styles count for each."""
    counts = collections.Counter()
    pattern = os.path.join(root, "unittest", "force-styles", "tests", "*pair*.yaml")
    for path in glob.glob(pattern):
        for line in read(path).splitlines():
            if line.startswith("pair_style:"):
                words = line.split(":", 1)[1].split()
                if not words:
                    break
                if words[0].startswith("hybrid"):
                    for word in words[1:]:
                        if word in base:
                            counts[word] += 1
                else:
                    counts[words[0]] += 1
                break
    return counts


def collect_usage(root, base):
    """Count example and bench inputs that use each base style (once per file)."""
    counts = collections.Counter()
    for sub in ("examples", "bench"):
        for path in glob.glob(os.path.join(root, sub, "**", "in.*"), recursive=True):
            seen = set()
            for line in read(path).splitlines():
                words = line.split()
                if words and words[0] == "pair_style":
                    for word in words[1:]:
                        for suffix in SUFFIXES:
                            word = word.replace("/" + suffix, "")
                        if word in base:
                            seen.add(word)
            for word in seen:
                counts[word] += 1
    return counts


def signals(header):
    source = header[:-2] + ".cpp"
    if not os.path.exists(source):
        return 0, "header-only"
    text = read(source)
    found = []
    checks = (
        ("comm", ("forward_comm", "reverse_comm")),
        ("map", ("atom->map(",)),
        ("full", ("REQ_FULL",)),
        ("ctab", ("ncoultablebits",)),
        ("nbody", ("ev_tally3", "ev_tally4", "v_tally")),
        ("bonus", ("bonus",)),
        ("fix", ("modify->add_fix", "modify->get_fix", "fix_")),
        ("rng", ("RanMars", "RanPark", "uniform()", "gaussian()")),
        ("torque", ("torque", "mu[")),
        ("rmass", ("rmass",)),
        ("single", ("::single(",)),
    )
    for name, needles in checks:
        if any(needle in text for needle in needles):
            found.append(name)
    return text.count("\n"), ",".join(found)


def package(root, header):
    rel = os.path.relpath(header, os.path.join(root, "src"))
    parts = rel.split(os.sep)
    return parts[0] if len(parts) > 1 else "core"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", help="LAMMPS repository root (default: auto)")
    ap.add_argument("--summary", action="store_true", help="print totals only")
    ap.add_argument(
        "--missing",
        choices=SUFFIXES,
        help="list only styles WITHOUT this accelerator variant",
    )
    ap.add_argument("--pkg", help="restrict to one package directory name")
    ap.add_argument(
        "--yaml", action="store_true", help="restrict to styles with a YAML reference"
    )
    args = ap.parse_args()

    root = args.root or find_root(os.path.dirname(__file__))
    base, variants = collect_styles(root)
    yaml = collect_yaml(root, base)
    usage = collect_usage(root, base)

    rows = []
    for name in sorted(base):
        header = base[name]
        have = variants.get(name, {})
        if args.missing and args.missing in have:
            continue
        pkg = package(root, header)
        if args.pkg and pkg != args.pkg:
            continue
        if args.yaml and not yaml[name]:
            continue
        flags = "".join(f if s in have else "-" for f, s in zip(FLAGS, SUFFIXES))
        lines, sig = signals(header)
        rows.append((name, pkg, flags, yaml[name], usage[name], lines, sig))

    if args.summary or not rows:
        total = len(base)
        print("base pair styles:", total)
        for f, s in zip(FLAGS, SUFFIXES):
            n = sum(1 for b in base if s in variants.get(b, {}))
            print("  with /%s: %d" % (s, n))
        print("  with a YAML reference:", sum(1 for b in base if yaml[b]))
        print("  without any variant:", sum(1 for b in base if not variants.get(b)))
        if args.summary:
            return
    print("name\tpkg\t%s\tyaml\tuse\tlines\tsignals" % FLAGS)
    for row in rows:
        print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    main()
