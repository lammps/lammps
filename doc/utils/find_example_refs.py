#!/usr/bin/env python3
"""find_example_refs.py -- Cross-reference LAMMPS doc pages with example
input scripts. See issue lammps/lammps#2680.

Modes:
  --report           Print a per-rst candidate report. Default mode.
  --apply FILE       Read a curated mapping and edit doc/src/*.rst in place,
                     inserting one "Example input scripts available: ..."
                     line right after the first Examples code block.
  --check            Verify every examples/ path mentioned in doc/src/*.rst
                     resolves to a real file or directory.

Mapping (curated) file format -- one entry per non-comment line:
    rst_basename: <path> [<path> ...]
  e.g.
    pair_charmm: examples/peptide/in.peptide
    compute_fep: examples/PACKAGES/fep

Conventions enforced:
  * Plain text path, no Sphinx :file: role, no URL.
  * Lists individual files when 1-3 are given; the curator switches to a
    folder when more would clutter the page.
  * The script is idempotent -- a page that already has an examples/
    reference (in any phrasing) is skipped.
"""

import argparse
import os
import re
import sys
from collections import defaultdict, Counter
from pathlib import Path


# Doc-index files that map style names to .rst page basenames.  Each line of
# interest looks like:  * :doc:`<style_name> [(suffix)] <rst_basename>`
INDEX_FILES = {
    'Commands_pair.rst':    'pair',
    'Commands_fix.rst':     'fix',
    'Commands_compute.rst': 'compute',
    'Commands_dump.rst':    'dump',
    'Commands_kspace.rst':  'kspace',
    'Commands_bond.rst':    'bond',     # multi-category; updated by heading
    'Commands_all.rst':     'command',
}

# Inside Commands_bond.rst, headings switch the active sub-category.
BOND_SECTIONS = {
    'Bond styles':     'bond',
    'Angle styles':    'angle',
    'Dihedral styles': 'dihedral',
    'Improper styles': 'improper',
}

# Token position of the style name for each style-declaring command.
# For 'fix', 'compute', 'dump' the style follows id+group, hence position 3.
DISPATCH = {
    'pair_style':     ('pair',     1),
    'bond_style':     ('bond',     1),
    'angle_style':    ('angle',    1),
    'dihedral_style': ('dihedral', 1),
    'improper_style': ('improper', 1),
    'kspace_style':   ('kspace',   1),
    'atom_style':     ('atom',     1),
    'fix':            ('fix',      3),
    'compute':        ('compute',  3),
    'dump':           ('dump',     3),
}

ACCEL_SUFFIXES = ('/gpu', '/intel', '/kk', '/omp', '/opt')

# Styles and commands that show up in nearly every input script and therefore
# would not point a reader to a meaningfully-illustrative example.  Adjust as
# needed during curation; the goal is to leave only entries where the example
# actually demonstrates the style/command.
TRIVIAL = {
    # Atom styles
    ('atom', 'atomic'), ('atom', 'molecular'), ('atom', 'charge'),
    ('atom', 'full'), ('atom', 'bond'), ('atom', 'hybrid'),

    # Pair "structural" styles
    ('pair', 'lj/cut'), ('pair', 'zero'), ('pair', 'none'),
    ('pair', 'hybrid'), ('pair', 'hybrid/overlay'), ('pair', 'hybrid/scaled'),
    ('pair', 'table'),

    # Bond/angle/dihedral/improper boilerplate
    ('bond', 'harmonic'), ('bond', 'none'), ('bond', 'zero'),
    ('bond', 'hybrid'),
    ('angle', 'harmonic'), ('angle', 'none'), ('angle', 'zero'),
    ('angle', 'hybrid'),
    ('dihedral', 'harmonic'), ('dihedral', 'none'), ('dihedral', 'zero'),
    ('dihedral', 'hybrid'),
    ('improper', 'harmonic'), ('improper', 'none'), ('improper', 'zero'),
    ('improper', 'hybrid'),

    # KSpace
    ('kspace', 'pppm'),

    # Fixes used everywhere
    ('fix', 'nve'), ('fix', 'nvt'), ('fix', 'npt'), ('fix', 'nph'),
    ('fix', 'langevin'), ('fix', 'viscous'), ('fix', 'temp/rescale'),
    ('fix', 'ave/time'), ('fix', 'ave/atom'), ('fix', 'ave/chunk'),
    ('fix', 'ave/correlate'), ('fix', 'ave/histo'),
    ('fix', 'recenter'), ('fix', 'momentum'), ('fix', 'balance'),
    ('fix', 'print'), ('fix', 'halt'),
    ('fix', 'store/state'), ('fix', 'store/force'),
    ('fix', 'addforce'), ('fix', 'setforce'),
    ('fix', 'property/atom'),

    # Computes used everywhere
    ('compute', 'temp'), ('compute', 'ke'), ('compute', 'pe'),
    ('compute', 'pressure'), ('compute', 'stress/atom'),
    ('compute', 'msd'), ('compute', 'rdf'),
    ('compute', 'displace/atom'),
    ('compute', 'com'), ('compute', 'gyration'),
    ('compute', 'chunk/atom'), ('compute', 'reduce'),
    ('compute', 'pair/local'), ('compute', 'property/atom'),

    # Dump styles used everywhere
    ('dump', 'atom'), ('dump', 'custom'), ('dump', 'xyz'), ('dump', 'dcd'),

    # Boilerplate top-level commands
    ('command', 'run'), ('command', 'read_data'), ('command', 'read_restart'),
    ('command', 'write_data'), ('command', 'write_restart'),
    ('command', 'write_dump'), ('command', 'restart'),
    ('command', 'thermo'), ('command', 'thermo_style'),
    ('command', 'thermo_modify'), ('command', 'dump_modify'),
    ('command', 'velocity'), ('command', 'neighbor'),
    ('command', 'neigh_modify'), ('command', 'timestep'),
    ('command', 'group'), ('command', 'region'), ('command', 'lattice'),
    ('command', 'mass'), ('command', 'boundary'), ('command', 'units'),
    ('command', 'dimension'), ('command', 'atom_modify'),
    ('command', 'comm_modify'), ('command', 'comm_style'),
    ('command', 'processors'), ('command', 'partition'),
    ('command', 'label'), ('command', 'jump'), ('command', 'next'),
    ('command', 'print'), ('command', 'shell'), ('command', 'suffix'),
    ('command', 'package'), ('command', 'log'), ('command', 'clear'),
    ('command', 'echo'), ('command', 'if'), ('command', 'quit'),
    ('command', 'variable'),
    ('command', 'kspace_modify'), ('command', 'pair_modify'),
    ('command', 'reset_timestep'), ('command', 'reset_atoms'),
    ('command', 'unfix'), ('command', 'uncompute'), ('command', 'undump'),
    ('command', 'fix_modify'), ('command', 'compute_modify'),
    ('command', 'min_modify'), ('command', 'min_style'),
    # Style-declaration verbs themselves (not the style values)
    ('command', 'pair_style'), ('command', 'bond_style'),
    ('command', 'angle_style'), ('command', 'dihedral_style'),
    ('command', 'improper_style'), ('command', 'kspace_style'),
    ('command', 'atom_style'),
    ('command', 'pair_coeff'), ('command', 'bond_coeff'),
    ('command', 'angle_coeff'), ('command', 'dihedral_coeff'),
    ('command', 'improper_coeff'),
}


DOC_ENTRY = re.compile(
    r'^\s*\*\s*:doc:`(?P<name>[^<]+?)\s*(?:\([^)]*\))?\s*'
    r'<(?P<page>[^>]+)>`'
)

EXAMPLES_REF = re.compile(r'examples/[A-Za-z0-9_./+-]+')

HEADING_UNDERLINE = re.compile(r'^[\"\-=^~+*#`\']+$')


def strip_accel(style):
    """Strip a trailing accelerator suffix from a style token, if any."""
    for s in ACCEL_SUFFIXES:
        if style.endswith(s):
            return style[:-len(s)]
    return style


def build_mapping(doc_dir):
    """Parse Commands_*.rst.

    Returns (mapping, pages):
      mapping  -- {(category, style_name): rst_basename}
      pages    -- {rst_basename: set of (category, style_name)}
    """
    mapping = {}
    pages = defaultdict(set)

    for fname, default_cat in INDEX_FILES.items():
        path = doc_dir / fname
        if not path.exists():
            continue
        category = default_cat
        prev_line = ''
        for raw in path.read_text(errors='ignore').splitlines():
            line = raw.rstrip()
            if fname == 'Commands_bond.rst' and line in BOND_SECTIONS:
                category = BOND_SECTIONS[line]
                prev_line = line
                continue
            m = DOC_ENTRY.match(line)
            if not m:
                prev_line = line
                continue
            name = m.group('name').strip()
            page = m.group('page').strip()
            # Skip cross-anchor links and external targets
            if '/' in page or '#' in page or '.' in page:
                prev_line = line
                continue
            # Drop accelerator-suffixed aliases (they redirect to base style)
            if any(name.endswith(s) for s in ACCEL_SUFFIXES):
                prev_line = line
                continue
            mapping[(category, name)] = page
            pages[page].add((category, name))
            prev_line = line
    return mapping, pages


def parse_input_script(path, mapping, seen=None, depth=0, max_depth=5):
    """Scan one input script.  Returns set of (category, style_name) tuples
    declared in the script (transitively through `include` directives)."""
    if seen is None:
        seen = set()
    path = Path(path).resolve()
    if path in seen or depth > max_depth or not path.exists():
        return set()
    seen.add(path)

    found = set()
    try:
        text = path.read_text(errors='ignore')
    except OSError:
        return found

    for raw in text.splitlines():
        line = raw.split('#', 1)[0].strip()
        if not line:
            continue
        toks = line.split()
        first = toks[0]

        if first == 'include' and len(toks) > 1:
            inc = path.parent / toks[1]
            if inc.exists():
                found |= parse_input_script(inc, mapping, seen,
                                            depth + 1, max_depth)
            continue

        if first in DISPATCH:
            cat, pos = DISPATCH[first]
            if pos < len(toks):
                style = toks[pos]
                if '$' in style:
                    continue
                base = strip_accel(style)
                if (cat, base) in mapping:
                    found.add((cat, base))
                # pair_style hybrid <sub> [args] <sub> [args] ...
                if first == 'pair_style' and base in ('hybrid',
                                                     'hybrid/overlay',
                                                     'hybrid/scaled'):
                    for tok in toks[2:]:
                        if '$' in tok:
                            continue
                        sub = strip_accel(tok)
                        if ('pair', sub) in mapping:
                            found.add(('pair', sub))
        elif ('command', first) in mapping:
            found.add(('command', first))
    return found


def find_input_scripts(examples_dir):
    """Yield every input-script path under examples_dir.  Matches in.* and
    in_*.lmp; ignores README, data files, log files, scripts, etc."""
    examples_dir = Path(examples_dir).resolve()
    for path in examples_dir.rglob('*'):
        if not path.is_file():
            continue
        name = path.name
        if name.startswith('in.') or (name.startswith('in_')
                                       and name.endswith('.lmp')):
            yield path


def scan_examples(examples_dir, mapping):
    """Return {rst_basename: [(script_path, line_count), ...]}."""
    refs = defaultdict(list)
    examples_dir = Path(examples_dir).resolve()
    for script in find_input_scripts(examples_dir):
        styles = parse_input_script(script, mapping)
        try:
            line_count = len(script.read_text(errors='ignore').splitlines())
        except OSError:
            line_count = 10**9
        rsts_added = set()
        for cat_name in styles:
            if cat_name in TRIVIAL:
                continue
            rst = mapping.get(cat_name)
            if not rst or rst in rsts_added:
                continue
            rsts_added.add(rst)
            refs[rst].append((script, line_count))
    return refs


def has_existing_ref(rst_path):
    """True if the .rst already mentions examples/ in any form."""
    try:
        text = Path(rst_path).read_text(errors='ignore')
    except OSError:
        return True
    return bool(EXAMPLES_REF.search(text))


def rank_candidates(candidates, repo_root):
    """Sort: non-PACKAGES first, shorter first, alphabetic last."""
    def key(c):
        path, lines = c
        rel = path.relative_to(repo_root)
        is_pkg = 1 if 'PACKAGES' in rel.parts else 0
        return (is_pkg, lines, str(rel))
    return sorted(candidates, key=key)


def select_paths(candidates, repo_root, max_files=3, dominance=0.75,
                 min_for_dominance=3):
    """Choose the path list for the proposed reference.

    Picks the dominant folder when 75%+ of candidates share a parent and that
    parent has at least 3 candidates; otherwise lists up to max_files
    individual scripts."""
    ranked = rank_candidates(candidates, repo_root)
    if not ranked:
        return []
    dirs = Counter()
    for path, _ in ranked:
        rel = path.relative_to(repo_root)
        dirs[str(rel.parent)] += 1
    best_dir, best_count = dirs.most_common(1)[0]
    if (best_count >= min_for_dominance
            and best_count / len(ranked) >= dominance):
        return [best_dir]
    return [str(p.relative_to(repo_root)) for p, _ in ranked[:max_files]]


def cmd_report(repo_root, output, list_cap=5):
    doc_dir = repo_root / 'doc' / 'src'
    ex_dir = repo_root / 'examples'

    mapping, pages = build_mapping(doc_dir)
    refs = scan_examples(ex_dir, mapping)

    out = []
    out.append('# Phase 2 candidate references for issue #2680.')
    out.append('# Generated by doc/utils/find_example_refs.py')
    out.append('# Format: rst_basename: <path> [<path> ...]')
    out.append('# Lines starting with # are comments. Edit, save, and pass')
    out.append('# this file to --apply.')
    out.append('')

    eligible = []
    skipped_existing = 0
    for rst in sorted(refs.keys()):
        rst_path = doc_dir / f'{rst}.rst'
        if not rst_path.exists():
            continue
        if has_existing_ref(rst_path):
            skipped_existing += 1
            continue
        eligible.append(rst)

    out.append(f'# Pages eligible: {len(eligible)}')
    out.append(f'# Pages already referencing examples/: {skipped_existing}')
    out.append('')

    for rst in eligible:
        cands = refs[rst]
        selected = select_paths(cands, repo_root)
        ranked = rank_candidates(cands, repo_root)[:list_cap]
        out.append(f'# --- {rst} ---')
        styles = sorted(f'{c}:{n}' for c, n in pages[rst])
        styles_line = ', '.join(styles)
        if len(styles_line) > 140:
            styles_line = styles_line[:137] + '...'
        out.append(f'# styles: {styles_line}')
        out.append(f'# candidates ({len(cands)} total):')
        for path, lines in ranked:
            rel = path.relative_to(repo_root)
            out.append(f'#   {rel}  [{lines}]')
        out.append(f'{rst}: {" ".join(selected)}')
        out.append('')

    text = '\n'.join(out) + '\n'
    if output:
        Path(output).write_text(text)
        sys.stderr.write(f'Wrote {len(eligible)} eligible pages '
                         f'to {output}\n')
    else:
        sys.stdout.write(text)


def insert_reference(rst_path, ref):
    """Insert the canonical reference line after the first Examples code
    block.  Returns True on success, False if no insertion point found."""
    with open(rst_path) as f:
        lines = f.readlines()

    ex_idx = None
    for i in range(len(lines) - 1):
        if lines[i].rstrip() in ('Examples', 'Example'):
            under = lines[i + 1].rstrip()
            if HEADING_UNDERLINE.match(under) and len(under) >= 4:
                ex_idx = i
                break
    if ex_idx is None:
        return False

    # Find first .. code-block:: in the Examples section.
    cb_idx = None
    for i in range(ex_idx + 2, len(lines)):
        stripped = lines[i].lstrip()
        if stripped.startswith('.. code-block::'):
            cb_idx = i
            break
        # If we hit another section heading first, give up.
        if (i + 1 < len(lines) and lines[i].rstrip() and
                HEADING_UNDERLINE.match(lines[i + 1].rstrip()) and
                len(lines[i + 1].rstrip()) >= 4):
            return False
    if cb_idx is None:
        return False

    # Walk past the directive options + body
    i = cb_idx + 1
    # Skip directive options (indented lines starting with ':')
    while i < len(lines) and lines[i].lstrip().startswith(':') \
            and lines[i].startswith((' ', '\t')):
        i += 1
    # Skip blank lines between directive and body
    while i < len(lines) and lines[i].strip() == '':
        i += 1
    # Walk through body (indented lines, plus blank lines inside body)
    while i < len(lines):
        line = lines[i]
        if line.strip() == '' or line.startswith((' ', '\t')):
            i += 1
            continue
        break
    # i now points to first non-indented non-blank line after the code block.
    # If the preceding line is blank, replace it; otherwise prepend a blank.
    if i > 0 and lines[i - 1].strip() == '':
        # Keep one blank above the insertion, one blank after.
        insertion = [f'Example input scripts available: {ref}\n', '\n']
    else:
        insertion = ['\n', f'Example input scripts available: {ref}\n', '\n']
    lines[i:i] = insertion

    with open(rst_path, 'w') as f:
        f.writelines(lines)
    return True


def cmd_apply(repo_root, mapping_file):
    doc_dir = repo_root / 'doc' / 'src'
    text = Path(mapping_file).read_text()

    applied = 0
    skipped = 0
    failed = 0
    bad_path = 0

    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        if ':' not in line:
            continue
        rst, paths_str = line.split(':', 1)
        rst = rst.strip()
        paths = paths_str.split()
        if not paths:
            continue
        rst_path = doc_dir / f'{rst}.rst'
        if not rst_path.exists():
            sys.stderr.write(f'WARN: {rst}.rst not found\n')
            failed += 1
            continue
        if has_existing_ref(rst_path):
            skipped += 1
            continue
        # Verify each cited example path exists.
        missing = [p for p in paths if not (repo_root / p).exists()]
        if missing:
            sys.stderr.write(f'WARN {rst}: missing paths {missing}\n')
            bad_path += 1
            continue
        ref = ', '.join(paths)
        if insert_reference(rst_path, ref):
            applied += 1
        else:
            sys.stderr.write(f'WARN {rst}: no Examples code-block found\n')
            failed += 1

    sys.stderr.write(f'Applied: {applied}  Skipped(existing): {skipped}  '
                     f'NoSection: {failed}  BadPath: {bad_path}\n')


def cmd_check(repo_root):
    doc_dir = repo_root / 'doc' / 'src'
    missing = []
    for rst in sorted(doc_dir.glob('*.rst')):
        text = rst.read_text(errors='ignore')
        for m in EXAMPLES_REF.finditer(text):
            p = m.group(0).rstrip('.,;:)')
            if not (repo_root / p).exists():
                missing.append((rst.name, p))
    for fname, p in missing:
        print(f'MISSING: {p}  (in {fname})')
    if missing:
        sys.exit(1)
    print('All examples/ paths in doc/src/*.rst resolve.')


def main():
    parser = argparse.ArgumentParser(
        description='Cross-reference doc pages with example input scripts '
                    '(LAMMPS issue #2680).')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--report', action='store_true',
                      help='Generate per-rst candidate report (default).')
    mode.add_argument('--apply', metavar='FILE',
                      help='Read curated mapping and edit doc/src/*.rst '
                           'in place.')
    mode.add_argument('--check', action='store_true',
                      help='Verify all examples/ paths in doc/src/*.rst '
                           'exist.')
    parser.add_argument('-o', '--output', metavar='FILE',
                        help='Write report to FILE (default stdout).')
    parser.add_argument('--repo', metavar='PATH',
                        help='Repo root (default: auto-detect from script '
                             'location).')
    parser.add_argument('--list-cap', type=int, default=5,
                        help='Max candidates shown per page in report '
                             '(default 5).')
    args = parser.parse_args()

    if args.repo:
        repo_root = Path(args.repo).resolve()
    else:
        repo_root = Path(__file__).resolve().parents[2]

    if args.check:
        cmd_check(repo_root)
    elif args.apply:
        cmd_apply(repo_root, args.apply)
    else:
        cmd_report(repo_root, args.output, args.list_cap)


if __name__ == '__main__':
    main()
