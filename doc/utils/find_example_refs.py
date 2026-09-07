#!/usr/bin/env python3
"""find_example_refs.py -- Cross-reference LAMMPS doc pages with example
input scripts. See issue lammps/lammps#2680.

Modes:
  --report           Print a per-rst candidate report. Default mode.
  --apply FILE       Read a curated mapping and edit doc/src/*.rst in place,
                     inserting one "Example input scripts available: ..."
                     line right after the first Examples code block.
  --check            Verify every examples/ path mentioned in doc/src/*.rst
                     resolves to a real file or directory and is not nested
                     more deeply than examples/PACKAGES/<package>/<case>
                     (see MAX_DEPTH below).  Problems are reported as
                     "file:line:col: warning: ..." so the output
                     can be used from an editor's compilation mode (e.g.
                     Emacs M-x compile) to jump to the offending location.
  --count            Count in how many example input scripts each style or
                     command from the Commands_*.rst index pages is used and
                     list them with the most frequently used first.  Useful
                     for tuning the TRIVIAL skip list below.
  --update           Like --check, but try to find where a missing example
                     file or folder went (git rename history, then the same
                     or an equivalently spelled name elsewhere in examples/)
                     and rewrite the doc/src/*.rst references accordingly.
                     Ambiguous or unresolvable references are reported as
                     warnings for manual follow-up.

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
import glob
import os
import re
import subprocess
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
    ('atom', 'full'), ('atom', 'bond'),

    # Pair "structural" styles
    ('pair', 'lj/cut'), ('pair', 'zero'), ('pair', 'none'), ('pair', 'hybrid'), ('pair', 'table'),

    # Bond/angle/dihedral/improper boilerplate
    ('bond', 'harmonic'), ('bond', 'none'), ('bond', 'zero'),
    ('angle', 'harmonic'), ('angle', 'none'), ('angle', 'zero'),
    ('dihedral', 'harmonic'), ('dihedral', 'none'), ('dihedral', 'zero'),
    ('improper', 'harmonic'), ('improper', 'none'), ('improper', 'zero'),

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
    ('command', 'write_dump'), ('command', 'restart'), ('command', 'run_style'),
    ('command', 'thermo'), ('command', 'thermo_style'),
    ('command', 'thermo_modify'), ('command', 'dump_modify'),
    ('command', 'velocity'), ('command', 'neighbor'), ('command', 'dielectric'),
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
    ('command', 'variable'), ('command', 'special_bonds'),
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

# Editor backup and patch leftovers; never proposed as replacement targets.
JUNK_SUFFIXES = ('~', '.bak', '.orig', '.rej')

# Deepest folder nesting a documented example path may have, counted in
# path components including "examples" itself: examples/PACKAGES/<package>/
# <test-case> or examples/<package>/<test-case>.  Deeper paths make the
# documentation hard to read; the examples tree should be flattened instead,
# typically by combining several folder names into one longer name.
MAX_DEPTH = 4

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
    """True if the .rst already mentions an examples/ path."""
    try:
        text = Path(rst_path).read_text(errors='ignore')
    except OSError:
        return True
    return any(True for line in text.splitlines() for _ in refs_in_line(line))


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
    out.append('# Candidates for example references generated by doc/utils/find_example_refs.py')
    out.append('# Format: rst_basename: <path> [<path> ...]')
    out.append('# Lines starting with # are comments. Edit, save, and pass this file to --apply.')
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
        # wrap each path string in double backticks
        ref = ', '.join([f"``{p}``" for p in paths])

        if insert_reference(rst_path, ref):
            applied += 1
        else:
            sys.stderr.write(f'WARN {rst}: no Examples code-block found\n')
            failed += 1

    sys.stderr.write(f'Applied: {applied}  Skipped(existing): {skipped}  '
                     f'NoSection: {failed}  BadPath: {bad_path}\n')


def refs_in_line(line):
    """Yield (start, end, ref) for every examples/ path mentioned in one
    line of text.  start/end delimit the reference within the line
    (0-based, half-open) so that it can be replaced in place.  A trailing
    '*' is kept (glob pattern); trailing sentence punctuation is dropped.
    Two kinds of matches are not paths in the examples/ tree and are
    skipped: paths inside URLs (e.g. links to the examples section of the
    LAMMPS website) and the python/examples/ folder of the LAMMPS Python
    module, which holds Python scripts and Jupyter notebooks."""
    for m in EXAMPLES_REF.finditer(line):
        word = re.split(r'[\s`<]', line[:m.start()])[-1]
        if ('://' in word) or word.endswith('python/'):
            continue
        ref = m.group(0)
        end = m.end()
        if end < len(line) and line[end] == '*':
            ref += '*'
            end += 1
        else:
            stripped = ref.rstrip('.,;:')
            end -= len(ref) - len(stripped)
            ref = stripped
        yield m.start(), end, ref


def iter_example_refs(doc_dir):
    """Yield (rst_path, lineno, start, end, ref) for every examples/ path
    mentioned in doc/src/*.rst; lineno is 1-based, the rest is as returned
    by refs_in_line()."""
    for rst in sorted(doc_dir.glob('*.rst')):
        try:
            lines = rst.read_text(errors='ignore').splitlines()
        except OSError:
            continue
        for lineno, line in enumerate(lines, 1):
            for start, end, ref in refs_in_line(line):
                yield rst, lineno, start, end, ref


def ref_exists(repo_root, ref):
    """True if the referenced file, folder, or glob pattern exists."""
    if '*' in ref:
        return bool(glob.glob(str(repo_root / ref)))
    return (repo_root / ref).exists()


def excess_depth(repo_root, ref):
    """Return the number of folder levels of an existing examples/ reference
    if it exceeds MAX_DEPTH, else 0.  A file name or glob pattern at the
    end does not count as a folder level."""
    parts = ref.rstrip('/').split('/')
    if '*' in parts[-1] or (repo_root / ref).is_file():
        parts = parts[:-1]
    return len(parts) if len(parts) > MAX_DEPTH else 0


def location(rst, lineno, start):
    """Compiler-style 'file:line:col' message prefix.  The file name is
    given relative to the current directory, which is what editors resolve
    such messages against (e.g. Emacs compilation-mode)."""
    try:
        fname = os.path.relpath(rst)
    except ValueError:
        fname = str(rst)
    return f'{fname}:{lineno}:{start + 1}'


def cmd_check(repo_root):
    doc_dir = repo_root / 'doc' / 'src'
    problems = 0
    for rst, lineno, start, end, ref in iter_example_refs(doc_dir):
        loc = location(rst, lineno, start)
        if not ref_exists(repo_root, ref):
            problems += 1
            print(f'{loc}: warning: example path not found: {ref}')
        elif depth := excess_depth(repo_root, ref):
            problems += 1
            print(f'{loc}: warning: example path has {depth} folder levels '
                  f'(max {MAX_DEPTH}), the examples tree should be '
                  f'flattened: {ref}')
    if problems:
        print(f'{problems} problem(s) with examples/ paths in '
              'doc/src/*.rst.')
        sys.exit(1)
    print('All examples/ paths in doc/src/*.rst resolve.')


def cmd_count(repo_root, output):
    doc_dir = repo_root / 'doc' / 'src'
    ex_dir = repo_root / 'examples'

    mapping, _ = build_mapping(doc_dir)
    counts = Counter()
    nscripts = 0
    for script in find_input_scripts(ex_dir):
        nscripts += 1
        counts.update(parse_input_script(script, mapping))

    out = []
    out.append(f'# Number of example input scripts (out of {nscripts}) '
               f'using each of the {len(mapping)}')
    out.append('# styles/commands listed in Commands_*.rst; most frequently '
               'used first.')
    out.append('# T marks entries in the TRIVIAL skip list of '
               'find_example_refs.py.')
    out.append(f'# {"count":>5}  {"category":<9} {"style":<32} page')
    for key in sorted(mapping, key=lambda k: (-counts[k], k)):
        cat, name = key
        flag = '  T' if key in TRIVIAL else ''
        out.append(f'{counts[key]:7d}  {cat:<9} {name:<32} '
                   f'{mapping[key]}{flag}')

    text = '\n'.join(out) + '\n'
    if output:
        Path(output).write_text(text)
        sys.stderr.write(f'Wrote counts for {len(mapping)} styles/commands '
                         f'to {output}\n')
    else:
        sys.stdout.write(text)


def git_renames(repo_root):
    """Return {old_path: new_path} for every rename git recorded below
    examples/; the most recent rename wins for each old path.  Empty when
    git or the history is not available."""
    cmd = ['git', '-C', str(repo_root), 'log', '--diff-filter=R',
           '--name-status', '--format=', '-M', '--', 'examples/']
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              check=True)
    except (OSError, subprocess.CalledProcessError):
        return {}
    renames = {}
    for line in proc.stdout.splitlines():
        parts = line.split('\t')
        if len(parts) == 3 and parts[0].startswith('R'):
            renames.setdefault(parts[1], parts[2])
    return renames


def normalized_name(name):
    """Fold the spelling variants seen in example renames: letter case,
    the separators '-', '_', '.', and a trailing '.lmp' extension."""
    if name.endswith('.lmp'):
        name = name[:-4]
    return re.sub(r'[-_.]+', '.', name.lower())


class ExampleTree:
    """Index of the examples/ tree plus git rename history, used to find
    the new location of a moved or renamed example file or folder."""

    def __init__(self, repo_root):
        self.root = repo_root
        self.by_name = defaultdict(list)
        self.by_norm = defaultdict(list)
        for dirpath, dirnames, filenames in os.walk(repo_root / 'examples'):
            rel = Path(dirpath).relative_to(repo_root)
            for name in dirnames + filenames:
                if name.endswith(JUNK_SUFFIXES):
                    continue
                self.by_name[name].append(str(rel / name))
                self.by_norm[normalized_name(name)].append(str(rel / name))
        self.renames = git_renames(repo_root)

    def resolve(self, ref):
        """Return (new_ref, candidates).  new_ref is the unique replacement
        for the missing reference or None; candidates lists the alternatives
        when the search was ambiguous (empty when nothing was found)."""
        trailing = '/' if ref.endswith('/') else ''
        path = ref.rstrip('/')
        for finder in (self._via_git, self._via_name):
            cands = finder(path)
            if len(cands) == 1:
                return cands[0] + trailing, []
            if cands:
                return None, cands
        return None, []

    def _via_git(self, path, depth=0):
        """Follow git rename records: directly for a file, or for a folder
        via the files that were renamed below it."""
        if depth > 5 or not self.renames:
            return []
        seen = set()
        p = path
        while p in self.renames and p not in seen:
            seen.add(p)
            p = self.renames[p]
        if p != path and (self.root / p).exists():
            return [p]
        prefix = path + '/'
        cands = Counter()
        for old, new in self.renames.items():
            if old.startswith(prefix):
                rest = old[len(prefix):]
                if new.endswith('/' + rest):
                    cands[new[:-len(rest) - 1]] += 1
        found = []
        for cand, _ in cands.most_common():
            if (self.root / cand).is_dir():
                found.append(cand)
            else:
                found.extend(self._via_git(cand, depth + 1))
        return list(dict.fromkeys(found))

    def _via_name(self, path):
        """Same file or folder name elsewhere in the examples/ tree; failing
        that, one that differs only in separators, case, or a .lmp suffix
        (e.g. in.bpm.fracture -> in.bpm-fracture)."""
        name = os.path.basename(path)
        cands = self.by_name.get(name)
        if not cands:
            cands = self.by_norm.get(normalized_name(name), [])
        return self._closest(path, cands)

    @staticmethod
    def _closest(path, cands):
        """Keep the candidates whose parent folders share the longest
        trailing run of folder names with the old path."""
        if len(cands) < 2:
            return cands
        old_parts = path.split('/')[:-1]

        def score(cand):
            parts = cand.split('/')[:-1]
            n = 0
            while (n < min(len(parts), len(old_parts))
                   and parts[-1 - n] == old_parts[-1 - n]):
                n += 1
            return n

        best = max(score(c) for c in cands)
        return [c for c in cands if score(c) == best]


def cmd_update(repo_root):
    doc_dir = repo_root / 'doc' / 'src'
    tree = ExampleTree(repo_root)
    if not tree.renames:
        sys.stderr.write('NOTE: git rename history not available, using '
                         'file name heuristics only\n')

    edits = defaultdict(list)
    updated = 0
    unresolved = 0
    for rst, lineno, start, end, ref in iter_example_refs(doc_dir):
        loc = location(rst, lineno, start)
        if ref_exists(repo_root, ref):
            if depth := excess_depth(repo_root, ref):
                unresolved += 1
                print(f'{loc}: warning: example path has {depth} folder '
                      f'levels (max {MAX_DEPTH}), the examples tree should '
                      f'be flattened: {ref}')
            continue
        if '*' in ref:
            print(f'{loc}: warning: no files match pattern: {ref}')
            unresolved += 1
            continue
        new, cands = tree.resolve(ref)
        if new and (depth := excess_depth(repo_root, new)):
            unresolved += 1
            print(f'{loc}: warning: {ref} moved to {new}, which has {depth} '
                  f'folder levels (max {MAX_DEPTH}); not updated, the '
                  'examples tree should be flattened instead')
        elif new:
            edits[rst].append((lineno, start, end, new))
            updated += 1
            print(f'{loc}: note: {ref} -> {new}')
        else:
            unresolved += 1
            hint = ''
            if cands:
                hint = '; candidates: ' + ', '.join(cands[:5])
                if len(cands) > 5:
                    hint += f' (+{len(cands) - 5} more)'
            print(f'{loc}: warning: example path not found: {ref}{hint}')

    for rst, items in edits.items():
        with open(rst, newline='') as f:
            lines = f.readlines()
        # Replace from the end so earlier offsets stay valid.
        for lineno, start, end, new in sorted(items, reverse=True):
            line = lines[lineno - 1]
            lines[lineno - 1] = line[:start] + new + line[end:]
        with open(rst, 'w', newline='') as f:
            f.writelines(lines)

    print(f'Updated: {updated} reference(s) in {len(edits)} file(s)  '
          f'Unresolved: {unresolved}')
    if unresolved:
        sys.exit(1)


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
                           'exist; report problems compiler-style as '
                           'file:line:col: warning: ...')
    mode.add_argument('--count', action='store_true',
                      help='Count the example input scripts using each '
                           'style/command; most frequently used first.')
    mode.add_argument('--update', action='store_true',
                      help='Find where missing examples/ paths in '
                           'doc/src/*.rst were moved or renamed to and '
                           'rewrite the references in place.')
    parser.add_argument('-o', '--output', metavar='FILE',
                        help='Write report or counts to FILE '
                             '(default stdout).')
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
    elif args.count:
        cmd_count(repo_root, args.output)
    elif args.update:
        cmd_update(repo_root)
    elif args.apply:
        cmd_apply(repo_root, args.apply)
    else:
        cmd_report(repo_root, args.output, args.list_cap)


if __name__ == '__main__':
    main()
