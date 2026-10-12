#!/usr/bin/env python3
# Validate the style/parameter tables in doc/src/fix_adapt.rst and
# doc/src/fix_adapt_fep.rst against the extract() implementations in src/.
#
# For every table row this checks that:
#   - the style name is a registered style (PairStyle/BondStyle/... macro)
#   - the style class (or a base class) implements extract()
#   - every listed parameter is an extract() keyword (exact, case-sensitive)
#   - the keyword's dim is compatible with the fix:
#       fix adapt:     pair dim 0 or 2; bond/angle/dihedral/improper dim 1
#       fix adapt/fep: pair dim 2 only
#   - (pair tables) the scope column matches the dim (type pairs <-> dim 2)
#
# It then reports the reverse direction: extract() keywords with a
# compatible dim and double type that are NOT listed in the tables,
# skipping accelerator (/gpu /intel /kk /omp /opt), internal (ALL-CAPS),
# hybrid, and deprecated styles.  Styles matching a pattern in
# SKIP_STYLES (or given via --skip) are also left out of that report;
# table rows that DO list such a style are still validated.
#
# Messages use the GNU compiler format "file:line: level: text" so
# tools parsing compiler output (e.g. the emacs compilation mode) can
# jump to the location: doc table problems point to the table row in
# the .rst file, missing-entry notes to the extract() implementation.

import argparse
import difflib
import fnmatch
import os
import re
import sys
from collections import defaultdict

KINDS = ('pair', 'bond', 'angle', 'dihedral', 'improper')
BASE_CLASSES = {'Pair', 'Bond', 'Angle', 'Dihedral', 'Improper'}
VALID_DIMS = {
    'fix_adapt.rst': {'pair': {0, 2}, 'bond': {1}, 'angle': {1},
                      'dihedral': {1}, 'improper': {1}},
    'fix_adapt_fep.rst': {'pair': {2}},
}
SUFFIXES = ('/gpu', '/intel', '/kk', '/kk/device', '/kk/host', '/omp', '/opt')
# fnmatch() glob patterns of style names to leave out of the
# missing-entry report (any style kind); extend ad hoc with --skip
SKIP_STYLES = (
    'oxdna*',                             # CG-DNA package
    'oxrna*',                             # CG-DNA package
    'buck6d*',                            # MOFFF package
    'amoeba*',                            # AMOEBA package
    'thole*',                             # DRUDE package
    '*/cs',                               # CORESHELL package
    '*/soft',                             # FEP package
    'zero'
)
NUMERIC_TYPES = ('double', 'float')
OTHER_TYPES = ('int', 'bigint', 'tagint', 'imageint', 'bool', 'char')

STYLE_RE = re.compile(r'^\s*(Pair|Bond|Angle|Dihedral|Improper)Style\(\s*([^,\s)]+)\s*,\s*(\w+)', re.M)
CLASS_RE = re.compile(r'\bclass\s+(\w+)\s*(?::\s*(?:public|protected|private)\s+(\w+))?[^;{]*\{')
EXTRACT_RE = re.compile(r'\bvoid\s*\*\s*(\w+)::extract\s*\(')
DIM_RE = re.compile(r'\bdim\s*=\s*(\d+)')
STRCMP_RE = re.compile(r'strcmp\(\s*str\s*,\s*"([^"]+)"\s*\)\s*[!=]=\s*0')
RETURN_ID_RE = re.compile(r'return\s*\(\s*void\s*\*\s*\)\s*&?\s*([A-Za-z_]\w*)')

SECTION_RE = re.compile(r'^The \*(\w+)\* keyword')
ROW_RE = re.compile(r'^\|\s*:doc:`([^<`]+)<[^>]*>`\s*\|([^|]*)\|([^|]*)\|')


def walk_files(top, ext):
    for root, dirs, files in os.walk(top):
        dirs[:] = [d for d in dirs if d not in ('STUBS', 'MAKE')]
        for f in sorted(files):
            if f.endswith(ext):
                yield os.path.join(root, f)


def scan_source(srcdir):
    """Collect style registrations, class hierarchy, and extract() keywords."""
    styles = {k: {} for k in KINDS}       # kind -> style name -> class
    parent = {}                           # class -> base class
    class_header = {}                     # class -> header path
    for path in walk_files(srcdir, '.h'):
        try:
            text = open(path, encoding='utf-8', errors='replace').read()
        except OSError:
            continue
        for m in STYLE_RE.finditer(text):
            kind, name, cls = m.group(1).lower(), m.group(2), m.group(3)
            styles[kind].setdefault(name, cls)
        for m in CLASS_RE.finditer(text):
            cls, base = m.group(1), m.group(2)
            if cls not in class_header:
                class_header[cls] = path
                if base:
                    parent[cls] = base

    extracts = {}                         # class -> (keywords, file, line)
    for path in walk_files(srcdir, '.cpp'):
        try:
            text = open(path, encoding='utf-8', errors='replace').read()
        except OSError:
            continue
        for m in EXTRACT_RE.finditer(text):
            cls = m.group(1)
            body = function_body(text, m.end())
            if body is not None:
                line = text.count('\n', 0, m.start()) + 1
                extracts[cls] = (parse_extract_body(body), path, line)
    return styles, parent, class_header, extracts


def function_body(text, pos):
    """Return the text of the brace-delimited body starting at/after pos."""
    start = text.find('{', pos)
    if start < 0:
        return None
    depth = 0
    for i in range(start, len(text)):
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0:
                return text[start:i + 1]
    return None


def parse_extract_body(body):
    """Map extract() keywords to (dim, returned identifier).

    Two 'dim = N' placements occur in the sources: a running assignment
    before a group of strcmp() tests, or an assignment inside the matched
    if-block between the strcmp() and its return.  For each keyword the
    in-block dim wins; otherwise the closest preceding top-level dim
    applies.  Keywords with neither (e.g. pair_table's early-return
    style) get the first dim assigned anywhere in the body.
    """
    kws = []                              # (start, end-of-return, keyword, ident)
    for m in STRCMP_RE.finditer(body):
        rpos = body.find('return', m.end())
        end = body.find(';', rpos) if rpos >= 0 else -1
        end = end if end >= 0 else m.end()
        rm = RETURN_ID_RE.search(body, m.end(), end + 1)
        kws.append((m.start(), end, m.group(1), rm.group(1) if rm else None))
    dims = [(m.start(), int(m.group(1))) for m in DIM_RE.finditer(body)]
    toplevel = [(p, d) for p, d in dims
                if not any(s < p <= e for s, e, _, _ in kws)]
    keywords = {}
    for start, end, kw, ident in kws:
        inblock = [d for p, d in dims if start < p <= end]
        before = [d for p, d in toplevel if p < start]
        if inblock:
            dim = inblock[-1]
        elif before:
            dim = before[-1]
        else:
            dim = dims[0][1] if dims else None
        keywords.setdefault(kw, (dim, ident))
    return keywords


def class_chain(cls, parent):
    chain = []
    while cls and cls not in BASE_CLASSES and cls not in chain:
        chain.append(cls)
        cls = parent.get(cls)
    return chain


def find_extract(cls, parent, extracts):
    """Return (implementing class, keywords, (file, line)) or (None, {}, None)."""
    for c in class_chain(cls, parent):
        if c in extracts:
            keywords, path, line = extracts[c]
            return c, keywords, (path, line)
    return None, {}, None


def member_type(ident, cls, parent, class_header, cache={}):
    """Best-effort C++ type of a class member; None if not found."""
    if ident is None:
        return None
    key = (cls, ident)
    if key not in cache:
        result = None
        decl = re.compile(r'^\s*(?:static\s+|const\s+)*(%s)\b[^;()=]*\b%s\b[^;()]*;'
                          % ('|'.join(NUMERIC_TYPES + OTHER_TYPES), re.escape(ident)), re.M)
        for c in class_chain(cls, parent):
            hdr = class_header.get(c)
            if not hdr:
                continue
            m = decl.search(open(hdr, encoding='utf-8', errors='replace').read())
            if m:
                result = m.group(1)
                break
        cache[key] = result
    return cache[key]


def parse_doc_tables(path):
    """Yield (kind, lineno, [styles], [params], scope) for each table row."""
    section = None
    rows = []
    for lineno, line in enumerate(open(path, encoding='utf-8', errors='replace'), 1):
        m = SECTION_RE.match(line)
        if m:
            section = m.group(1) if m.group(1) in KINDS else section
            continue
        m = ROW_RE.match(line)
        if m and section in KINDS:
            names = [s.strip() for s in m.group(1).split(',') if s.strip()]
            params = [p.strip() for p in m.group(2).split(',') if p.strip()]
            rows.append((section, lineno, names, params, m.group(3).strip()))
    return rows


def check_doc_file(rst, rows, styles, parent, extracts, class_header, report):
    valid_dims = VALID_DIMS[os.path.basename(rst)]
    documented = defaultdict(set)         # (kind, style) -> set(params)

    for kind, lineno, names, params, scope in rows:
        for name in names:
            documented[(kind, name)].update(params)
            loc = (rst, lineno)
            cls = styles[kind].get(name)
            if cls is None:
                sugg = suggest(name, styles[kind])
                report('error:', loc, '%s style "%s" is not a registered style%s'
                       % (kind, name, sugg))
                continue
            impl, keywords, _ = find_extract(cls, parent, extracts)
            if impl is None:
                report('error:', loc, '%s style "%s" (%s) has no extract() method '
                       'but is listed with: %s' % (kind, name, cls, ','.join(params)))
                continue
            for p in params:
                if p not in keywords:
                    sugg = suggest(p, keywords)
                    report('error:', loc, '%s style "%s": param "%s" is not an extract() '
                           'keyword of %s%s' % (kind, name, p, impl, sugg))
                    continue
                dim, ident = keywords[p]
                mtype = member_type(ident, cls, parent, class_header)
                if mtype in OTHER_TYPES:
                    report('warning:', loc, '%s style "%s": param "%s" is backed by '
                           'variable "%s" of type %s; fix adapt writes double '
                           'values through this pointer' % (kind, name, p, ident, mtype))
                if dim not in valid_dims[kind]:
                    report('error:', loc, '%s style "%s": param "%s" has dim %s; '
                           '%s supports dim %s for %s styles'
                           % (kind, name, p, dim,
                              'fix adapt/fep' if 'fep' in rst else 'fix adapt',
                              '/'.join(map(str, sorted(valid_dims[kind]))), kind))
                elif kind == 'pair':
                    is_global = 'global' in scope
                    if dim == 2 and is_global:
                        report('warning:', loc, 'pair style "%s": param "%s" is per-type '
                               '(dim 2) but scope column says "%s"' % (name, p, scope))
                    if dim == 0 and not is_global:
                        report('warning:', loc, 'pair style "%s": param "%s" is global '
                               '(dim 0) but scope column says "%s"' % (name, p, scope))
    return documented


def report_missing(rst, documented, styles, parent, extracts, class_header,
                   skip, report):
    valid_dims = VALID_DIMS[os.path.basename(rst)]
    groups = defaultdict(list)            # (kind, impl, missing) -> [styles]
    skipped = []
    for kind in valid_dims:
        for name, cls in sorted(styles[kind].items()):
            if (name.endswith(SUFFIXES) or re.fullmatch(r'[A-Z0-9_/]+', name)
                    or name.startswith('hybrid') or 'Deprecated' in cls):
                continue
            if any(fnmatch.fnmatchcase(name, pat) for pat in skip):
                skipped.append(name)
                continue
            impl, keywords, loc = find_extract(cls, parent, extracts)
            if impl is None:
                continue
            adaptable = set()
            for kw, (dim, ident) in keywords.items():
                mtype = member_type(ident, cls, parent, class_header)
                if dim in valid_dims[kind] and mtype in NUMERIC_TYPES + (None,):
                    adaptable.add((kw, dim))
            missing = {(kw, dim) for kw, dim in adaptable
                       if kw not in documented.get((kind, name), set())}
            if missing:
                groups[(kind, impl, loc, frozenset(missing))].append(name)
    for (kind, impl, loc, missing), names in sorted(groups.items(), key=lambda g: g[1]):
        kws = ', '.join('%s (dim %d)' % (kw, dim) for kw, dim in sorted(missing))
        report('note:', loc, '%s table in %s has no entry for style%s %s (%s): %s'
               % (kind, os.path.basename(rst), 's' if len(names) > 1 else '',
                  ', '.join(names), impl, kws))
    if skipped:
        print('%d styles matching skip patterns (%s) were not checked '
              'for missing entries' % (len(skipped), ', '.join(skip)))


def suggest(word, candidates):
    folded = [c for c in candidates if c.lower() == word.lower()]
    match = folded or difflib.get_close_matches(word, list(candidates), 1, 0.5)
    return '; did you mean "%s"?' % match[0] if match else ''


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('--root', help='LAMMPS repository root (default: script location)', required=True)
    argp.add_argument('--no-missing', action='store_true',
                      help='skip the code -> doc missing-entry report')
    argp.add_argument('--skip', action='append', default=[], metavar='PATTERN',
                      help='add a style name glob pattern to leave out of the '
                      'missing-entry report (repeatable); table rows listing '
                      'a matching style are still validated')
    args = argp.parse_args()
    skip = list(SKIP_STYLES) + args.skip

    srcdir = os.path.join(args.root, 'src')
    docdir = os.path.join(args.root, 'doc', 'src')
    if not (os.path.isdir(srcdir) and os.path.isdir(docdir)):
        sys.exit('cannot find src/ and doc/src/ under ' + args.root)

    print('scanning %s ...' % srcdir)
    styles, parent, class_header, extracts = scan_source(srcdir)
    print('found %s registered styles, %d extract() implementations'
          % ('/'.join(str(len(styles[k])) for k in KINDS), len(extracts)))

    counts = defaultdict(int)
    for rst in ('fix_adapt.rst', 'fix_adapt_fep.rst'):
        path = os.path.join(docdir, rst)
        rows = parse_doc_tables(path)
        print('\n===== %s (%d table rows) =====' % (rst, len(rows)))

        def report(level, loc, msg):
            counts[level] += 1
            print('%s:%d: %s %s' % (os.path.relpath(loc[0]), loc[1], level, msg))

        documented = check_doc_file(path, rows, styles, parent, extracts,
                                    class_header, report)
        if not args.no_missing:
            report_missing(path, documented, styles, parent, extracts,
                           class_header, skip, report)

    print('\nsummary: %d errors, %d warnings, %d missing entries'
          % (counts['error:'], counts['warning:'], counts['note:']))
    return 1 if counts['error:'] or counts['note:'] else 0


if __name__ == '__main__':
    sys.exit(main())
