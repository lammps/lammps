// unittest/commands/test_fix_bond_react_topology.cpp -*- c++ -*-
// ----------------------------------------------------------------------
// LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
// https://www.lammps.org/, Sandia National Laboratories
// ----------------------------------------------------------------------
//
// Unit test for the REACTER (fix bond/react) TopologyMatcher refactor.
//
// Rather than hand-authoring a data file that duplicates a reaction's
// pre-reaction topology, this test inserts the *actual* pre-reaction
// molecule template into the simulation via `create_atoms ... mol`, lets
// one reaction fire, and compares the resulting bond/angle/dihedral/
// improper topology (including types) and atom types directly against
// the *actual* post-reaction molecule template -- read through LAMMPS's
// own Molecule class, not retyped by hand. It also reads the initiator
// atoms' separation directly out of the pre-reaction template's own
// Coords (via the map file's InitiatorIDs), rather than requiring each
// test case to hand-specify Rmin/Rmax values matching whatever geometry
// happens to be baked into that particular template.
//
// This makes the test data-driven and self-discovering: TEST_P is
// instantiated over every {stem}.pre/{stem}.post/{stem}.map triple found
// directly inside a subdirectory of bond-react-topology/ (see
// discover_reaction_library() below), plus every such triple found one
// level further down inside a top-level NEGATIVES/ folder, if one
// exists. A folder can hold more than one reaction this way -- e.g.
// `tiny_nylon/rxn1_stp1.{pre,post,map}` and a later
// `tiny_nylon/rxn2_stp1.{pre,post,map}` both become their own `TEST_P`
// case. Adding a new reaction to the library is just dropping a new
// {pre,post,map} triple somewhere under bond-react-topology/ (under
// NEGATIVES/ if it's a counterexample that must legitimately produce
// zero reactions) -- no C++ changes, and no separate list to keep in
// sync. Each reaction's bond/angle/dihedral/improper type counts, and
// whether it uses plain numeric types or type labels, are read directly
// out of its own template files rather than hand-specified per entry.
//
// IMPORTANT -- symmetric/indistinguishable atoms: this test does NOT
// assume simulation atom tag N corresponds to post-reaction atom-role N
// in the resulting topology. It only assumes that during *insertion*,
// tag N corresponds to *pre*-reaction template atom N (guaranteed by
// construction -- see below). Once TopologyMatcher has run, comparison
// against the post-reaction template is done via topology_isomorphic():
// a bounded search for *any* type- and connectivity-consistent bijection
// between the live simulation's atoms and the post-reaction template's
// atoms. This matters because when a reaction template contains
// genuinely graph-symmetric atoms (e.g. the two H's on a spectator -CH2-
// group, or the three on a -CH3), TopologyMatcher's search has no reason
// to prefer one of the symmetric candidates over another, and different
// (but equally valid) choices are possible depending on internal search
// order; a comparison that hard-codes one specific labeling would fail
// on an equally-correct alternate outcome. See the bundled
// "SymmetricSpectatorGroup" example for a worked case where the two
// pre-reaction-symmetric atoms are deliberately given *different*
// post-reaction bond types, so this distinction is actually exercised
// (see bond-react-topology/README.md for the full walkthrough).
//
// Scope / assumptions:
//  - runs on any number of MPI ranks. The inserted template is centered
//    on the origin, so with more than one rank it straddles subdomain
//    boundaries and the reaction has to be assembled from owned + ghost
//    atoms on several ranks. The live topology is gathered from all
//    ranks (see graph_from_live_atoms()) before it is compared, and the
//    comparison runs identically on every rank. The ghost cutoff is
//    widened to span a whole template (see build_system()), since fix
//    bond/react needs every template atom of a candidate site visible
//    on the rank that owns the initiator atom.
//  - the reaction creates and deletes no atoms, i.e. every atom in the
//    post-reaction template also appears in the pre-reaction template
//    per the map file's Equivalences section (the map file must still
//    be authored correctly for the *simulation* to behave correctly,
//    even though this test's own verification no longer reads that
//    section -- see discover_reaction_library() below).
// Reactions using CreateIDs/DeleteIDs are not yet supported; see
// discover_reaction_library() for how to extend this.

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "molecule.h"
#include "variable.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <mpi.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifndef BOND_REACT_TEST_DIR
// Fallback for standalone compilation; the real build should define this
// via target_compile_definitions(), pointing at
// unittest/commands/bond-react-topology. Unlike most test_*.cpp fixtures
// in this directory (which construct their inputs inline via
// input->one() and need no auxiliary files), this one is built around a
// growable *library* of molecule template/map files -- see "Extending
// the library" in bond-react-topology/README.md -- so its data files
// get their own subfolder to avoid cluttering unittest/commands as that
// library grows, rather than living beside the .cpp like a one-off
// fixture's inputs would.
#define BOND_REACT_TEST_DIR "unittest/commands/bond-react-topology"
#endif

using LAMMPS_NS::LAMMPS;
using LAMMPS_NS::Molecule;
using LAMMPS_NS::tagint;

// ---------------------------------------------------------------------
// One reaction, discovered from a {stem}.pre/{stem}.post/{stem}.map
// triple somewhere under bond-react-topology/ -- see
// discover_reaction_library() below for how every field here is filled
// in automatically from those three files: nobody hand-maintains type
// counts or type labels.
// ---------------------------------------------------------------------
struct ReactionTestCase {
    std::string name;                                   // folder name, plus stem if >1 per folder
    std::string pre_file, post_file, map_file;           // relative to BOND_REACT_TEST_DIR
    int natomtypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes;
    std::vector<std::string> atom_type_labels;           // empty => plain numeric types
    std::vector<std::string> bond_type_labels;
    std::vector<std::string> angle_type_labels;
    std::vector<std::string> dihedral_type_labels;
    std::vector<std::string> improper_type_labels;
    // create_box's extra/{bond,angle,dihedral,improper,special}/per/atom
    // headroom, sized per example (see max_atom_degree()) rather than a
    // single fixed guess -- a small hand-picked template can get away
    // with a generic constant, but a denser real-world one (many bonds
    // per atom, branching, rings) needs more, and a fixed constant
    // either wastes it on small examples or is too tight for large ones.
    int extra_bond_per_atom, extra_angle_per_atom, extra_dihedral_per_atom;
    int extra_improper_per_atom, extra_special_per_atom;
    // true (the default) for an ordinary example: within the initiator
    // atoms' own cutoff, the reaction must fire and match the post-
    // template. false for a deliberate counterexample -- an insertion
    // that's well within range but must still legitimately produce zero
    // reactions (e.g. a constraint in the map file that's never
    // satisfied, or atom types that don't qualify as bonding partners)
    // -- signaled by the reaction living under a top-level NEGATIVES/
    // folder; see discover_reaction_library().
    bool expect_reaction = true;
    // false (the default): the system this reaction fires against is
    // just a verbatim insertion of the pre-reaction template itself
    // (build_system() creates the box and does `create_atoms ... mol
    // mol_pre ...`). true if a `{stem}.system` file was found alongside
    // the usual pre/post/map triple -- see "Custom systems for
    // counterexamples" in bond-react-topology/README.md for why this
    // exists: some real counterexamples (e.g. a match that's ambiguous
    // only because of a *pre-existing* topology the template alone
    // doesn't capture, like two candidate insertion points overlapping
    // on a shared backbone) can't be represented by inserting the
    // template as-is, since the template trivially matches itself. When
    // true, build_system() sources system_script_file verbatim instead
    // of auto-building the box, and MatchesExpectedOutcome compares the
    // live system's own before/after topology instead of the abstract
    // template's -- so custom_system reactions are only meaningful as
    // counterexamples (expect_reaction == false; enforced at runtime).
    bool custom_system = false;
    std::string system_script_file;                     // relative to BOND_REACT_TEST_DIR
    double custom_rmin = 0.0, custom_rmax = 0.0;         // from the paired `{stem}.range` file
};

// Lets GoogleTest print a readable test parameter (e.g. in --gtest_list_tests
// and failure messages) instead of a raw 344-byte hex dump of the struct.
void PrintTo(const ReactionTestCase &tc, std::ostream *os)
{
    *os << tc.name << " (" << tc.pre_file << ")";
}

namespace {

using BondKey = std::tuple<int, int, int>;              // type, atomlo, atomhi
using AngleKey = std::tuple<int, int, int, int>;         // type, central, outerlo, outerhi
using DihedralKey = std::tuple<int, int, int, int, int>; // type, a, b, c, d (canonical direction)
using ImproperKey = std::tuple<int, int, int, int, int>; // type, a, b, c, d (as-stored order)

std::string test_file(const std::string &name)
{
    return std::string(BOND_REACT_TEST_DIR) + "/" + name;
}

std::string trim(const std::string &s)
{
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::vector<std::string> read_lines(const std::string &path)
{
    std::ifstream in(path);
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);) lines.push_back(line);
    return lines;
}

// ---- auto-discovery: turning a bond-react-topology/<Example>/ folder
// (with its *.pre/*.post/*.map files) into a full ReactionTestCase, with
// no hand-authored type counts or labelmap lists ------------------------

struct TypeInfo {
    int count = 0;
    std::vector<std::string> labels; // empty => plain numeric types 1..count
};

// Reads the second whitespace-separated token off every entry line of
// the named section (e.g. "Bonds") -- that's always the type column, for
// Types/Bonds/Angles/Dihedrals/Impropers alike -- across both template
// files, and merges what's found. A section that's absent from a file
// (LAMMPS mol files omit e.g. "Impropers" entirely when that count is 0)
// contributes nothing, which is exactly what's wanted.
//
// If every token found anywhere parses as a plain nonnegative integer,
// this is a numeric-typed category: `count` is the highest value seen
// and `labels` stays empty. Otherwise every token -- even ones that
// happen to look numeric, since LAMMPS forbids mixing labels and numbers
// within one category -- is treated as a type label and collected into
// `labels` in first-seen order (pre-template first, then post). That
// order is arbitrary but consistent, which is all maybe_issue_labelmap()
// needs: it assigns labels[i] to numeric type i+1 and issues the
// `labelmap` command that makes the template files' own label strings
// resolve correctly, so the order need not match anything in the
// source files themselves.
TypeInfo scan_type_column(const std::string &pre_path, const std::string &post_path,
                           const std::string &section)
{
    std::vector<std::string> tokens;
    for (const std::string &path : {pre_path, post_path}) {
        auto lines = read_lines(path);
        for (size_t i = 0; i < lines.size(); ++i) {
            if (trim(lines[i]) != section) continue;
            size_t j = i + 1;
            while (j < lines.size() && trim(lines[j]).empty()) ++j; // header -> blank -> entries
            for (; j < lines.size(); ++j) {
                std::string t = trim(lines[j]);
                if (t.empty()) break; // blank line ends the section
                std::istringstream iss(t);
                std::string id_tok, type_tok;
                if (iss >> id_tok >> type_tok) tokens.push_back(type_tok);
            }
            break; // at most one instance of this section keyword per file
        }
    }

    TypeInfo info;
    if (tokens.empty()) return info; // category unused by this example: count=0

    bool all_numeric = true;
    int max_val = 0;
    for (auto &tok : tokens) {
        if (tok.empty() || !std::all_of(tok.begin(), tok.end(), [](unsigned char c) {
                return std::isdigit(c);
            })) {
            all_numeric = false;
            break;
        }
        max_val = std::max(max_val, std::stoi(tok));
    }
    if (all_numeric) {
        info.count = max_val;
        return info;
    }

    std::set<std::string> seen;
    for (auto &tok : tokens)
        if (seen.insert(tok).second) info.labels.push_back(tok);
    info.count = (int) info.labels.size();
    return info;
}

// The highest number of times any single atom ID appears across the
// atom-index columns (every column after id+type) of the named section
// in one file -- e.g. for "Bonds", how many bonds any one atom is party
// to. This is a safe upper bound on the create_box extra/*/per/atom
// headroom that atom needs, regardless of which specific atom LAMMPS
// internally designates as a given bond/angle/dihedral/improper's
// "owner" for storage purposes: no atom can need more storage than
// interactions it actually appears in, so reserving against total
// appearances is always enough even if not perfectly tight.
int max_atom_degree(const std::string &path, const std::string &section)
{
    std::map<int, int> degree;
    auto lines = read_lines(path);
    for (size_t i = 0; i < lines.size(); ++i) {
        if (trim(lines[i]) != section) continue;
        size_t j = i + 1;
        while (j < lines.size() && trim(lines[j]).empty()) ++j;
        for (; j < lines.size(); ++j) {
            std::string t = trim(lines[j]);
            if (t.empty()) break;
            std::istringstream iss(t);
            std::string id_tok, type_tok, atom_tok;
            if (!(iss >> id_tok >> type_tok)) continue;
            while (iss >> atom_tok) {
                try {
                    ++degree[std::stoi(atom_tok)];
                } catch (const std::exception &) {
                    // not a plain integer atom id -- shouldn't happen for
                    // InitiatorIDs/Coords-style columns, only ever seen
                    // here for the id/type columns already consumed above
                }
            }
        }
        break;
    }
    int mx = 0;
    for (auto &[id, d] : degree) mx = std::max(mx, d);
    return mx;
}

// Total atom count from a template's own header line ("N atoms").
int read_natoms_header(const std::string &path)
{
    for (auto &l : read_lines(path)) {
        std::istringstream iss(l);
        int n;
        std::string kw;
        if (iss >> n >> kw && kw == "atoms") return n;
    }
    return 0;
}

// Reads "rmin rmax" (whitespace-separated) from a custom-system
// reaction's `{stem}.range` file -- see ReactionTestCase::custom_system.
// initiator_gap_distance() can't supply these for a custom system, since
// it derives Rmin/Rmax from the *template's* own Coords, which have no
// relationship to the real system's geometry in this case.
bool parse_range_file(const std::string &path, double &rmin, double &rmax)
{
    std::ifstream in(path);
    return static_cast<bool>(in >> rmin >> rmax);
}

// GoogleTest parameterized-test names must be valid identifier
// fragments; a folder name may contain characters ('-', etc.) that
// aren't, so anything not alnum/underscore is replaced.
std::string sanitize_test_name(std::string s)
{
    for (auto &c : s)
        if (!std::isalnum((unsigned char) c) && c != '_') c = '_';
    if (s.empty() || std::isdigit((unsigned char) s[0])) s = "Example_" + s;
    return s;
}

// Returns every file directly inside `dir` with extension `ext` (".pre",
// ".post", or ".map"), keyed by stem (filename without extension) --
// e.g. "rxn1_stp1.pre" contributes {"rxn1_stp1": ".../rxn1_stp1.pre"}.
// A folder can hold more than one reaction this way: each stem that has
// all three extensions present is one reaction (see
// load_reactions_in_dir() below); a stem missing one of the three is
// silently not a reaction (lets a folder carry incidental extra files,
// like `in.rxntest`, without confusing the scan).
std::map<std::string, std::string> collect_by_extension(const std::filesystem::path &dir,
                                                          const std::string &ext)
{
    std::map<std::string, std::string> found;
    for (auto &entry : std::filesystem::directory_iterator(dir)) {
        if (!entry.is_regular_file() || entry.path().extension() != ext) continue;
        found[entry.path().stem().string()] = entry.path().string();
    }
    return found;
}

// Builds one ReactionTestCase per {stem}.pre/{stem}.post/{stem}.map
// triple found directly inside `dir` (a folder can hold several
// reactions this way -- e.g. `tiny_nylon/rxn1_stp1.*` and
// `tiny_nylon/rxn1_stp2.*`). All three files of a triple must share the
// exact same stem -- rename them if a source example doesn't already
// (see "Layout" in bond-react-topology/README.md). `root` is
// BOND_REACT_TEST_DIR --
// pre/post/map paths are stored relative to it regardless of how deep
// under it `dir` sits (see NEGATIVES handling in
// discover_reaction_library()). `name_prefix` is prepended to every
// test name, so examples with the same folder name under different
// parents (e.g. NEGATIVES/) don't collide. `expect_reaction` is fixed
// for every reaction found in `dir` -- true for an ordinary example,
// false for everything under NEGATIVES/.
std::vector<ReactionTestCase> load_reactions_in_dir(const std::filesystem::path &dir,
                                                     const std::filesystem::path &root,
                                                     const std::string &name_prefix,
                                                     bool expect_reaction)
{
    std::vector<ReactionTestCase> found;

    std::map<std::string, std::string> pres = collect_by_extension(dir, ".pre");
    std::map<std::string, std::string> posts = collect_by_extension(dir, ".post");
    std::map<std::string, std::string> maps = collect_by_extension(dir, ".map");
    // Optional, only for a reaction that needs a custom system -- see
    // ReactionTestCase::custom_system. Both files must be present
    // together for a stem to be treated as custom-system; if only one
    // shows up, it's silently ignored and that stem falls back to the
    // ordinary verbatim-template-insertion path (most likely surfacing
    // as a confusing extra/unused file rather than a crash, which is
    // preferable to guessing at a half-specified custom system).
    std::map<std::string, std::string> systems = collect_by_extension(dir, ".system");
    std::map<std::string, std::string> ranges = collect_by_extension(dir, ".range");

    std::vector<std::string> stems;
    for (auto &[stem, path] : pres)
        if (posts.count(stem) && maps.count(stem)) stems.push_back(stem);
    std::sort(stems.begin(), stems.end()); // deterministic test order

    std::string dir_name = sanitize_test_name(dir.filename().string());

    for (auto &stem : stems) {
        const std::string &pre = pres[stem];
        const std::string &post = posts[stem];
        const std::string &map = maps[stem];

        ReactionTestCase tc;
        // Always folder name + file stem, unconditionally -- regardless
        // of how many reactions this particular folder happens to hold
        // right now. This keeps every test's name stable as siblings are
        // added or removed (no renaming when a second reaction shows up
        // in a folder that used to have only one), and it means the file
        // stem itself is free to be anything the author wants -- it
        // never needs to relate to the folder name, since the folder
        // name is always part of the test name regardless.
        tc.name = name_prefix + dir_name + "_" + sanitize_test_name(stem);
        tc.pre_file = std::filesystem::relative(pre, root).string();
        tc.post_file = std::filesystem::relative(post, root).string();
        tc.map_file = std::filesystem::relative(map, root).string();

        auto fill = [&](const char *section, int &count, std::vector<std::string> &labels) {
            TypeInfo info = scan_type_column(pre, post, section);
            count = info.count;
            labels = info.labels;
        };
        fill("Types", tc.natomtypes, tc.atom_type_labels);
        fill("Bonds", tc.nbondtypes, tc.bond_type_labels);
        fill("Angles", tc.nangletypes, tc.angle_type_labels);
        fill("Dihedrals", tc.ndihedraltypes, tc.dihedral_type_labels);
        fill("Impropers", tc.nimpropertypes, tc.improper_type_labels);

        // +2 headroom on top of the tightest-sufficient value: cheap
        // insurance against an off-by-one in exactly how LAMMPS internally
        // assigns bond/angle/dihedral/improper "ownership" for storage,
        // since max_atom_degree() only bounds it from the reservation
        // side, not from having replicated that internal logic exactly.
        auto degree_extra = [&](const char *section) {
            return std::max(max_atom_degree(pre, section), max_atom_degree(post, section)) + 2;
        };
        tc.extra_bond_per_atom = degree_extra("Bonds");
        tc.extra_angle_per_atom = degree_extra("Angles");
        tc.extra_dihedral_per_atom = degree_extra("Dihedrals");
        tc.extra_improper_per_atom = degree_extra("Impropers");
        // Special-neighbor (1-2/1-3/1-4) slots aren't bounded by a simple
        // per-section degree the way bonds/angles/etc. are -- they grow
        // with local connectivity density (branching, rings). A
        // template this size can't have more special neighbors than
        // atoms - 1 though, so that's used as a simple, always-safe cap
        // rather than trying to derive a tighter bound from the bond
        // graph's shape.
        tc.extra_special_per_atom =
            std::max(1, std::max(read_natoms_header(pre), read_natoms_header(post)) - 1);

        tc.expect_reaction = expect_reaction;

        if (systems.count(stem) && ranges.count(stem)) {
            tc.custom_system = true;
            tc.system_script_file = std::filesystem::relative(systems[stem], root).string();
            // A malformed `.range` file (not two whitespace-separated
            // numbers) leaves Rmin/Rmax at 0.0/0.0 -- fix bond/react will
            // then simply never see anything in range, which surfaces as
            // an easy-to-diagnose "expected 0 reactions, got 0, but for
            // the wrong reason" rather than a crash here.
            parse_range_file(ranges[stem], tc.custom_rmin, tc.custom_rmax);
        }

        found.push_back(std::move(tc));
    }
    return found;
}

// Scans BOND_REACT_TEST_DIR for reaction examples, two ways:
//  - every {stem}.pre/{stem}.post/{stem}.map triple found directly in
//    any immediate subdirectory (other than NEGATIVES itself) is added
//    as an ordinary example, expecting one reaction. A folder may
//    contribute more than one reaction this way -- e.g. `tiny_nylon/`
//    can hold both `rxn1_stp1.{pre,post,map}` and
//    `rxn2_stp1.{pre,post,map}`, each becoming its own `TEST_P` case.
//  - every such triple found one level inside a top-level NEGATIVES/
//    folder, if one exists, is added as a counterexample
//    (expect_reaction = false) instead -- NEGATIVES/ is where every
//    example that must legitimately produce zero reactions lives, kept
//    together and apart from the ordinary examples above.
// Either way, if the same stem also has a `.system` file (and a paired
// `.range` file), that reaction gets ReactionTestCase::custom_system =
// true -- see that field's comment, and "Custom systems for
// counterexamples" in bond-react-topology/README.md.
// Everything else about a ReactionTestCase -- type counts, and whether
// each category uses plain numbers or type labels -- is derived straight
// from the two template files via scan_type_column(), so growing the
// library (see bond-react-topology/README.md) is purely a matter of
// adding new files; this file never needs to change.
std::vector<ReactionTestCase> discover_reaction_library()
{
    std::vector<ReactionTestCase> lib;
    std::filesystem::path root(BOND_REACT_TEST_DIR);
    if (!std::filesystem::is_directory(root)) return lib;

    std::vector<std::filesystem::path> dirs;
    std::filesystem::path negatives_dir;
    for (auto &entry : std::filesystem::directory_iterator(root)) {
        if (!entry.is_directory()) continue;
        if (entry.path().filename() == "NEGATIVES") {
            negatives_dir = entry.path();
            continue; // handled separately below, not as an example itself
        }
        dirs.push_back(entry.path());
    }
    std::sort(dirs.begin(), dirs.end()); // deterministic test order

    for (auto &dir : dirs) {
        auto reactions = load_reactions_in_dir(dir, root, /*name_prefix=*/"",
                                                /*expect_reaction=*/true);
        for (auto &tc : reactions) lib.push_back(std::move(tc));
    }

    if (!negatives_dir.empty() && std::filesystem::is_directory(negatives_dir)) {
        std::vector<std::filesystem::path> neg_dirs;
        for (auto &entry : std::filesystem::directory_iterator(negatives_dir))
            if (entry.is_directory()) neg_dirs.push_back(entry.path());
        std::sort(neg_dirs.begin(), neg_dirs.end());

        for (auto &dir : neg_dirs) {
            auto reactions = load_reactions_in_dir(dir, root, /*name_prefix=*/"NEGATIVES_",
                                                    /*expect_reaction=*/false);
            for (auto &tc : reactions) lib.push_back(std::move(tc));
        }
    }

    return lib;
}

// Canonicalize a-b-c-d against its exact reverse d-c-b-a so a dihedral
// is only counted once regardless of which end it was listed from. This
// is a fixed, mapping-independent canonicalization of *direction*
// (always valid for any dihedral, symmetric atoms or not) -- separate
// from, and composable with, the atom-relabeling search below.
DihedralKey canonical_dihedral(int type, int a, int b, int c, int d)
{
    std::tuple<int, int, int, int> fwd(a, b, c, d), rev(d, c, b, a);
    auto &chosen = (fwd <= rev) ? fwd : rev;
    return DihedralKey(type, std::get<0>(chosen), std::get<1>(chosen), std::get<2>(chosen),
                        std::get<3>(chosen));
}

// ---- a topology, in some atom-id space of its own -------------------

struct TopologyGraph {
    std::map<int, int> atom_type; // atom-id -> type, in this graph's own id space
    std::set<BondKey> bonds;
    std::set<AngleKey> angles;
    std::set<DihedralKey> dihedrals;
    std::set<ImproperKey> impropers;
};

// Collects the live system's topology from *all* MPI ranks: each rank
// flattens what it stores for its own atoms into one int buffer, the
// buffers are allgathered, and every rank rebuilds the identical global
// graph from them. With newton_bond on each bond/angle/dihedral/improper
// is stored on exactly one owning atom; with newton_bond off it is stored
// on several, and the std::set containers below drop those duplicates --
// so this is correct either way.
TopologyGraph graph_from_live_atoms(LAMMPS *lmp)
{
    enum { ATOM = 0, BOND, ANGLE, DIHEDRAL, IMPROPER };
    auto *atom = lmp->atom;
    std::vector<int> mine;
    auto put = [&](std::initializer_list<long long> v) {
        for (auto x : v) mine.push_back((int) x);
    };
    for (int i = 0; i < atom->nlocal; ++i) {
        tagint itag = atom->tag[i];
        put({ATOM, itag, atom->type[i]});
        for (int k = 0; k < atom->num_bond[i]; ++k)
            put({BOND, atom->bond_type[i][k], itag, atom->bond_atom[i][k]});
        for (int k = 0; k < atom->num_angle[i]; ++k)
            put({ANGLE, atom->angle_type[i][k], atom->angle_atom1[i][k], atom->angle_atom2[i][k],
                 atom->angle_atom3[i][k]});
        for (int k = 0; k < atom->num_dihedral[i]; ++k)
            put({DIHEDRAL, atom->dihedral_type[i][k], atom->dihedral_atom1[i][k],
                 atom->dihedral_atom2[i][k], atom->dihedral_atom3[i][k],
                 atom->dihedral_atom4[i][k]});
        for (int k = 0; k < atom->num_improper[i]; ++k)
            put({IMPROPER, atom->improper_type[i][k], atom->improper_atom1[i][k],
                 atom->improper_atom2[i][k], atom->improper_atom3[i][k],
                 atom->improper_atom4[i][k]});
    }

    MPI_Comm world = lmp->world;
    int nprocs = 1;
    MPI_Comm_size(world, &nprocs);
    int nmine = (int) mine.size();
    std::vector<int> counts(nprocs), displs(nprocs, 0);
    MPI_Allgather(&nmine, 1, MPI_INT, counts.data(), 1, MPI_INT, world);
    for (int p = 1; p < nprocs; ++p) displs[p] = displs[p - 1] + counts[p - 1];
    std::vector<int> all(displs[nprocs - 1] + counts[nprocs - 1]);
    MPI_Allgatherv(mine.data(), nmine, MPI_INT, all.data(), counts.data(), displs.data(),
                   MPI_INT, world);

    TopologyGraph g;
    for (size_t n = 0; n < all.size();) {
        switch (all[n]) {
            case ATOM:
                g.atom_type[all[n + 1]] = all[n + 2];
                n += 3;
                break;
            case BOND:
                g.bonds.emplace(all[n + 1], std::min(all[n + 2], all[n + 3]),
                                std::max(all[n + 2], all[n + 3]));
                n += 4;
                break;
            case ANGLE:
                g.angles.emplace(all[n + 1], all[n + 3], std::min(all[n + 2], all[n + 4]),
                                 std::max(all[n + 2], all[n + 4]));
                n += 5;
                break;
            case DIHEDRAL:
                g.dihedrals.insert(canonical_dihedral(all[n + 1], all[n + 2], all[n + 3],
                                                      all[n + 4], all[n + 5]));
                n += 6;
                break;
            case IMPROPER:
                g.impropers.emplace(all[n + 1], all[n + 2], all[n + 3], all[n + 4], all[n + 5]);
                n += 6;
                break;
            default:
                ADD_FAILURE() << "corrupt topology buffer in graph_from_live_atoms()";
                return g;
        }
    }
    return g;
}

// Builds a TopologyGraph directly from a Molecule template, using its
// own 1-indexed template-local atom IDs as this graph's id space. Types
// come from tmpl->type[]/bond_type[]/... which LAMMPS has already
// resolved to plain numeric types regardless of whether the template
// file used labels or numbers, so labels need no special handling here.
TopologyGraph graph_from_molecule(Molecule *tmpl)
{
    TopologyGraph g;
    for (int i = 0; i < tmpl->natoms; ++i) {
        int id = i + 1; // template files are 1-indexed
        g.atom_type[id] = tmpl->type[i];
        for (int k = 0; k < tmpl->num_bond[i]; ++k) {
            int other = (int) tmpl->bond_atom[i][k];
            g.bonds.emplace(tmpl->bond_type[i][k], std::min(id, other), std::max(id, other));
        }
        for (int k = 0; k < tmpl->num_angle[i]; ++k) {
            int a1 = (int) tmpl->angle_atom1[i][k];
            int a2 = (int) tmpl->angle_atom2[i][k];
            int a3 = (int) tmpl->angle_atom3[i][k];
            g.angles.emplace(tmpl->angle_type[i][k], a2, std::min(a1, a3), std::max(a1, a3));
        }
        for (int k = 0; k < tmpl->num_dihedral[i]; ++k) {
            g.dihedrals.insert(canonical_dihedral(
                tmpl->dihedral_type[i][k], (int) tmpl->dihedral_atom1[i][k],
                (int) tmpl->dihedral_atom2[i][k], (int) tmpl->dihedral_atom3[i][k],
                (int) tmpl->dihedral_atom4[i][k]));
        }
        for (int k = 0; k < tmpl->num_improper[i]; ++k) {
            g.impropers.emplace(tmpl->improper_type[i][k], (int) tmpl->improper_atom1[i][k],
                                (int) tmpl->improper_atom2[i][k], (int) tmpl->improper_atom3[i][k],
                                (int) tmpl->improper_atom4[i][k]);
        }
    }
    return g;
}

bool verify_relations(const TopologyGraph &a, const TopologyGraph &b,
                       const std::map<int, int> &assign)
{
    for (auto &[type, lo, hi] : a.bonds) {
        int mlo = std::min(assign.at(lo), assign.at(hi));
        int mhi = std::max(assign.at(lo), assign.at(hi));
        if (!b.bonds.count({type, mlo, mhi})) return false;
    }
    for (auto &[type, central, lo, hi] : a.angles) {
        int mc = assign.at(central);
        int mlo = std::min(assign.at(lo), assign.at(hi));
        int mhi = std::max(assign.at(lo), assign.at(hi));
        if (!b.angles.count({type, mc, mlo, mhi})) return false;
    }
    for (auto &[type, x1, x2, x3, x4] : a.dihedrals) {
        // a's dihedrals were canonicalized in a's own id space; the
        // canonical direction can change once mapped into b's id space
        // (it depends on the numeric id values), so re-canonicalize.
        auto mapped =
            canonical_dihedral(type, assign.at(x1), assign.at(x2), assign.at(x3), assign.at(x4));
        if (!b.dihedrals.count(mapped)) return false;
    }
    for (auto &[type, x1, x2, x3, x4] : a.impropers) {
        ImproperKey mapped(type, assign.at(x1), assign.at(x2), assign.at(x3), assign.at(x4));
        if (!b.impropers.count(mapped)) return false;
    }
    return true;
}

// Iterative Weisfeiler-Leman-style color refinement, run jointly over
// both graphs' atoms (via their bond adjacency only) so that a shared
// refined color in the two returned maps genuinely means "the same
// local bonding neighborhood, to whatever depth distinguishes it" --
// not just "the same bare force-field atom type". Two atoms that end
// up in different colors can *never* be part of any valid
// type/topology-preserving bijection between a and b (each round's
// color is a function of the previous round's colors and the bond
// types connecting them, so any true isomorphism maps equal colors to
// equal colors, by induction on the round). This makes it a strictly
// sound way to shrink topology_isomorphic()'s per-atom candidate
// lists: it can only ever split a same-type candidate set into smaller
// pieces, never wrongly exclude a genuine candidate. For a real
// molecule this typically collapses a large same-type group (e.g. many
// tens of chemically identical aromatic CH atoms, which raw atom type
// alone can't tell apart) down to just the atoms that are also
// genuinely symmetric by ring position/substituent -- usually a
// handful at most -- which is what actually keeps the backtracking
// search below tractable on real reaction templates rather than only
// small hand-designed ones.
std::pair<std::map<int, int>, std::map<int, int>> refine_colors(const TopologyGraph &a,
                                                                  const TopologyGraph &b)
{
    auto adjacency_of = [](const TopologyGraph &g) {
        std::map<int, std::vector<std::pair<int, int>>> adj; // atom -> [(bond_type, neighbor)]
        for (auto &[type, lo, hi] : g.bonds) {
            adj[lo].emplace_back(type, hi);
            adj[hi].emplace_back(type, lo);
        }
        return adj;
    };
    auto adj_a = adjacency_of(a);
    auto adj_b = adjacency_of(b);

    std::map<int, int> color_a, color_b;
    for (auto &[id, t] : a.atom_type) color_a[id] = t;
    for (auto &[id, t] : b.atom_type) color_b[id] = t;

    // A WL-style partition can only refine (never coarsen) at each
    // round, and there are at most N atoms total, so it's guaranteed to
    // stabilize within N rounds; stopping as soon as a round changes
    // nothing keeps this fast on real (mostly-asymmetric) molecules,
    // where it typically stabilizes in just 2-3 rounds.
    size_t max_rounds = a.atom_type.size() + b.atom_type.size() + 1;
    using Signature = std::pair<int, std::vector<std::pair<int, int>>>;

    for (size_t round = 0; round < max_rounds; ++round) {
        auto signature_of = [](int id, const std::map<int, int> &color,
                                const std::map<int, std::vector<std::pair<int, int>>> &adj) {
            std::vector<std::pair<int, int>> neighbor_colors;
            auto it = adj.find(id);
            if (it != adj.end())
                for (auto &[bond_type, neighbor] : it->second)
                    neighbor_colors.emplace_back(bond_type, color.at(neighbor));
            std::sort(neighbor_colors.begin(), neighbor_colors.end());
            return Signature(color.at(id), std::move(neighbor_colors));
        };

        // One shared signature->id table for both graphs this round, so
        // an atom in `a` and an atom in `b` with identical signatures
        // are assigned the exact same new color.
        std::map<Signature, int> signature_to_color;
        int next_color = 0;
        std::map<int, int> new_color_a, new_color_b;
        for (auto &[id, unused] : a.atom_type) {
            Signature sig = signature_of(id, color_a, adj_a);
            auto it = signature_to_color.find(sig);
            new_color_a[id] =
                (it != signature_to_color.end()) ? it->second : (signature_to_color[sig] = next_color++);
        }
        for (auto &[id, unused] : b.atom_type) {
            Signature sig = signature_of(id, color_b, adj_b);
            auto it = signature_to_color.find(sig);
            new_color_b[id] =
                (it != signature_to_color.end()) ? it->second : (signature_to_color[sig] = next_color++);
        }

        bool changed = (new_color_a != color_a) || (new_color_b != color_b);
        color_a = std::move(new_color_a);
        color_b = std::move(new_color_b);
        if (!changed) break;
    }

    return {color_a, color_b};
}

// Searches for *some* bijection between a's atom-ids and b's atom-ids
// under which every one of a's bonds/angles/dihedrals/impropers maps
// onto one of b's (types included) -- i.e. whether a and b describe the
// same topology up to a consistent relabeling. Candidates for each atom
// are narrowed by refine_colors() before any backtracking starts (see
// its own comment for why that's sound), and partial bond-consistency
// is checked incrementally during the search itself, both purely for
// pruning; for a fully asymmetric template (or a symmetric group of a
// handful of atoms) this is effectively instant, and for a real
// molecule with many tens of same-type atoms (e.g. aromatic CH's),
// color refinement typically leaves only genuinely symmetric groups of
// a handful of atoms each for the backtracking to actually branch on.
bool topology_isomorphic(const TopologyGraph &a, const TopologyGraph &b, std::string *why)
{
    auto fail = [&](const std::string &msg) {
        if (why) *why = msg;
        return false;
    };

    if (a.atom_type.size() != b.atom_type.size())
        return fail("different atom counts (" + std::to_string(a.atom_type.size()) + " vs " +
                    std::to_string(b.atom_type.size()) + ")");
    if (a.bonds.size() != b.bonds.size())
        return fail("different bond counts (" + std::to_string(a.bonds.size()) + " vs " +
                    std::to_string(b.bonds.size()) + ")");
    if (a.angles.size() != b.angles.size())
        return fail("different angle counts (" + std::to_string(a.angles.size()) + " vs " +
                    std::to_string(b.angles.size()) + ")");
    if (a.dihedrals.size() != b.dihedrals.size())
        return fail("different dihedral counts (" + std::to_string(a.dihedrals.size()) + " vs " +
                    std::to_string(b.dihedrals.size()) + ")");
    if (a.impropers.size() != b.impropers.size())
        return fail("different improper counts (" + std::to_string(a.impropers.size()) + " vs " +
                    std::to_string(b.impropers.size()) + ")");

    std::map<int, int> a_type_counts, b_type_counts;
    for (auto &[id, t] : a.atom_type) a_type_counts[t]++;
    for (auto &[id, t] : b.atom_type) b_type_counts[t]++;
    if (a_type_counts != b_type_counts) return fail("different atom-type multiset");

    auto [color_a, color_b] = refine_colors(a, b);
    std::map<int, std::vector<int>> b_by_color;
    for (auto &[id, c] : color_b) b_by_color[c].push_back(id);

    std::map<int, std::vector<int>> candidates;
    for (auto &[id, t] : a.atom_type) {
        candidates[id] = b_by_color[color_a.at(id)];
        if (candidates[id].empty())
            return fail("atom " + std::to_string(id) +
                        " has no counterpart with the same bonding environment");
    }

    // Search order: start from the most-constrained atom, then always
    // take next the unordered atom with the most bonds back to atoms
    // already ordered (ties: fewer candidates, then lower id). Growing
    // the assignment outward along the bond graph like this means every
    // bond/angle/dihedral/improper becomes fully determined -- and so
    // checkable -- as early as possible, instead of only once some
    // unrelated far-away atom happens to be reached.
    std::map<int, std::vector<int>> adj;
    for (auto &[type, lo, hi] : a.bonds) {
        adj[lo].push_back(hi);
        adj[hi].push_back(lo);
    }
    std::vector<int> a_atoms;
    std::map<int, int> order_of; // a atom id -> position in a_atoms
    {
        std::map<int, int> links; // unordered atom -> #bonds into ordered set
        for (auto &[id, t] : a.atom_type) links[id] = 0;
        while (!links.empty()) {
            auto best = links.begin();
            for (auto it = links.begin(); it != links.end(); ++it) {
                if (it->second != best->second) {
                    if (it->second > best->second) best = it;
                    continue;
                }
                if (candidates[it->first].size() < candidates[best->first].size()) best = it;
            }
            int id = best->first;
            links.erase(best);
            order_of[id] = (int) a_atoms.size();
            a_atoms.push_back(id);
            for (int nb : adj[id]) {
                auto it = links.find(nb);
                if (it != links.end()) ++it->second;
            }
        }
    }

    // Bucket every relation by the search depth at which its *last* atom
    // gets assigned; at that depth (and only then) it is checked. Every
    // relation is thus checked exactly once, the moment it's fully
    // determined, so a wrong choice for a symmetric atom is rejected
    // right where it's made rather than after the whole remaining search
    // tree below it has been enumerated. (Checking only bonds during the
    // search and angles/dihedrals/impropers only at the very end made
    // this exponential on large templates with several symmetric groups,
    // e.g. benzoxazine's CH2/CH3 hydrogens: each bond-consistent but
    // otherwise wrong combination of those groups was fully enumerated.)
    const size_t n = a_atoms.size();
    auto depth_of = [&](std::initializer_list<int> ids) {
        int d = 0;
        for (int id : ids) d = std::max(d, order_of.at(id));
        return (size_t) d;
    };
    std::vector<std::vector<BondKey>> bonds_at(n);
    std::vector<std::vector<AngleKey>> angles_at(n);
    std::vector<std::vector<DihedralKey>> dihedrals_at(n);
    std::vector<std::vector<ImproperKey>> impropers_at(n);
    for (auto &r : a.bonds) bonds_at[depth_of({std::get<1>(r), std::get<2>(r)})].push_back(r);
    for (auto &r : a.angles)
        angles_at[depth_of({std::get<1>(r), std::get<2>(r), std::get<3>(r)})].push_back(r);
    for (auto &r : a.dihedrals)
        dihedrals_at[depth_of({std::get<1>(r), std::get<2>(r), std::get<3>(r), std::get<4>(r)})]
            .push_back(r);
    for (auto &r : a.impropers)
        impropers_at[depth_of({std::get<1>(r), std::get<2>(r), std::get<3>(r), std::get<4>(r)})]
            .push_back(r);

    std::map<int, int> assign;
    std::set<int> used_b;

    auto consistent_at = [&](size_t idx) {
        for (auto &[type, lo, hi] : bonds_at[idx]) {
            int x = assign.at(lo), y = assign.at(hi);
            if (!b.bonds.count({type, std::min(x, y), std::max(x, y)})) return false;
        }
        for (auto &[type, central, lo, hi] : angles_at[idx]) {
            int x = assign.at(lo), y = assign.at(hi);
            if (!b.angles.count({type, assign.at(central), std::min(x, y), std::max(x, y)}))
                return false;
        }
        for (auto &[type, x1, x2, x3, x4] : dihedrals_at[idx]) {
            // re-canonicalize: the canonical direction depends on the
            // numeric id values, which change under the mapping
            if (!b.dihedrals.count(canonical_dihedral(type, assign.at(x1), assign.at(x2),
                                                      assign.at(x3), assign.at(x4))))
                return false;
        }
        for (auto &[type, x1, x2, x3, x4] : impropers_at[idx]) {
            if (!b.impropers.count(
                    {type, assign.at(x1), assign.at(x2), assign.at(x3), assign.at(x4)}))
                return false;
        }
        return true;
    };

    // Hard cap on search effort, so that a pathological case turns into
    // a clear test failure instead of an apparently hung test binary.
    const long max_steps = 20000000;
    long steps = 0;
    bool budget_exceeded = false;

    std::function<bool(size_t)> rec = [&](size_t idx) -> bool {
        if (idx == n) return verify_relations(a, b, assign); // cheap final sanity check
        int ai = a_atoms[idx];
        for (int bi : candidates[ai]) {
            if (used_b.count(bi)) continue;
            if (++steps > max_steps) {
                budget_exceeded = true;
                return false;
            }
            assign[ai] = bi;
            used_b.insert(bi);
            if (consistent_at(idx) && rec(idx + 1)) return true;
            used_b.erase(bi);
            assign.erase(ai);
            if (budget_exceeded) return false;
        }
        return false;
    };

    if (rec(0)) return true;
    if (budget_exceeded)
        return fail("atom relabeling search gave up after " + std::to_string(max_steps) +
                    " steps without finding a match (template too symmetric for this checker?)");
    if (rec(0)) return true;
    return fail("no consistent atom relabeling found matching the template "
                "(topology or types genuinely differ, not just a symmetric relabeling)");
}

// ---- raw-text helpers over this harness's own template/map files -----
//
// These only ever parse atom IDs and coordinates (InitiatorIDs, Coords),
// never type columns, so they are unaffected by whether a template uses
// numeric types or type labels.

// Reads the two atom IDs following the mandatory "InitiatorIDs" section.
std::vector<int> parse_initiator_ids(const std::string &path)
{
    auto lines = read_lines(path);
    std::vector<int> ids;
    for (size_t i = 0; i < lines.size(); ++i) {
        if (trim(lines[i]) != "InitiatorIDs") continue;
        for (size_t j = i + 1; j < lines.size() && ids.size() < 2; ++j) {
            std::string t = trim(lines[j]);
            if (t.empty()) continue;
            std::istringstream iss(t);
            int id;
            if (iss >> id) ids.push_back(id);
        }
        break;
    }
    return ids;
}

struct Vec3 {
    double x = 0, y = 0, z = 0;
};

// Reads the header's "N atoms" count, then exactly N "id x y z" lines
// following the "Coords" section keyword.
std::map<int, Vec3> parse_template_coords(const std::string &path)
{
    auto lines = read_lines(path);

    int natoms = -1;
    for (auto &l : lines) {
        std::istringstream iss(l);
        int n;
        std::string kw;
        if (iss >> n >> kw && kw == "atoms") {
            natoms = n;
            break;
        }
    }

    std::map<int, Vec3> coords;
    if (natoms < 0) return coords;

    for (size_t i = 0; i < lines.size(); ++i) {
        if (trim(lines[i]) != "Coords") continue;
        int read = 0;
        for (size_t j = i + 1; j < lines.size() && read < natoms; ++j) {
            std::string t = trim(lines[j]);
            if (t.empty()) continue;
            std::istringstream iss(t);
            int id;
            double x, y, z;
            if (iss >> id >> x >> y >> z) {
                coords[id] = {x, y, z};
                ++read;
            }
        }
        break;
    }
    return coords;
}

// Distance (template/simulation distance units) between the two
// InitiatorIDs atoms, as placed in the pre-reaction template's own
// Coords -- i.e. exactly the gap fix bond/react's Rmin/Rmax will see
// once the template is inserted unchanged via create_atoms.
double initiator_gap_distance(const ReactionTestCase &tc)
{
    std::vector<int> ids = parse_initiator_ids(test_file(tc.map_file));
    if (ids.size() != 2)
        ADD_FAILURE() << "map file " << tc.map_file << " did not yield exactly 2 InitiatorIDs";
    std::map<int, Vec3> coords = parse_template_coords(test_file(tc.pre_file));
    const Vec3 &p1 = coords.at(ids[0]);
    const Vec3 &p2 = coords.at(ids[1]);
    double dx = p1.x - p2.x, dy = p1.y - p2.y, dz = p1.z - p2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Largest distance between any two atoms in the pre-reaction template's
// own Coords. fix bond/react walks the whole template outward from the
// initiator atoms, so every atom of a candidate site must be available
// (owned or ghost) on the rank doing the matching; a ghost cutoff of at
// least this span guarantees that regardless of how the site is split
// across subdomains.
double template_span(const ReactionTestCase &tc)
{
    std::map<int, Vec3> coords = parse_template_coords(test_file(tc.pre_file));
    double maxsq = 0.0;
    for (auto &[i, p] : coords)
        for (auto &[j, q] : coords) {
            double dx = p.x - q.x, dy = p.y - q.y, dz = p.z - q.z;
            maxsq = std::max(maxsq, dx * dx + dy * dy + dz * dz);
        }
    return std::sqrt(maxsq);
}

// Issues `labelmap <kind> 1 label1 2 label2 ...` if any labels were
// given; a no-op (matching plain numeric types) if labels is empty.
void maybe_issue_labelmap(LAMMPS *lmp, const char *kind, const std::vector<std::string> &labels)
{
    if (labels.empty()) return;
    std::string cmd = std::string("labelmap ") + kind;
    for (size_t i = 0; i < labels.size(); ++i) cmd += " " + std::to_string(i + 1) + " " + labels[i];
    lmp->input->one(cmd.c_str());
}

} // namespace

// The reaction library, built once at static-initialization time by
// scanning bond-react-topology/ -- see discover_reaction_library().
const std::vector<ReactionTestCase> kReactionLibrary = discover_reaction_library();

class FixBondReactTopologyTest : public ::testing::TestWithParam<ReactionTestCase> {
protected:
    LAMMPS *lmp = nullptr;

    void SetUp() override
    {
        const char *args[] = {"lammps", "-log", "none", "-echo", "none", "-nocite"};
        int nargs = sizeof(args) / sizeof(char *);
        lmp = new LAMMPS(nargs, (char **) args, MPI_COMM_WORLD);
    }

    void TearDown() override
    {
        delete lmp;
        lmp = nullptr;
    }

    // Builds the system, inserts/reads the atoms it will react against,
    // and defines the fix bond/react. Does not run any timesteps.
    // rmin/rmax control whether the reaction is expected to fire.
    //
    // Ordinarily (tc.custom_system == false) the "system" is just the
    // pre-reaction template itself, inserted verbatim via create_atoms,
    // and rmin/rmax are derived from the initiator atoms' own separation
    // baked into that template (see initiator_gap_distance()).
    //
    // For a custom-system reaction (tc.custom_system == true --
    // see ReactionTestCase::custom_system and "Custom systems for
    // counterexamples" in bond-react-topology/README.md), everything
    // this method would otherwise auto-generate -- units, atom_style,
    // boundary, the box and its atoms (via read_data or its own
    // create_box+create_atoms), pair/bond/angle/dihedral/improper
    // styles, and special_bonds -- is instead read verbatim from
    // tc.system_script_file, and rmin/rmax are whatever the caller
    // passed in (from tc.custom_rmin/custom_rmax, read from the
    // reaction's `.range` file) rather than derived from a template.
    void build_system(const ReactionTestCase &tc, double rmin, double rmax)
    {
        lmp->input->one("atom_modify map array");

        if (tc.custom_system) {
            // The script's working directory is wherever ctest/the test
            // binary itself is running from, not the example's own
            // folder, so it can't just say e.g. "read_data foo.data" and
            // expect that to resolve -- expose the example's own folder
            // as a LAMMPS string variable it can build paths from, e.g.
            // "read_data ${exampledir}/foo.data" (see "Custom systems
            // for counterexamples" in bond-react-topology/README.md).
            std::string script_path = test_file(tc.system_script_file);
            std::string example_dir = std::filesystem::path(script_path).parent_path().string();
            lmp->input->one(("variable exampledir string \"" + example_dir + "\"").c_str());
            lmp->input->file(script_path.c_str());
        } else {
            lmp->input->one("units real");
            lmp->input->one("atom_style molecular");
            lmp->input->one("boundary p p p");
            lmp->input->one("region simbox block -50 50 -50 50 -50 50");

            // extra/*/per/atom headroom is sized per example (see
            // max_atom_degree() in discover_reaction_library()), not a
            // fixed guess -- a small hand-picked template can get away
            // with a generic constant, but a real, denser template
            // (more bonds per atom, branching, rings, like a distilled
            // tiny_nylon crosslink) needs more, and LAMMPS rejects
            // create_atoms/fix bond/react outright ("Molecule
            // topology/atom exceeds system topology/atom") if a
            // template's own per-atom bond/angle/dihedral/improper/
            // special-neighbor count exceeds what create_box reserved.
            char box_cmd[512];
            std::snprintf(box_cmd, sizeof(box_cmd),
                           "create_box %d simbox bond/types %d angle/types %d dihedral/types %d"
                           " improper/types %d extra/bond/per/atom %d extra/angle/per/atom %d"
                           " extra/dihedral/per/atom %d extra/improper/per/atom %d"
                           " extra/special/per/atom %d",
                           tc.natomtypes, tc.nbondtypes, tc.nangletypes, tc.ndihedraltypes,
                           tc.nimpropertypes, tc.extra_bond_per_atom, tc.extra_angle_per_atom,
                           tc.extra_dihedral_per_atom, tc.extra_improper_per_atom,
                           tc.extra_special_per_atom);
            lmp->input->one(box_cmd);
        }

        // molecule files cannot define their own label map (only data
        // files can, via an Atom Type Labels section); if a template
        // uses type labels, they must already be defined before it's
        // read -- see scan_type_column()/discover_reaction_library().
        // Applies equally to a custom system: the labelmap just needs
        // to come after the box exists (created above, one way or the
        // other) and before the `molecule` commands below.
        maybe_issue_labelmap(lmp, "atom", tc.atom_type_labels);
        maybe_issue_labelmap(lmp, "bond", tc.bond_type_labels);
        maybe_issue_labelmap(lmp, "angle", tc.angle_type_labels);
        maybe_issue_labelmap(lmp, "dihedral", tc.dihedral_type_labels);
        maybe_issue_labelmap(lmp, "improper", tc.improper_type_labels);

        if (!tc.custom_system) {
            lmp->input->one("mass * 12.0");

            char pair_cmd[64];
            std::snprintf(pair_cmd, sizeof(pair_cmd), "pair_style zero %g", rmax + 5.0);
            lmp->input->one(pair_cmd);
            lmp->input->one("pair_coeff * *");
            // *_coeff * (or * *) expands to an explicit 1..N type range,
            // and LAMMPS rejects that range when N is 0 ("Invalid range
            // string: *"), so each bond/angle/dihedral/improper
            // style+coeff pair is only issued when that reaction
            // library entry's template actually declares types of that
            // kind (a template's post-reaction side can introduce a
            // bonded-interaction category the pre-reaction side has
            // none of, e.g. an all-linear pre-reaction template with
            // zero improper types -- see "SymmetricSpectatorGroup").
            if (tc.nbondtypes > 0) {
                lmp->input->one("bond_style zero");
                lmp->input->one("bond_coeff * 0.0");
            }
            if (tc.nangletypes > 0) {
                lmp->input->one("angle_style zero");
                lmp->input->one("angle_coeff *");
            }
            if (tc.ndihedraltypes > 0) {
                lmp->input->one("dihedral_style zero");
                lmp->input->one("dihedral_coeff *");
            }
            if (tc.nimpropertypes > 0) {
                lmp->input->one("improper_style zero");
                lmp->input->one("improper_coeff *");
            }
            lmp->input->one("special_bonds lj/coul 0.0 0.0 0.0");
        }

        lmp->input->one(("molecule mol_pre " + test_file(tc.pre_file)).c_str());
        lmp->input->one(("molecule mol_post " + test_file(tc.post_file)).c_str());

        if (!tc.custom_system) {
            // create_atoms places the template's geometric center at the
            // given point, so inserting at the origin centers it in the
            // box -- and, with more than one MPI rank, right across the
            // subdomain boundaries (a default processor grid splits the
            // box at 0 along every direction with an even rank count),
            // so the reaction genuinely has to be assembled across ranks.
            lmp->input->one("create_atoms 0 single 0.0 0.0 0.0 mol mol_pre 12345 units box");
        }
        // For a custom system, tc.system_script_file already put every
        // atom the reaction will actually be tested against into place
        // (e.g. via read_data); mol_pre/mol_post above are just the
        // abstract stencils fix bond/react matches against, not
        // themselves inserted into the simulation.

        // Ghost atoms must reach across a whole template (see
        // template_span()); on a single rank this is harmless. A custom
        // system's real geometry needn't resemble the template's at all
        // (e.g. NEGATIVES/crosslink_negative's chain has 16-18 unit long
        // bonds), so there the live system's own extent is covered too.
        {
            double cut = template_span(tc) + rmax + 5.0;
            if (tc.custom_system) cut = std::max(cut, live_system_span() + 1.0);
            char comm_cmd[128];
            std::snprintf(comm_cmd, sizeof(comm_cmd), "comm_modify cutoff %g", cut);
            lmp->input->one(comm_cmd);
        }

        char rxn_cmd[512];
        std::snprintf(rxn_cmd, sizeof(rxn_cmd),
                       "fix myrxn all bond/react react rxntest all 1 %g %g mol_pre mol_post %s",
                       rmin, rmax, test_file(tc.map_file).c_str());
        lmp->input->one(rxn_cmd);

        lmp->input->one("fix 1 all nve");
        lmp->input->one("timestep 1.0");

        // Trim thermo down to just the step and this fix's own reaction
        // count (its compute_vector()/get_thermo_colname() are what
        // supplies the "f_myrxn[1]" column's label), so anyone running
        // this test's LAMMPS instance directly -- e.g. via -V/--gtest
        // output, or the standalone in.rxntest scripts -- can confirm
        // success from the run log at a glance instead of picking it out
        // of the full default thermo block.
        lmp->input->one("thermo 1");
        lmp->input->one("thermo_style custom step f_myrxn[1]");
    }

    // Largest minimum-image distance between any two atoms of the live
    // system, over all ranks. Only used for (small) custom systems.
    double live_system_span()
    {
        auto *atom = lmp->atom;
        std::vector<double> mine;
        for (int i = 0; i < atom->nlocal; ++i)
            for (int k = 0; k < 3; ++k) mine.push_back(atom->x[i][k]);
        int nprocs = 1;
        MPI_Comm_size(lmp->world, &nprocs);
        int nmine = (int) mine.size();
        std::vector<int> counts(nprocs), displs(nprocs, 0);
        MPI_Allgather(&nmine, 1, MPI_INT, counts.data(), 1, MPI_INT, lmp->world);
        for (int p = 1; p < nprocs; ++p) displs[p] = displs[p - 1] + counts[p - 1];
        std::vector<double> all(displs[nprocs - 1] + counts[nprocs - 1]);
        MPI_Allgatherv(mine.data(), nmine, MPI_DOUBLE, all.data(), counts.data(),
                       displs.data(), MPI_DOUBLE, lmp->world);
        double maxsq = 0.0;
        for (size_t i = 0; i < all.size(); i += 3)
            for (size_t j = i + 3; j < all.size(); j += 3) {
                double dx = all[i] - all[j], dy = all[i + 1] - all[j + 1],
                       dz = all[i + 2] - all[j + 2];
                lmp->domain->minimum_image(FLERR, dx, dy, dz);
                maxsq = std::max(maxsq, dx * dx + dy * dy + dz * dz);
            }
        return std::sqrt(maxsq);
    }

    double reaction_count()
    {
        lmp->input->one("variable nrxn equal f_myrxn[1]");
        int ivar = lmp->input->variable->find("nrxn");
        double count = lmp->input->variable->compute_equal(ivar);
        lmp->input->one("variable nrxn delete");
        return count;
    }

    Molecule *find_molecule(const char *name)
    {
        int idx = lmp->atom->find_molecule(name);
        return (idx >= 0) ? lmp->atom->molecules[idx] : nullptr;
    }
};

// ---------------------------------------------------------------------
// The pre-reaction template is inserted as-is, with Rmax set to
// comfortably exceed the initiator-atom separation baked into the
// template's own Coords (auto-detected -- no per-case tuning), so the
// initiator atoms are always within range. This test is purely about
// TopologyMatcher's own matching logic -- not Rmin/Rmax cutoff
// enforcement, which is a separate, already-well-covered concern this
// harness isn't trying to re-test -- so what "correct" means depends on
// the example, via ReactionTestCase::expect_reaction:
//
//  - ordinary example (expect_reaction == true, the default): the
//    reaction must fire exactly once, and the resulting topology
//    (bonds, angles, dihedrals, impropers, and every one of their
//    types, plus atom types) must be isomorphic to the post-reaction
//    template -- derived directly from the template file, not retyped
//    in this test, and tolerant of any genuinely symmetric atoms (see
//    file header).
//  - counterexample (the reaction lives under a top-level NEGATIVES/
//    folder): being within range is not sufficient for a match -- e.g.
//    a map-file constraint that's never satisfied, or atom types that
//    don't qualify as bonding partners -- so the reaction must fire
//    *zero* times, and the topology must be isomorphic to the unchanged
//    *pre*-reaction template. This is what lets a new example assert
//    "TopologyMatcher correctly declines this" instead of only ever
//    being able to assert "matches."
//
// A counterexample can instead be a *custom-system* reaction
// (tc.custom_system -- see "Custom systems for counterexamples" in
// bond-react-topology/README.md), for the case above's assumption
// doesn't hold: sometimes a candidate match is only ambiguous/invalid
// because of pre-existing topology the abstract pre-reaction template
// alone can't represent (e.g. two candidate insertion points that
// physically overlap on a shared backbone). There, "the pre-reaction
// template" isn't the system at all -- tc.system_script_file built a
// real, different system instead -- so "unchanged" means the *live
// system's own* topology from just before `run 1`, snapshotted here,
// not the template's.
// ---------------------------------------------------------------------
TEST_P(FixBondReactTopologyTest, MatchesExpectedOutcome)
{
    const ReactionTestCase &tc = GetParam();

    // A custom system exists specifically to exercise a match that a
    // verbatim template insertion can never produce (see above), so it
    // only makes sense as a counterexample; an ordinary (expect_reaction
    // == true) custom-system reaction would need an independent
    // "expected post-reaction system" this harness has no way to check
    // against, and is presumably a misplaced file rather than an
    // intentional new case.
    if (tc.custom_system && tc.expect_reaction) {
        GTEST_SKIP() << "'" << tc.name << "' has a custom system (" << tc.system_script_file
                     << ") but isn't under NEGATIVES/ -- custom systems are only supported "
                        "for counterexamples; see \"Custom systems for counterexamples\" in "
                        "bond-react-topology/README.md";
    }

    double rmin, rmax;
    if (tc.custom_system) {
        rmin = tc.custom_rmin;
        rmax = tc.custom_rmax;
    } else {
        double gap = initiator_gap_distance(tc);
        rmin = 0.0;
        rmax = gap * 1.25 + 0.05;
    }
    build_system(tc, rmin, rmax);

    TopologyGraph expected;
    if (tc.custom_system) {
        // The live system's own starting topology, not a template's --
        // see the comment above TEST_P.
        expected = graph_from_live_atoms(lmp);
    } else {
        Molecule *expected_tmpl = find_molecule(tc.expect_reaction ? "mol_post" : "mol_pre");
        ASSERT_NE(expected_tmpl, nullptr);
        expected = graph_from_molecule(expected_tmpl);
    }

    lmp->input->one("run 1");

    EXPECT_DOUBLE_EQ(reaction_count(), tc.expect_reaction ? 1.0 : 0.0);

    TopologyGraph actual = graph_from_live_atoms(lmp);
    std::string why;
    EXPECT_TRUE(topology_isomorphic(actual, expected, &why)) << why;
}

INSTANTIATE_TEST_SUITE_P(BondReactLibrary, FixBondReactTopologyTest,
                         ::testing::ValuesIn(kReactionLibrary),
                         [](const testing::TestParamInfo<ReactionTestCase> &info) {
                             return info.param.name;
                         });

// GoogleMock's stock main() never calls MPI_Init(), which a real MPI
// library rejects outright (OpenMPI aborts on the first MPI call), so
// this test brings its own -- the same one the other MPI-aware LAMMPS
// unit tests use; it also routes gtest output through rank 0 only.
#include "../testing/test_mpi_main.h"
