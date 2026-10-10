# fix bond/react topology unit test

Exercises `TopologyMatcher` (via `fix bond/react`) end-to-end: for each
reaction in this folder, it inserts the real pre-reaction molecule
template into a live simulation, lets one reaction fire, and checks the
resulting bonds/angles/dihedrals/impropers (and their types) directly
against the real post-reaction template -- read through LAMMPS's own
`Molecule` object, not retyped by hand. It complements
`unittest/force-styles`, which checks force-field numerics, not topology
bookkeeping.

The library is entirely file-driven: `discover_reaction_library()` in
`../test_fix_bond_react_topology.cpp` scans this folder at test-startup
time and builds one `TEST_P` case per reaction found. **Adding, removing,
or changing an example never requires touching the C++.** This file
documents that file-driven contract; it deliberately does not describe
which examples currently exist, since that's expected to keep changing.

## Adding an example

1. Make a subfolder anywhere under `bond-react-topology/` (any name --
   it becomes part of the test name; see "Test naming" below).
2. Put a pre-reaction molecule template, a post-reaction molecule
   template, and a `fix bond/react` map file directly inside it, named
   with `.pre`, `.post`, and `.map` extensions and **sharing the exact
   same filename stem** -- e.g. `reaction.pre`/`reaction.post`/
   `reaction.map`. All three are required. If a source example uses a
   different name per file (e.g. `{id}_unreacted.data`/
   `{id}_reacted.data`/`{id}.map`, as LAMMPS's own bundled REACTION
   examples do), rename the `.pre`/`.post` files to match the `.map`
   file's stem -- this harness matches by exact stem only, it does not
   guess at a looser match.
3. That's it. The next test run picks it up automatically:
   - Bond/angle/dihedral/improper/atom type counts, and whether each
     category uses plain numeric types or [type
     labels](https://docs.lammps.org/Howto_type_labels.html), are read
     directly from the `.pre`/`.post` files (see "What's auto-detected"
     below).
   - `Rmin`/`Rmax` for `fix bond/react` are derived from the initiator
     atoms' own separation in the `.pre` template's `Coords` (via the
     `.map` file's `InitiatorIDs`) -- no manual tuning needed.
   - The reaction is expected to fire exactly once, and the resulting
     topology must be isomorphic to the `.post` template (see "How a
     reaction is checked" below for what "isomorphic" means here and
     why).
4. A folder can hold more than one reaction: drop in another
   `{stem2}.pre`/`{stem2}.post`/`{stem2}.map` triple (different stem,
   same folder) for a multi-step example. Each triple found becomes its
   own `TEST_P` case.

## Adding a counterexample

A counterexample is a reaction that must legitimately produce **zero**
reactions even though it's well within `fix bond/react`'s cutoff -- e.g.
a map-file constraint that's never satisfied, or atom types that don't
qualify as bonding partners.

Build it exactly like an ordinary example (a `.pre`/`.post`/`.map`
triple with a shared stem), but put its folder one level inside a
top-level `NEGATIVES/` folder instead of directly under
`bond-react-topology/`:

```
bond-react-topology/
└── NEGATIVES/
    └── YourCounterexampleName/
        ├── stem.pre
        ├── stem.post   (still required -- see below)
        └── stem.map
```

That's the only way to mark a counterexample -- there's no marker file
or flag. Everything found one level inside `NEGATIVES/` is a
counterexample (reaction must fire zero times, topology must come back
unchanged); everything else in the library is an ordinary example. A
`NEGATIVES/` folder can itself hold more than one reaction triple, same
as any other folder.

A `.post` file is still required even though the reaction should never
actually produce it: `fix bond/react` always takes a product-template
argument, and this harness still needs it to read type counts/labels
from. Write it as whatever the reaction *would* produce if it fired.

## Custom systems for counterexamples

An ordinary counterexample still works like every ordinary example
does: the system under test is just the `.pre` template itself,
inserted verbatim. That's enough for a rejection the template and map
file alone can express (a distance/angle/custom constraint, a type
mismatch) -- but not for a rejection that depends on *pre-existing*
topology the abstract template can't represent by itself (for example,
two candidate insertion points whose neighbor searches physically
overlap on the same real atoms, which `fix bond/react`'s own internal
consistency checks are meant to catch). Inserting the template as-is
can't reproduce that: it would just trivially match itself, with no
ambiguity left.

For that case, a reaction under `NEGATIVES/` can supply its own system
instead, with two more files alongside the usual triple:

```
bond-react-topology/
└── NEGATIVES/
    └── YourCounterexampleName/
        ├── stem.pre       (abstract stencil fix bond/react matches against)
        ├── stem.post
        ├── stem.map
        ├── stem.data      (the real system -- referenced by stem.system)
        ├── stem.system    (LAMMPS commands that build that real system)
        └── stem.range     ("Rmin Rmax" for this reaction)
```

- **`stem.system`** is a LAMMPS input-script fragment containing
  everything that would otherwise be auto-generated from the templates
  alone: `units`, `atom_style`, `boundary`, the box and its atoms (via
  `read_data`, or its own `create_box`/`create_atoms`), the
  pair/bond/angle/dihedral/improper styles, and `special_bonds`. It's
  sourced verbatim in place of all of that; the harness still issues the
  `molecule`, `fix ... bond/react ...`, and thermo commands itself.
  Reference any sibling file (like `stem.data`) as
  `${exampledir}/stem.data`, not a bare relative filename -- the test
  binary's own working directory isn't this folder, but the harness
  defines `${exampledir}` to this folder's path before sourcing the
  script.
- **`stem.range`** holds this reaction's own `Rmin`/`Rmax` as two
  whitespace-separated numbers (e.g. `0.0 1.5`), since there's no
  template geometry here to derive them from automatically.

Both files must be present together for a reaction to use a custom
system; a custom-system reaction is only meaningful as a counterexample
and is skipped at runtime (with an explanatory message) if found outside
`NEGATIVES/`. Its "topology comes back unchanged" check compares against
the live system's own starting topology (snapshotted right before the
reaction runs), not the `.pre` template -- the two aren't the same thing
here.

## How a reaction is checked

The test does not hand-author a data file duplicating a reaction's
topology, does not hand-type the expected result in C++, and does not
assume a fixed correspondence between simulation atom tags and template
atom IDs. Instead, for each reaction:

1. The `.pre` template is inserted directly via
   `create_atoms 0 single ... mol mol_pre ...` (or, for a custom-system
   counterexample, the real system is built from `stem.system` instead).
2. `fix bond/react` runs for one timestep.
3. The resulting topology is compared against the expected one (the
   `.post` template for an ordinary example; the pre-run topology,
   unchanged, for a counterexample) via a bounded search for **any**
   type- and connectivity-consistent relabeling between the live atoms
   and the expected template's atoms -- not a single hard-coded one.

The relabeling search matters because a reaction template can contain
genuinely graph-symmetric atoms (e.g. the hydrogens on a spectator
`-CH2-` or `-CH3` group). When that happens, there's no basis for
`TopologyMatcher` to prefer one symmetric candidate over another, so two
equally correct runs can associate a given real atom with different
template roles. A comparison that insists on one fixed labeling would
fail on an equally valid alternate outcome; the relabeling search
accepts either. It only ever pairs same-type atoms and prunes
incrementally as bonds become determined, so for realistic reaction
sizes it's effectively instant, and a fully asymmetric template (every
atom uniquely typed) degenerates to a direct check with no real
backtracking.

## What's auto-detected

Nothing about a reaction's shape is hand-specified anywhere in this
harness -- all of it is read from the `.pre`/`.post` files themselves:

- **Type counts and labels.** For each category (atom/bond/angle/
  dihedral/improper), the type column is read from both files. If every
  value found is a plain nonnegative integer, that's a numeric-typed
  category and the count is just the highest value seen. Otherwise
  every value is treated as a [type
  label](https://docs.lammps.org/Howto_type_labels.html): the distinct
  labels are assigned numeric types in first-seen order (pre-file
  first, then post), and the matching `labelmap` command is issued
  before either `molecule` command runs, since a molecule template file
  can't define its own label map the way a data file can.
- **`Rmin`/`Rmax`.** Derived from the initiator atoms' separation in the
  `.pre` template's own `Coords`, via the `.map` file's `InitiatorIDs`
  (except for a custom-system reaction, which supplies these explicitly
  via `stem.range` instead).
- **`extra/{bond,angle,dihedral,improper}/per/atom` and
  `extra/special/per/atom`.** `create_box` needs to know, per atom, how
  much extra topology-storage headroom to reserve beyond what the
  inserted template already uses -- get this wrong and LAMMPS fails
  outright with `ERROR: Molecule topology/atom exceeds system
  topology/atom`. This is computed by counting, for each category, the
  most times any single atom appears across that category's data in
  either template, plus a small safety margin (`extra/special/per/atom`
  is instead capped at the template's own atom count minus one, since
  special-neighbor density doesn't reduce to a simple per-category
  count).

`InitiatorIDs`, `Coords`, and `Equivalences` are always plain integers
regardless of whether other sections use type labels.

## Multiple reactions per folder and test naming

A folder is a unit of *organization*, not of *reaction count*: it can
hold any number of `{stem}.pre`/`{stem}.post`/`{stem}.map` triples, each
becoming its own `TEST_P` case. This is meant for a multi-step
real-world example where several distinct reactions naturally belong
together.

Every test's name is `<folder>_<stem>` (e.g. filter with
`ctest -R FixBondReactTopology`, or `--gtest_filter=*<folder>*` on the
built binary directly), with a `NEGATIVES_` prefix added for anything
under `NEGATIVES/`. This naming is unconditional -- it doesn't matter
whether a folder holds one reaction or several -- so a test's name never
shifts just because a sibling reaction is later added to or removed from
its folder. The stem itself can be anything; it doesn't need to relate
to the folder name in any way. Non-identifier characters in either the
folder name or the stem are replaced with `_`.

## MPI

The test runs on any number of ranks, and is registered with ctest both
serially and (on MPI builds) on 3 and 4 ranks -- see
`CMakeLists_bond_react_topology_snippet.cmake`. It brings its own
`main()` (`unittest/testing/test_mpi_main.h`) because GoogleMock's stock
`main()` never calls `MPI_Init()`, and it prints results from rank 0 only.

- The live topology is gathered from all ranks before comparison, and
  the comparison runs identically on every rank.
- The inserted template is centered on the origin, so with 4 ranks it
  straddles subdomain boundaries; with 3 ranks it sits entirely inside
  a rank other than 0. Both layouts exercise different parallel paths
  in `fix bond/react`.
- The ghost cutoff is widened to span a whole template (and, for a
  custom system, the live system itself), since `fix bond/react` needs
  every atom of a candidate site visible on the rank doing the matching.

## Scope boundaries

- The reaction must create and delete no atoms -- every atom in the
  `.post` template must also appear in the `.pre` template (per the
  `.map` file's `Equivalences` section, which must still be authored
  correctly for the simulation itself to behave correctly, even though
  this test's own verification doesn't read it directly). Reactions
  using `CreateIDs`/`DeleteIDs` aren't supported.
- Before backtracking, the relabeling search narrows each atom's
  candidates with a Weisfeiler-Leman-style color-refinement pass over
  bond adjacency (not just raw atom type), so a large same-type group in
  a real molecule (e.g. many tens of chemically identical aromatic CH
  atoms) is usually collapsed down to just the atoms that are also
  genuinely symmetric by bonding position. This refinement is sound --
  it can only split a same-type group into smaller pieces, never
  wrongly exclude a genuine candidate -- but it isn't a full graph
  isomorphism solver, so a pathological template with very deep,
  perfectly repeating symmetry could in principle still leave a large
  group for the backtracking to branch on.
