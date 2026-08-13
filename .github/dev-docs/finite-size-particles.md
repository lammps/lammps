# Extended (Finite-Size) Particles and Inertia / Angular-Momentum Diagnostics

Lessons from making `compute inertia`, `compute inertia/chunk`, `compute angmom/chunk`,
`compute omega/chunk`, and the `inertia()`/`angmom()`/`omega()` variable functions
account for finite-size particles instead of treating every atom as a point mass.

## Where the shape data lives

**The extended atom styles are CORE, not in a package.** `src/atom_vec_ellipsoid.{h,cpp}`,
`atom_vec_line`, `atom_vec_tri`, and `atom_vec_body` live directly in `src/`, so any core
file may `#include` them and resolve a particle's shape data via
`dynamic_cast<AtomVec*>(atom->style_match("ellipsoid"|"line"|"tri"|"body"))` plus the
per-atom index arrays `atom->ellipsoid/line/tri/body` (>= 0 means a shape is assigned).
No package guard is needed (e.g. `set.cpp` already does this).  The body *styles*
(`body nparticle`, ...) are in the BODY package, but the atom-vec interface is core.

## Per-particle spin inertia

Use `MathExtra::inertia_*` (see `fix rigid`).  Reference implementation is
`src/RIGID/fix_rigid.cpp` (`setup_bodies_static` for inertia, `setup_bodies_dynamic`
for angular momentum).  Prefactors: sphere `SINERTIA = 2/5`, line `LINERTIA = 1/12`.
Helpers return the tensor as a 6-vector in **Voigt order `{Ixx,Iyy,Izz,Iyz,Ixz,Ixy}`**:

- `inertia_ellipsoid(shape, quat, mass, ivec)` -- regular ellipsoid from semi-axes.
- `inertia_ellipsoid(idiag, quat, ivec)` -- rotate stored principal moments to the
  space frame; used for **superellipsoids** (`bonus_super[].inertia`) and **body**
  particles (`bonus[].inertia`), which both store mass-weighted principal moments +
  quaternion.
- `inertia_line(length, theta, mass, ivec)`, `inertia_triangle(idiag, quat, mass, ivec)`.

**Voigt-order vs output-column mismatch is an easy bug.** `compute inertia/chunk` and
the `inertia()` variable order columns `{Ixx,Iyy,Izz,Ixy,Iyz,Ixz}`, but `MathExtra`
returns `{Ixx,Iyy,Izz,Iyz,Ixz,Ixy}`.  Map `col3+=ivec[5]` (Ixy), `col4+=ivec[3]` (Iyz),
`col5+=ivec[4]` (Ixz).  Off-diagonals are 0 for axis-aligned shapes, so a wrong mapping
only shows up once a body is rotated -- test with a non-identity quaternion.

## Dispatch order and storage traps

- **Check shape indices BEFORE the sphere (radius) branch.** Ellipsoid, superellipsoid,
  line, triangle, and body particles also carry a finite `atom->radius` (their bounding
  radius) -- since the spheroid support for lines/tris (Feb 2026), the line/tri atom
  vecs ALWAYS store it when bonus data is set -- so a leading `radius[i] > 0` test
  silently mis-handles them as spheres.  Dispatch on the bonus indices
  (`ellipsoid/line/tri/body[i] >= 0`) first and on `radius` last; spheroid particles
  have bonus index -1 with `radius > 0` and still classify as spheres.  Radius-first
  cascades are a real regressed bug class (they froze tri/line orientations in
  `fix rigid`/`fix rigid/small` and made `fix srd` do sphere collisions); grep for the
  pattern whenever touching finite-size dispatch.
- **Superellipsoids use `bonus_super`, not `bonus`.** `AtomVecEllipsoid` allocates
  exactly one of the two depending on `atom->superellipsoid_flag`; the regular `bonus`
  pointer is null in a superellipsoid system, so dereferencing it crashes.  Branch on
  the flag.
- **2d rigid-body frames must keep `ez = +z`.**  For 2d bodies the in-plane inertia
  moments are near-degenerate, so `MathExtra::jacobi3` eigenvector signs are
  solver/build dependent, and the right-handedness enforcement afterwards can leave
  `ez = -z`.  The scalar-theta orientation math for LINE particles assumes a pure +z
  rotation (a flipped frame turns the in-plane map into a reflection, giving wrong
  per-rod orientations); TRIANGLE quaternion math is immune.  After the handedness
  fix-up, if `ez[2] < 0` in 2d, negate BOTH `ey` and `ez` (still a right-handed
  eigenbasis) -- see `setup_bodies_static` in `fix rigid`/`fix rigid/small`.
- **Angular momentum: OMEGA-type vs ANGMOM-type.** Sphere and line particles carry an
  angular velocity (`atom->omega`); their spin angular momentum is `I_spin * omega`
  (sphere: all 3 components; line: z only).  Ellipsoid, superellipsoid, triangle, and
  body particles store their angular momentum directly in `atom->angmom` -- add it
  as-is.  Guard on the `atom->omega`/`atom->angmom` pointers (null for atom styles
  lacking them).

## Do NOT add finite-size terms to the shared Group::inertia / Group::angmom

Those methods are also called by physics routines that remove bulk rotation from the
*translational* velocities -- `fix momentum` (angular), `velocity ... zero angular`,
`compute temp/rotate`, `fix addtorque/group`, `fix tfmc`.  Those need the point-mass
(orbital) tensor/vector: they form `omega = I^-1 L` and subtract `omega x r` from `v`,
which is self-consistent only if both `I` and `L` are orbital.  Adding spin to the
shared method (while `L` stays orbital) makes `omega` wrong and breaks them for
finite-size systems.  Instead the diagnostics call the orbital method and then add an
**additive** finite-size term: `Group::inertia_extended()` / `Group::angmom_extended()`
(each with a region variant).  The spin contribution is center-of-mass-independent, so
it is a clean additive term.  Per-chunk computes inline the same per-particle
contribution.

## Test setup for finite-size particles (`unittest/commands/`)

Build tiny systems with `create_atoms ... single X Y Z units box` so the
inertia/angmom is analytic:

- `set shape Sx Sy Sz` stores **semi-axes = 0.5 * the given diameters** and initializes
  the quaternion to identity.
- `atom_style ellipsoid superellipsoid` + default blockiness `(2,2)` makes a
  superellipsoid reduce exactly to an ellipsoid (same analytic result), exercising the
  `bonus_super` path.  The stored inertia is computed **at `set shape` time**, so set
  the mass *before* the shape.
- `set ... omega Wx Wy Wz` (needs `omega_flag`, e.g. sphere) and `set ... angmom Lx Ly
  Lz` (needs `angmom_flag`, e.g. ellipsoid) seed the spin terms.
- Body particles need a `read_data` file with a `Bodies` section (`atom_style body
  nparticle 1 1`; one body, `ndouble = 6 + 3*nsub`, first 6 doubles are the inertia
  tensor that `data_body` diagonalizes).  A single body at the origin round-trips its
  stored inertia tensor through `compute inertia`.
- Point-particle regressions (e.g. the `fourmol` chunk-compute tests) must stay
  unchanged, which confirms finite-size handling is inert for point particles.
