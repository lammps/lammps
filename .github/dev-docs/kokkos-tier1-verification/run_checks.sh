#!/bin/bash
# Drive the KOKKOS sync detectors over the Tier-1 styles.
#   ./run_checks.sh step0     thermo: plain CPU styles vs KOKKOS styles, same binary
#   ./run_checks.sh detect    watch/stale detectors on the KOKKOS runs
# Every KOKKOS run uses the settings a GPU would pick, per
# .github/dev-docs/kokkos-sync-debugging.md.
set -u
D=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$D/../../.." && pwd)
LMP=${LMP:-$ROOT/build-sync/lmp}
export LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS:-$ROOT/potentials}
DATA="-var tdir $ROOT/unittest/force-styles/tests"
OUT=${OUT:-$D/out}
MPIRUN=${MPIRUN:-"mpirun --allow-run-as-root --oversubscribe"}
GPUPK="neigh full newton off comm device sort device atom/map device gpu/aware on"
mkdir -p "$OUT"

# deck|variable assignments (empty for none)|label
CASES=$(cat <<'EOC'
in.fep|ps lj/cut/soft
in.fep|ps coul/cut/soft
in.fep|ps lj/cut/coul/cut/soft
in.fep|ps lj/class2/soft
in.fep|ps morse/soft
in.fep|ps coul/long/soft
in.fep|ps lj/cut/coul/long/soft
in.fep|ps lj/charmm/coul/long/soft
in.cs|ps coul/long/cs
in.cs|ps born/coul/long/cs
in.cs|ps buck/coul/long/cs
in.cs|ps lj/cut/coul/long/cs
in.cs|ps lj/class2/coul/long/cs
in.cs|ps coul/wolf/cs
in.cs|ps born/coul/wolf/cs
in.cs|ps born/coul/dsf
in.cs|ps born/coul/dsf/cs
in.single|ps lj/smooth/linear
in.single|ps lj/sf
in.single|ps nm/cut/split
in.single|ps coul/slater/cut
in.single|ps lj/relres
in.single|ps lj/charmmfsw/coul/charmmfsh
in.bonded|ts linear
in.bonded|ts spline
in.bonded_cut|ts linear
in.bonded_cut|ts spline
in.manybody|ps sw/mod pot Si.sw|neigh half newton on comm device sort device atom/map device gpu/aware on
in.manybody|ps tersoff/mod/c pot Si.tersoff.modc|neigh half newton on comm device sort device atom/map device gpu/aware on
in.sphere_pair|ps lj/expand/sphere
in.gran|ps gran/hooke|neigh half newton on comm device sort device atom/map device gpu/aware on
in.gran|ps gran/hertz/history|neigh half newton on comm device sort device atom/map device gpu/aware on
in.fixes|
in.fixes_sphere|
in.piston|
in.computes|
in.regions|
in.asphere|
in.min|ms sd
in.min|ms quickmin
in.min|ms fire
in.min|ms cg
EOC
)

vargs() { local v=($1); local o=""; local i=0
  while [ $i -lt ${#v[@]} ]; do o="$o -var ${v[$i]} ${v[$((i+1))]}"; i=$((i+2)); done
  echo "$o"; }

label() { echo "$1$(echo " $2" | tr ' /' '_-')"; }

step0() {
  printf '%-46s %-12s %s\n' STYLE VERDICT TOL_RATIO
  echo "$CASES" | while IFS='|' read -r deck vars pkov; do
    case "$deck" in ""|\#*) continue;; esac
    L=$(label "$deck" "$vars"); VA=$(vargs "$vars"); PK=${pkov:-$GPUPK}
    ( cd "$D" && $LMP -in "$deck" $DATA $VA -log none -screen "$OUT/$L.cpu" >/dev/null 2>&1 )
    ( cd "$D" && $LMP -in "$deck" $DATA $VA -log none -screen "$OUT/$L.kk" \
        -k on -sf kk -pk kokkos $PK >/dev/null 2>&1 )
    python3 "$D/cmp.py" "$OUT/$L.cpu" "$OUT/$L.kk" "$L"
  done
}

detect() {
  echo "$CASES" | while IFS='|' read -r deck vars pkov; do
    case "$deck" in ""|\#*) continue;; esac
    L=$(label "$deck" "$vars"); VA=$(vargs "$vars"); PK=${pkov:-$GPUPK}
    ( cd "$D" && LMP_KOKKOS_WATCH= LMP_KOKKOS_STALE= LMP_KOKKOS_STALE_STRICT=1 \
        $LMP -in "$deck" $DATA $VA -log none -screen "$OUT/$L.det" \
        -k on -sf kk -pk kokkos $PK > "$OUT/$L.det.err" 2>&1 )
    n=$(grep -cE '^\[(stale|watch)\]' "$OUT/$L.det.err" 2>/dev/null)
    fin=$(grep -c "Total wall time" "$OUT/$L.det" 2>/dev/null)
    printf '%-46s reports=%-5s finished=%s\n' "$L" "$n" "$fin"
  done
}


# same comparison under two MPI ranks, which exercises the exchange and border
# communication the serial pass never reaches
mpi() {
  printf '%-46s %-12s %s\n' STYLE VERDICT TOL_RATIO
  echo "$CASES" | while IFS='|' read -r deck vars pkov; do
    case "$deck" in ""|\#*) continue;; esac
    L=$(label "$deck" "$vars"); VA=$(vargs "$vars"); PK=${pkov:-$GPUPK}
    ( cd "$D" && $MPIRUN -np 2 $LMP -in "$deck" $DATA $VA -log none -screen "$OUT/$L.mcpu" >/dev/null 2>&1 )
    ( cd "$D" && $MPIRUN -np 2 $LMP -in "$deck" $DATA $VA -log none -screen "$OUT/$L.mkk" \
        -k on -sf kk -pk kokkos $PK >/dev/null 2>&1 )
    python3 "$D/cmp.py" "$OUT/$L.mcpu" "$OUT/$L.mkk" "$L"
  done
}

# watch/stale detectors under two MPI ranks
detect_mpi() {
  echo "$CASES" | while IFS='|' read -r deck vars pkov; do
    case "$deck" in ""|\#*) continue;; esac
    L=$(label "$deck" "$vars"); VA=$(vargs "$vars"); PK=${pkov:-$GPUPK}
    R=$(mktemp -d "${TMPDIR:-/tmp}/rcout.XXXXXX")
    ( cd "$D" && LMP_KOKKOS_WATCH= LMP_KOKKOS_STALE= LMP_KOKKOS_STALE_STRICT=1 \
        $MPIRUN -np 2 --output-filename "$R" $LMP -in "$deck" $DATA $VA -log none \
        -screen "$OUT/$L.mdet" -k on -sf kk -pk kokkos $PK >/dev/null 2>&1 )
    cat "$R"/1/rank.*/std* > "$OUT/$L.mdet.err" 2>/dev/null
    rm -rf "$R"
    n=$(grep -cE '^\[(stale|watch)\]' "$OUT/$L.mdet.err" 2>/dev/null)
    fin=$(grep -c "Total wall time" "$OUT/$L.mdet" 2>/dev/null)
    printf '%-46s reports=%-5s finished=%s\n' "$L" "$n" "$fin"
  done
}


# poison mode: the only detector that sees a read through a plain pointer.
# needs an executable built with -D KOKKOS_DEBUG_SYNC_ASAN=on, which LMP must
# point at (build-poison/lmp at the top of the repository by default).
poison() {
  echo "$CASES" | while IFS='|' read -r deck vars pkov; do
    case "$deck" in ""|\#*) continue;; esac
    L=$(label "$deck" "$vars"); VA=$(vargs "$vars"); PK=${pkov:-$GPUPK}
    ( cd "$D" && LMP_KOKKOS_POISON=1 ASAN_OPTIONS=detect_leaks=0 \
        ${POISON_LMP:-$ROOT/build-poison/lmp} -in "$deck" $DATA $VA -log none \
        -screen "$OUT/$L.psn" -k on -sf kk -pk kokkos $PK > "$OUT/$L.psn.err" 2>&1 )
    n=$(grep -c "ERROR: AddressSanitizer" "$OUT/$L.psn.err" 2>/dev/null)
    fin=$(grep -c "Total wall time" "$OUT/$L.psn" 2>/dev/null)
    printf '%-46s asan=%-5s finished=%s\n' "$L" "$n" "$fin"
  done
}

"$@"
