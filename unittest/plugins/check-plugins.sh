#!/bin/sh
# Load plugins into LAMMPS and check that each registers the styles that its source file
# registers (plugin.style and plugin.name in lammpsplugin_init()).  A plugin that cannot be
# loaded only gets a warning from LAMMPS, so its styles are missing from "plugin list".  A
# plugin file that does not exist (e.g. after a failed compilation) is reported, and the others
# are still loaded.
#
# Usage: check-plugins.sh <lmp executable> <plugin source file> <plugin file> [<source> <plugin> ...]
# e.g.   check-plugins.sh build/lmp examples/PACKAGES/pace/plugin/paceplugin.cpp build/paceplugin/paceplugin.so

if [ $# -lt 3 ] || [ $(($# % 2)) -ne 1 ]; then
  echo "usage: $0 <lmp executable> <plugin source file> <plugin file> [<source> <plugin> ...]" >&2
  exit 2
fi
lmp=$1
shift
tmp=$(mktemp -d) || exit 2
trap 'rm -rf "$tmp"' EXIT

files=
: > "$tmp/expected"
while [ $# -gt 0 ]; do
  [ -f "$1" ] || { echo "no plugin source file $1" >&2; exit 2; }
  awk -F'"' '/plugin\.style *=/ {style = $2} /plugin\.name *=/ {print style, $2}' "$1" >> "$tmp/expected"
  if [ -f "$2" ]; then
    files="$files $2"
  else
    echo "no plugin file $2" >&2
  fi
  shift 2
done

: > "$tmp/output"
if [ -n "$files" ]; then
  printf 'plugin load%s\nplugin list\n' "$files" | "$lmp" -nocite -log none -echo screen > "$tmp/output" 2>&1
  status=$?
  cat "$tmp/output"
  if [ $status -ne 0 ]; then
    echo "LAMMPS failed with exit status $status" >&2
    exit 1
  fi
fi
sed -n 's/^ *[0-9]*: \(.*\) style plugin \(.*\)$/\1 \2/p' "$tmp/output" | sort > "$tmp/loaded"
sort -o "$tmp/expected" "$tmp/expected"
if ! cmp -s "$tmp/expected" "$tmp/loaded"; then
  echo "styles registered by the sources (<) and loaded (>) differ:" >&2
  diff "$tmp/expected" "$tmp/loaded" >&2
  exit 1
fi
echo "all $(wc -l < "$tmp/loaded") styles of the plugins are loaded"
