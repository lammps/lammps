#!/bin/bash

if [ $# -lt 1 ]
then
    echo "usage: $0 <list of directory trees to be indexed>"
    exit 1
fi

for d in $(find "$@" -type d -print)
do \
    dir=$(basename $d)
    # index file header
    cat > "$d/index.html" <<EOF
<html>
  <head><title>LAMMPS-ICMS RPM Repository: $d</title></head>
  <body style="font-family: mono, fixed, sans; margin-bottom: 2%; margin-top: 2%; margin-left: 5%; margin-right: 5%;">
    <h1>LAMMPS-ICMS RPM Repository: $d</h1>
    <h2>Contents of $dir</h2>
    <hr>
    <ui>
      <li> <a href="../index.html">[DIR] (Up one level)</a></li>
EOF
    # first all directory entries
    for e in $(/bin/ls -1 "$d")
    do \
        if [ -d "$d/$e" ]
          then
          echo '<li>[DIR]' >> "$d/index.html"
          ls -lhd "$d/$e" | cut -d \  -f  6-8 >> "$d/index.html"
          echo "<a href=\"$e/index.html\">$e</a></li>" >> "$d/index.html"
      fi
    done
    # then all file entries
    for e in $(/bin/ls -1 "$d")
    do \
        if [ -f "$d/$e" ]
          then
          echo '<li>' >> "$d/index.html"
          ls -lhd "$d/$e" | cut -d \  -f  5-8 >> "$d/index.html"
          echo "<a href=\"$e\">$e</a></li>" >> "$d/index.html"
      fi
    done
    # footer
    cat >> "$d/index.html" <<EOF
    <ui>
  </body>
</html>
EOF
done

