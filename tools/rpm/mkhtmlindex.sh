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
    <table>
    <tbody>
      <tr><th>[DIR]</th><th></th><th align="left"> &nbsp; <a href="../index.html">(Up one level)</a></th></tr>
EOF
    # first all directory entries
    for e in $(/bin/ls -1 "$d")
      do [ "$e" = "repodata" ] && continue
      [ "$e" = "drpms" ] && continue
      if [ -d "$d/$e" ]
      then
        echo '<tr><td>[DIR]</td><td></td>' >> "$d/index.html"
        echo "<td> &nbsp; <a href=\"$e/index.html\"> $e</a></td></tr>" >> "$d/index.html"
      fi
    done
    # then all file entries
    for e in $(/bin/ls -1 "$d")
      do [ "$e" = "index.html" ] && continue
      if [ -f "$d/$e" ]
      then
        echo '<tr><td>' >> "$d/index.html"
        ls -lh "$d/$e" | cut -d \  -f  6-7 >> "$d/index.html"
        echo "</td><td> &nbsp; <a href=\"$e\"> $e</a></td></tr>" >> "$d/index.html"
      fi
    done
    # footer
    cat >> "$d/index.html" <<EOF
    </table>
  </body>
</html>
EOF
done

