#!/bin/bash

if [ $# -lt 1 ]
then
    echo "usage: $0 <list of directory trees to be indexed>"
    exit 1
fi

for d in $(find "$@" -type d -print)
do \
    dir=$(basename $d)
    uplink=index.html
    [ $dir == windows ] && uplink=windows.html
    # index file header
    cat > "$d/index.html" <<EOF
<html>
  <head><title>LAMMPS-ICMS Binaries Repository: $d</title>

    <script type="text/javascript">
      var _gaq = _gaq || [];
      _gaq.push(['_setAccount', 'UA-12573382-3']);
      _gaq.push(['_trackPageview']);

      (function() {
      var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
      ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
      var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
      })();
    </script>


  </head>
  <body style="font-family: mono, fixed, sans; margin-bottom: 2%; margin-top: 2%; margin-left: 5%; margin-right: 5%;">
    <h1>LAMMPS-ICMS Binaries Repository: $d</h1>
    <h2>Contents of $dir</h2>
    <hr>
    <table>
    <tbody>
      <tr><th>[DIR]</th><th></th><th align="left"> &nbsp; <a href="../${uplink}">(Up one level)</a></th></tr>
EOF
    # first all directory entries
    for e in $(/bin/ls -1 "$d")
      do [ "$e" = "repodata" ] && continue
      [ "$e" = "drpms" ] && continue
      [ "$e" = "cache" ] && continue
      if [ -d "$d/$e" ]
      then
        echo '<tr><td>[DIR]</td><td></td>' >> "$d/index.html"
        echo "<td> &nbsp; <a href=\"$e/index.html\"> $e</a></td></tr>" >> "$d/index.html"
      fi
    done
    # then all file entries
    for e in $(/bin/ls -1r "$d")
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

