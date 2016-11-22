
proc lerpcolor { col1 col2 alpha } {
  set dc [vecsub $col2 $col1]
  set nc [vecadd $col1 [vecscale $dc $alpha]]
  return $nc  
}

proc coltogs { col } {
  foreach {r g b} $col {}
  set gray [expr ($r + $g + $b) / 3.0]
  return [list $gray $gray $gray]
}

proc jet_tricolor_scale {} {
  display update off
  set mincolorid [expr [colorinfo num] - 1]
  set maxcolorid [expr [colorinfo max] - 1]
  set colrange [expr $maxcolorid - $mincolorid]
  set colhalf [expr $colrange / 2]
  for {set i $mincolorid} {$i < $maxcolorid} {incr i} {
    set colpcnt [expr ($i - $mincolorid) / double($colrange)]

    # The following color definitions for "jet" sort-of came from:
    #http://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale
    # but it was missing "green", so I inserted a green somewhere in the middle.

    # darkblue/violet 0.0
    set color0 { 0.08 0.0 0.77 }

    # blue 0.19
    set color1 { 0.0 0.0 1.0 }

    # cyan 0.34
    set color2 { 0.0 1.0 1.0 }

    # turquoise 0.4001
    set color3 { 0.0 1.0 0.78 }

    # green 0.445
    set color4 { 0.0 1.0 0.0 }

    # chartreuse 0.535
    set color5 { 0.875 1.0 0.0 }

    # yellow 0.69
    set color6 { 1.0 1.0 0.0 }

    # orange 0.73
    set color7 { 1.0 0.25 0.0 }

    # red 0.755
    set color8 { 1.0 0.0 0.0 }

    # darkred 1.0
    set color9 { 0.93 0.0 0.0 }

    if { $colpcnt < 0.19 } {
      set nc [lerpcolor $color0 $color1 [expr $colpcnt/(0.19-0.0)]]
    } elseif { $colpcnt < 0.34 } {
      set nc [lerpcolor $color1 $color2 [expr ($colpcnt-0.19)/(0.34-0.19)]]
    } elseif { $colpcnt < 0.4001 } {
      set nc [lerpcolor $color2 $color3 [expr ($colpcnt-0.34)/(0.4001-0.34)]]
    } elseif { $colpcnt < 0.445 } {
      set nc [lerpcolor $color2 $color3 [expr ($colpcnt-0.4001)/(0.445-0.4001)]]
    } elseif { $colpcnt < 0.535 } {
      set nc [lerpcolor $color3 $color4 [expr ($colpcnt-0.445)/(0.535-0.445)]]
    } elseif { $colpcnt < 0.69 } {
      set nc [lerpcolor $color4 $color5 [expr ($colpcnt-0.535)/(0.69-0.535)]]
    } elseif { $colpcnt < 0.73} {
      set nc [lerpcolor $color5 $color6 [expr ($colpcnt-0.69)/(0.73-0.69)]]
    } elseif { $colpcnt < 0.755} {
      set nc [lerpcolor $color6 $color7 [expr ($colpcnt-0.73)/(0.755-0.73)]]
    } else {
      set nc [lerpcolor $color7 $color8 [expr ($colpcnt-0.755)/(1.0-0.755)]]
    }

    #    set nc [coltogs $nc]
    foreach {r g b} $nc {}
    puts "index: $i $r $g $b   -- $colpcnt"
    display update ui
    color change rgb $i $r $g $b 
  }
  display update on
}

jet_tricolor_scale

