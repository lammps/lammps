#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 64
}

camera {
	location <15,10,-30>
	right 0.04*x*image_width/image_height
	up 0.04*y
	look_at <0,0,0>
}

background{rgb 1}

light_source{<-4,35,-20> color rgb <0.72,0.7,0.7>}
light_source{<20,5,-10> color rgb <0.4,0.42,0.42>}

#declare r=0.01;

union{
#include "degenerate_v.pov"
	texture{T_Gold_3C}
}

sphere{<0,0,0>,0.5-r
	pigment{rgb <0.1,0.15,0.4>}
	finish{reflection 0.25 specular 0.4 ambient 0.6}
}
