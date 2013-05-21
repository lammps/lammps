#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 64
}

camera {
	location <20,30,-50>
	right 0.26*x*image_width/image_height
	up 0.26*y
	look_at <0,-2,0>
}

background{rgb 1}

light_source{<-8,30,-20> color rgb <0.77,0.75,0.75>}
light_source{<20,5,-15> color rgb <0.38,0.40,0.40>}

#declare r=0.025;
#declare s=0.07;

union{
#include "tetrahedron_p.pov"
	scale 5
	rotate <0,-65,0>
	texture{T_Silver_3C}
}

union{
#include "tetrahedron_v.pov"
	scale 5
	rotate <0,-65,0>
	pigment{rgb <0.3,0.3,0.9>} finish{phong 0.9 ambient 0.42 reflection 0.1}
}
