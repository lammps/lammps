#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 64
}

camera {
	location <1,2.5,-5>
	right 0.23*x*image_width/image_height
	up 0.23*y
	look_at <0,0,0>
}

background{rgb 1}

light_source{<-0.8,3,-2> color rgb <0.77,0.75,0.65>}
light_source{<2.5,1.2,-1.2> color rgb <0.38,0.30,0.40>}

#declare r=0.006;

union{
#include "superellipsoid_m.pov"
	texture{T_Gold_3C}
}

union{
#include "superellipsoid_v.pov"
	texture{T_Silver_3C}
}
