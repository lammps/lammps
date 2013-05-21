#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 64
}

camera {
	location <30,-50,25>
	sky z
	up 0.15*z
	right -0.15*x*image_width/image_height
	look_at <0,0,2.8>
}

background{rgb 1}

light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-12,12> color rgb <0.43,0.45,0.45>}

#declare r=0.06;
#declare rr=0.08;

#declare f1=finish{reflection 0.15 specular 0.3 ambient 0.42}

#declare t1=texture{pigment{rgbft <0.3,0.7,0.9,0,0.4>} finish{f1}}
#declare t2=texture{pigment{rgb <0.3,0.4,0.9>} finish{f1}}

union{
#include "polygons6_v.pov"
}

union{
#include "polygons_d.pov"
	texture{T_Silver_4B}
}
