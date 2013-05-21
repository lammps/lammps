#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

camera {
	location <0,-50,15>
	sky z
	right -0.25*x*image_width/image_height
	up 0.25*y
	look_at <0,0,0>
}

background{rgb 1}

light_source{<-16,-30,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-16,8> color rgb <0.43,0.45,0.45>}

#declare r=0.05;
#declare s=0.5;

#declare particles=union {
#include "irregular_p.pov"
	pigment{rgb 0.9} finish{reflection 0.2 specular 0.25 ambient 0.28 metallic}
}

union{particles translate <-4,0,0>}
union{particles translate <4,0,0>}

union{
#include "irregular_v.pov"
	translate <4,0,0>
	pigment{rgb <0.7,0.3,0.9>} finish{specular 0.5 ambient 0.42}
}
