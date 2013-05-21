#version 3.6;

#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 64
}

camera {
	location <15,-30,10>
	sky z
	right -0.033*x*image_width/image_height
	up 0.033*y
	look_at <0,0,0>
}

background{rgb 1}

light_source{<-4,-20,35> color rgb <0.72,0.7,0.7>}
light_source{<20,-10,5> color rgb <0.4,0.42,0.42>}

#declare r=0.01;

#declare f1=finish{reflection 0.17 specular 0.3 ambient 0.42}
#declare t_odd=texture{pigment{rgb <0.8,0.8,0.8>} finish{f1}}
#declare t_even=texture{pigment{rgb <0.2,0.2,0.2>} finish{f1}}

intersection{
#include "odd_even_pl.pov"
}

union{
#include "odd_even_v.pov"
	pigment{rgb <0.9,0.1,0.15>}
	finish{reflection 0.25 specular 0.4 ambient 0.3}
}
