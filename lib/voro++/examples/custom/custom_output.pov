#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

// Increase tracing level for more accurate reflections
global_settings {
	max_trace_level 64
}

// Right-handed coordinate system where the z-axis points upwards
camera {
	location <30,-50,25>
	sky z
	right -0.15*x*image_width/image_height
	up 0.15*z
	look_at <0,0,2.8>
}

// White background
background{rgb 1}

// Two lights with slightly different colors
light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-12,12> color rgb <0.43,0.45,0.45>}

// The radius of the cylinders to be used when drawing the Voronoi cells
#declare r=0.06;

// The radius of the particles
#declare s=0.5;

// Different colors for the particles
#declare f1=finish{reflection 0.17 specular 0.3 ambient 0.42}
#declare t8=texture{pigment{rgb <0.6,0.6,0.6>} finish{f1}}
#declare t9=texture{pigment{rgb <0.7,0.3,1>} finish{f1}}
#declare t10=texture{pigment{rgb <0.3,0.7,1>} finish{f1}}
#declare t11=texture{pigment{rgb <0.2,0.9,0.9>} finish{f1}}
#declare t12=texture{pigment{rgb <0.3,1,0.7>} finish{f1}}
#declare t13=texture{pigment{rgb <0.7,1,0.3>} finish{f1}}
#declare t14=texture{pigment{rgb <0.9,0.9,0.2>} finish{f1}}
#declare t15=texture{pigment{rgb <1,0.7,0.3>} finish{f1}}
#declare t16=texture{pigment{rgb <1,0.3,0.7>} finish{f1}}
#declare t17=texture{pigment{rgb <0.9,0.2,0.9>} finish{f1}}
#declare t18=texture{pigment{rgb <1,1,1>} finish{f1}}

// The particle packing
union{
#include "custom_output_p.pov"
}

// The computed Voronoi cells
union{
#include "pack_six_cube_v.pov"
	texture{T_Silver_4B}
}
