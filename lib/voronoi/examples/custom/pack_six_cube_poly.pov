#version 3.6;

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
light_source{<25,-12> color rgb <0.43,0.45,0.45>}

// The radius of the cylinders to be used when drawing the Voronoi cells
#declare r=0.08;

// The polydisperse particle packing
union{
#include "pack_six_cube_poly_p.pov"
	pigment{rgb <0.92,0.65,1>} finish{reflection 0.17 specular 0.3 ambient 0.42}
}

// The Voronoi cells for the packing, computed using the radical Voronoi
// tessellation
union{
#include "pack_six_cube_poly_v.pov"
	pigment{rgb <0.7,0.95,1>} finish{specular 0.5 ambient 0.42}
}
