#version 3.6;

// Right-handed coordinate system in which the z-axis points upwards
camera {
	location <20,-50,30>
	sky z
	right -0.3*x*image_width/image_height
	up 0.3*z
	look_at <0,0,3>
}

// White background
background{rgb 1}

// Two lights with slightly different colors
light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<20,-15,5> color rgb <0.38,0.40,0.40>}

// Radius of the Voronoi cell network, and the particle radius
#declare r=0.025;
#declare s=0.1;

// Particles
union{
#include "frustum_p.pov"
	scale 10
	pigment{rgb <0.4,0.85,0.95>} finish{reflection 0.1 specular 0.3 ambient 0.42 metallic}
}

// Voronoi cells
union{
#include "frustum_v.pov"
	scale 10
	pigment{rgb <0.5,0.5,0.51>} finish{specular 0.3 ambient 0.42 reflection 0.4 metallic}
}
