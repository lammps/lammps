#version 3.6;

// Increase the trace level for more accurate reflections
global_settings {
	max_trace_level 64
}

// Right-handed coordinate system in which the z-axis points upwards
camera {
	location <0,-40,34>
	sky z
	right -0.4*x*image_width/image_height
	up 0.4*z
	look_at <0,0,-1.31>
}

// White background
background{rgb 1}

// Two lights with slightly different colors
light_source{<-16,-20,43> color rgb <0.72,0.69,0.69>}
light_source{<30,-15,12> color rgb <0.34,0.37,0.37>}

// Radius of the Voronoi cell network, and the particle radius
#declare r=0.06;
#declare s=0.5;

// Particles
union{
	#include "torus_p.pov"
	pigment{rgb 0.97}
	finish{reflection 0.1 ambient 0.30 specular 0.3}
}

// Voronoi cells, using a radial pigment map
union{
	#include "torus_v.pov"
	pigment{radial pigment_map {
			[0 rgb <0.5,0.7,1>]
			[0.25 rgb <0.38,0.82,0.92>]
			[0.5 rgb <0.5,0.7,1>]
			[0.75 rgb <0.65,0.4,1>]
			[1 rgb <0.5,0.7,1>]}
		rotate <270,0,0>}
	finish{specular 0.3 ambient 0.3 reflection 0.1}
}
