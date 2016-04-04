// Inputs
ratio=6;
n_r = 60;
n_angle=n_r/ratio;
r_in = 1;
r_out = 2;
stet=0.86602540378;
ctet=0.5;
prog=1.1;

// Geometry
Point(1) = {0, r_in, 0};
Point(2) = {0, r_out, 0};
Point(3) = {0, r_out*ctet, r_out*stet};
Point(4) = {0, r_in*ctet, r_in*stet};
Point(5) = {0, 0, 0};

Line(1) = {1,2};
Circle(2) = {2, 5, 3};
Line(3) = {3, 4};
Circle(4) = {4, 5, 1};

Line Loop(5) = {-1, -2, -3, -4}; 	
Ruled Surface(6) = {5};


// Surface 6 structured

Transfinite Line {2, 4} = n_angle;
Transfinite Line {1} = n_r Using Progression prog;
Transfinite Line {3} = n_r Using Progression 1/prog;
Transfinite Surface {6};
Recombine Surface {6};

surf[] = Extrude{{0,0,1}, {0,0,0}, -3.14159265359/2}{
Surface{6};
Layers{n_angle};
Recombine;
};

//
// surfaceVector contains in the following order:
// [0]	- front surface (opposed to source surface)
// [1] - extruded volume
// [2] - bottom surface (belonging to 1st line in "Line Loop (6)")
// [3] - right surface (belonging to 2nd line in "Line Loop (6)")
// [4] - top surface (belonging to 3rd line in "Line Loop (6)")
// [5] - left surface (belonging to 4th line in "Line Loop (6)") 
Physical Surface(2) = surf[0];
Physical Volume(1) = surf[1];
Physical Surface(3) = surf[2];
Physical Surface(5) = surf[3];
Physical Surface(6) = surf[4];
Physical Surface(4) = surf[5];
Physical Surface(1) = {6};

