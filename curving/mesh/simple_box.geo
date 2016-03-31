// Inputs
side = 1; 
height = side ; 
cl = side / 10;

// Geometry
Point(1) = {0, 0, 0, cl};
Point(2) = {side, 0, 0, cl};
Point(3) = {side, side, 0, cl};
Point(4) = {0, side, 0, cl};
Line(1) = {1, 2};				// bottom line
Line(2) = {2, 3};				// right line
Line(3) = {3, 4};				// top line
Line(4) = {4, 1};				// left line
Line Loop(5) = {1, 2, 3, 4}; 	
Plane Surface(6) = {5};

// Transfinite Line
Transfinite Line {1,2,3,4} = n_layer+1;

//Transfinite surface:
Transfinite Surface {6};
Recombine Surface {6};

surf[] = Extrude {0, 0, height} {
Surface{6};
Layers{n_layer};
Recombine;
};
// surfaceVector contains in the following order:
// [0]	- front surface (opposed to source surface)
// [1] - extruded volume
// [2] - bottom surface (belonging to 1st line in "Line Loop (6)")
// [3] - right surface (belonging to 2nd line in "Line Loop (6)")
// [4] - top surface (belonging to 3rd line in "Line Loop (6)")
// [5] - left surface (belonging to 4th line in "Line Loop (6)") 
// Physical Surface("front") = surfaceVector[0];
// Physical Volume("internal") = surfaceVector[1];
// Physical Surface("bottom") = surfaceVector[2];
// Physical Surface("right") = surfaceVector[3];
// Physical Surface("top") = surfaceVector[4];
// Physical Surface("left") = surfaceVector[5];
// Physical Surface("back") = {6};
Physical Surface(100) = {surf[0], surf[2], surf[3], surf[4], surf[5], 6};
Physical Volume(200) = surf[1];
