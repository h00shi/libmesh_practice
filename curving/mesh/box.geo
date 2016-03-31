// Inputs
side = 200; 
height = side ; 
cl = side / 10;
n_layer = 30;

// Geometry
Point(1) = {-side/2, -side/2, 0, cl};
Point(2) = {side/2, -side/2, 0, cl};
Point(3) = {side/2, side/2, 0, cl};
Point(4) = {-side/2, side/2, 0, cl};
Line(1) = {1, 2};				// bottom line
Line(2) = {2, 3};				// right line
Line(3) = {3, 4};				// top line
Line(4) = {4, 1};				// left line
Line Loop(5) = {1, 2, 3, 4}; 	
Plane Surface(6) = {5};

// Transfinite Line
Transfinite Line {1,2,3,4} = n_layer Using Bump 1/10;

//Transfinite surface:
Transfinite Surface {6};
Recombine Surface {6};

surfaceVector[] = Extrude {0, 0, height} {
Surface{6};
//Layers{{n_layer/3, n_layer/3, n_layer/3}, {0.2, 0.8, 1}};
Layers{n_layer} Using Bump 1/10;
Recombine;
};
// surfaceVector contains in the following order:
// [0]	- front surface (opposed to source surface)
// [1] - extruded volume
// [2] - bottom surface (belonging to 1st line in "Line Loop (6)")
// [3] - right surface (belonging to 2nd line in "Line Loop (6)")
// [4] - top surface (belonging to 3rd line in "Line Loop (6)")
// [5] - left surface (belonging to 4th line in "Line Loop (6)") 
Physical Surface("front") = surfaceVector[0];
Physical Volume("internal") = surfaceVector[1];
Physical Surface("bottom") = surfaceVector[2];
Physical Surface("right") = surfaceVector[3];
Physical Surface("top") = surfaceVector[4];
Physical Surface("left") = surfaceVector[5];
Physical Surface("back") = {6};
