/****************************************************************************
 *                                    Input                                 *
 ****************************************************************************/
v_max = Pi/6;
r0 = 1;
h0 = r0 * Cos(v_max);
l = 5 * r0 * Sin(v_max);
w = l/2;
g = l/3;
b= 2;
c= 4;

/****************************************************************************
 *                           Intermediate Values                            *
 ****************************************************************************/
// back
u_back= 0;
ux_back= l/2;
uy_back= 0;
r_back= r0;
h_back= h0;
p_back= Sqrt(r_back*r_back - h_back * h_back);

// front
u_front= w/c;
ux_front= l/2 + w/c;
uy_front= w;
r_front= r0 - b * w / c;
h_front= h0 - h0/r0*b * w / c;
p_front= Sqrt(r_front*r_front - h_front * h_front);


/****************************************************************************
 *                                Size Control                              *
 ****************************************************************************/
cl = 1e-1;

//nx1 = 8;
progx1 = 1;

nx2 = nx1 / (l/2 - p_back) * 2 * p_back;
If (nx2 < 3)
nx2 = 3;
EndIf
progx2 = 1;

nx3 = nx1 / 1.5;
progx3 = 1;
If (nx3 < 3)
nx3 = 3;
EndIf

ny = nx1/2 + 2;
If (ny > 8)
ny = 8;
EndIf

nz = nx1*3;
If (nz > 20)
nz = 20;
EndIf
progz = 1.4;

Printf("%.0f %.0f %.0f %.0f %.0f", nx1, nx2, nx3, ny, nz);

/****************************************************************************
 *                               Points                                     *
 ****************************************************************************/
// back - down
Point(1)  = {0, uy_back, 0, cl};
Point(2)  = {ux_back - p_back, uy_back, 0,cl};
Point(3)  = {ux_back, uy_back, -h_back,cl};
Point(4)  = {ux_back + p_back, uy_back, 0,cl};
Point(5)  = {l, uy_back, 0,cl};

// front - down
Point(6)  = {0, uy_front, 0, cl};
Point(7)  = {ux_front - p_front, uy_front, 0,cl};
Point(8)  = {ux_front, uy_front, -h_front,cl};
Point(9)  = {ux_front + p_front, uy_front, 0,cl};
Point(10)  = {l, uy_front, 0,cl};

// top - back
Point(11)  = {0, uy_back, g, cl};
Point(12)  = {ux_back - p_back, uy_back, g,cl};
Point(13)  = {ux_back + p_back, uy_back, g,cl};
Point(14)  = {l, uy_back, g,cl};

// top - front
Point(15)  = {0, uy_front, g, cl};
Point(16)  = {ux_front - p_front, uy_front, g,cl};
Point(17)  = {ux_front + p_front, uy_front, g,cl};
Point(18)  = {l, uy_front, g,cl};

/****************************************************************************
 *                                 Lines                                    *
 ****************************************************************************/

// Lower Surface
Line(1) = {10, 5};
Line(2) = {5, 4};
Circle(3) = {4, 3, 2};
Line(4) = {2, 1};
Line(5) = {1, 6};
Line(6) = {6, 7};
Circle(7) = {7, 8, 9};
Line(8) = {9, 10};
Line(9) = {9, 4};
Line(10) = {2, 7};

// Upper surface
Line(11) = {14, 18};
Line(12) = {18, 17};
Line(13) = {17, 16};
Line(14) = {16, 15};
Line(15) = {15, 11};
Line(16) = {11, 12};
Line(17) = {12, 13};
Line(18) = {13, 14};
Line(27) = {13, 17};
Line(28) = {16, 12};

// Vertical Lines
Line(19) = {10, 18};
Line(20) = {9, 17};
Line(21) = {7, 16};
Line(22) = {6, 15};
Line(23) = {1, 11};
Line(24) = {2, 12};
Line(25) = {4, 13};
Line(26) = {5, 14};

// Transfinite 

// vertical stuff
Transfinite Line {19, 20, 21, 22, 23, 24, 25, 26} = nz Using Progression progz;

// in y direction
Transfinite Line {5, 10, -9, -1} = ny;
Transfinite Line {-15, -28, 27, 11} = ny;

// in x direction, from left to right
Transfinite Line {-4, 6, 16, -14} = nx1 Using Progression progx1;
Transfinite Line {-3, 7, 17, -13} = nx2 Using Progression progx2;
Transfinite Line {-2, 8, 18, -12} = nx3 Using Progression progx3;


/****************************************************************************
 *                                 Surfaces                                 *
 ****************************************************************************/

//back
Line Loop(29) = {-16, -23, -4, 24};
Plane Surface(30) = {29};
Line Loop(31) = {-17, -24, -3, 25};
Plane Surface(32) = {31};
Line Loop(33) = {26, -18, -25, -2};
Plane Surface(34) = {33};

//front
Line Loop(35) = {-12, -19, -8, 20};
Plane Surface(36) = {35};
Line Loop(37) = {-13, -20, -7, 21};
Plane Surface(38) = {37};
Line Loop(39) = {22, -14, -21, -6};
Plane Surface(40) = {39};
 
//bottom
Line Loop(41) = {4, 5, 6, -10};
Plane Surface(42) = {41};
Line Loop(43) = {10, 7, 9, 3};
Ruled Surface(44) = {43};
Line Loop(45) = {8, 1, 2, -9};
Plane Surface(46) = {45};
 
//top
Line Loop(47) = {14, 15, 16, -28};
Plane Surface(48) = {47};
Line Loop(49) = {13, 28, 17, 27};
Ruled Surface(50) = {49};
Line Loop(51) = {12, -27, 18, 11};
Plane Surface(52) = {51};

// middle separators
Line Loop(53) = {5, 22, 15, -23};
Plane Surface(54) = {53};
Line Loop(55) = {10, 21, 28, -24};
Ruled Surface(56) = {55};
Line Loop(57) = {-9, 20, -27, -25};
Plane Surface(58) = {57};
Line Loop(60) = {19, -11, -26, -1};
Plane Surface(61) = {60};

// Recombine and Transfinite
Transfinite Surface{30, 32, 34};
Transfinite Surface{36, 38, 40};
Transfinite Surface{42, 44, 46};
Transfinite Surface{48, 50, 52};
Transfinite Surface{54, 56, 58, 61};

Recombine Surface{30, 32, 34};
Recombine Surface{36, 38, 40};
Recombine Surface{42, 44, 46};
Recombine Surface{48, 50, 52};
Recombine Surface{54, 56, 58, 61};


/****************************************************************************
 *                                 Volumes                                  *
 ****************************************************************************/
// left
Surface Loop(62) = {42, 30, 48, 40, -54, 56};
Volume(63) = {62};

// mid
Surface Loop(64) = {44, 32, 50, 38, -56, 58};
Volume(65) = {64};

// right
Surface Loop(66) = {46, 34, 52, 36, -58, 61}; 
Volume(67) = {66};

// transfinite stuff
Transfinite Volume{63, 65, 67};

/****************************************************************************
 *                             Physical Lables                              *
 ****************************************************************************/
Physical Volume(100) = {63, 65, 67};

// perpendicular to x
Physical Surface(1) = {61, 54};

// perpendicular to y
Physical Surface(2) = {30, 32, 34, 40, 38, 36};

// perpendicular to z
Physical Surface(3) = {42, 46, 48, 50, 52};

// must be curved
Physical Surface(4) = {44};

