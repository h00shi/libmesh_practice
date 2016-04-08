/****************************************************************************
 *                                    Input                                 *
 ****************************************************************************/
v_max = Pi/6;
r0 = 1;
h0 = r0 * Cos(v_max);
l = 5 * r0 * Sin(v_max);
w = l/2;
g = l/3;
b= 1;
c= 4;

cl = 1e-1;
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
h_front= h0 - b * w / c;
p_front= Sqrt(r_front*r_front - h_front * h_front);

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
