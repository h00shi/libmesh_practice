cl__1 = 1e+22;
Point(1) = {0, 0, 0, 1e+22};
Point(2) = {0.25, 0, 0, 1e+22};
Point(3) = {1, 0, 0, 1e+22};
Point(4) = {0, 0.25, 0, 1e+22};
Point(5) = {0, 1, 0, 1e+22};
Circle(1) = {3, 1, 5};
Line(2) = {5, 4};
Circle(3) = {4, 1, 2};
Line(4) = {2, 3};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};
Physical Surface(501) = {6};
