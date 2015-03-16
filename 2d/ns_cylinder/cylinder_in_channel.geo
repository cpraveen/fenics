D    = 0.1;   // diameter
r    = D/2.0; // radius

h1   =  0.003; // on cylinder
h2   =  0.05; // on boundary
h3   =  0.01; // at inlet
xmin = -1.5;  // left end
xmax =  2.2;  // right end
H    =  0.2;  // half channel height
xc   =  0.25; // center of cylinder

Point(1) = {xmin, 0, 0, h3};
Point(2) = {-r+xc,0, 0, h1};
Point(3) = { xc,  r, 0, h1};
Point(4) = { r+xc,0, 0, h1};
Point(5) = {xmax, 0, 0, h2};
Point(6) = {xmax, H, 0, h2};
Point(7) = {xmin, H, 0, h3};

Point(8) = {xc, 0, 0, h1};

Line(1)   = {1, 2};
Circle(2) = {2, 8, 3};
Circle(3) = {3, 8, 4};
Line(4)   = {4, 5};
Line(5)   = {5, 6};
Line(6)   = {6, 7};
Line(7)   = {7, 1};

Line Loop(1) = {1,2,3,4,5,6,7};
Plane Surface(1) = {1};
Symmetry{0,1,0,0}{ Duplicata{Surface{1};} }

Physical Surface(1000000) = {1, 8};
Physical Line(1000001) = {7, 15}; // inlet
Physical Line(1000002) = {6, 14}; // side walls
Physical Line(1000003) = {5, 13}; // outlet
Physical Line(1000005) = {2, 3, 10, 11}; // cylinder
