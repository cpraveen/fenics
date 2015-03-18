Mesh.CharacteristicLengthFromPoints = 0;

D    = 0.1;   // diameter
r    = D/2.0; // radius

theta = 60;
phi = 20;

theta = (Pi/180.0)*theta;
phi =  (Pi/180.0)*phi;

h1   =  0.003; // on cylinder
h2   =  0.06; // on boundary
h3   =  0.01; // at inlet
xmin = -1.5;  // left end
xmax =  2.2;  // right end
H    =  0.2;  // half channel height
xc   =  0.25; // center of cylinder

Point(1) = {xmin, 0, 0, h3};
Point(2) = {-r+xc,0, 0, h1};
Point(3) = { xc + r*Cos(theta + 0.5*phi),  r*Sin(theta + 0.5*phi), 0, h1};
Point(4) = { xc + r*Cos(theta - 0.5*phi),  r*Sin(theta - 0.5*phi), 0, h1};
Point(5) = { r+xc,0, 0, h1};
Point(6) = {xmax, 0, 0, h2};
Point(7) = {xmax, H, 0, h2};
Point(8) = {xmin, H, 0, h3};

Point(9) = {xc, 0, 0, h1};

Line(1)   = {1, 2};
Circle(2) = {2, 9, 3};
Circle(3) = {3, 9, 4};
Circle(4) = {4, 9, 5};
Line(5)   = {5, 6};
Line(6)   = {6, 7};
Line(7)   = {7, 8};
Line(8)   = {8, 1};

Line Loop(1) = {1,2,3,4,5,6,7,8};
Plane Surface(1) = {1};
Symmetry{0,1,0,0}{ Duplicata{Surface{1};} }

Field[1] = MathEval;
Field[1].F = "(x+1.5)*0.5/3.5 - (x-2)*0.005/3.5";

Field[2] = MathEval;
Field[2].F = "0.0025*Sqrt((x-0.25)^2+y*y)/0.05";

Field[3] = MathEval;
Field[3].F = "0.02";

Field[4] = Min;
Field[4].FieldsList = {1,2,3};

Background Field = 4;

// Don't extend the elements sizes from the boundary inside the domain
Mesh.CharacteristicLengthExtendFromBoundary = 0;

Physical Surface(1000000) = {1, -9};
Physical Line(1) = {8, 17};        // inlet
Physical Line(2) = {7, 16};        // side walls
Physical Line(3) = {6, 15};        // outlet
Physical Line(4) = {2, 4, 11, 13}; // cylinder
Physical Line(5) = {3};            // feedback boundary
Physical Line(6) = {12};           // feedback boundary
