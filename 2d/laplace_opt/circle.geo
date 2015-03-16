n1 = 70;
r  = 5;
h1 = 2*Pi*r/n1;

n2=35;
h2= 2*Pi/n2;

Point(1) = { 0,  0, 0};
Point(2) = { r,  0, 0};
Point(3) = { 0,  r, 0};
Point(4) = {-r,  0, 0};
Point(5) = { 0, -r, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Field[1] = Cylinder;
Field[1].VIn  = h2;
Field[1].VOut = h1;
Field[1].XCenter = 0;
Field[1].YCenter = 0;
Field[1].ZCenter = 0;
Field[1].ZAxis   = 1;
Field[1].Radius  = 1.5;

Field[2] = Cylinder;
Field[2].VIn  = h2;
Field[2].VOut = h1;
Field[2].XCenter = -3;
Field[2].YCenter = 0;
Field[2].ZCenter = 0;
Field[2].ZAxis   = 1;
Field[2].Radius  = 1.5;

Field[3] = Cylinder;
Field[3].VIn  = h2;
Field[3].VOut = h1;
Field[3].XCenter = 0;
Field[3].YCenter = -3;
Field[3].ZCenter = 0;
Field[3].ZAxis   = 1;
Field[3].Radius  = 1.5;

Field[4] = Min;
Field[4].FieldsList = {1,2,3};

Background Field = 4;
