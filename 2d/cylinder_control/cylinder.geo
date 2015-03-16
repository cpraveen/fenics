// Authors: Anant Diwakar, Praveen. C

// Uncomment following two lines to generate 8 node quadrilaterals
//Mesh.ElementOrder = 2;
//Mesh.SecondOrderIncomplete = 1;

lc = 0.02;
lo = 0.05;  // outer boundary
lw = 0.05;  // wake outlet

d = 0.1;
r = 0.5*d;
R = 5*r;

assym = 0.0*d;

h = 15*d;
l1 = 15*d;
l2 = 35*d;

n = 2;
nl1 = 4;  // along x before cylinder
nl2 = 150; // along wake
nr  = 40;  // along radius
nw  = 10; // normal to side wall

nc = 15; // around cylinder

//nz = 30;
nz = 1;

Point(1) = {0,0,0,lc};
Point(2) = {0,r,0,lc};
Point(3) = {r,0,0,lc};
Point(4) = {-r,0,0,lc};
Point(5) = {-l1,0,0,0.2*lo};
Point(6) = {l2,0,0,lw};
Point(7) = {-l1,h,0,lo};
Point(8) = {l2,h,0,lo};
Point(9) = {R,0.0,0,lc};
Point(10) = {-R,0.0,0,lc};
Point(11) = {0.0,R,0,lc};
Point(12) = {0.0,h,0,lo};
Point(13) = {r*Sin(Pi/4),r*Sin(Pi/4),0,lc};
Point(14) = {-r*Sin(Pi/4),r*Sin(Pi/4),0,lc};
Point(15) = {R*Sin(Pi/4),R*Sin(Pi/4),0,lc};
Point(16) = {-R*Sin(Pi/4),R*Sin(Pi/4),0,lc};
Point(17) = {l2,R*Sin(Pi/4),0,lw};
Point(18) = {-l1,R*Sin(Pi/4),0,lo};
Point(19) = {R*Sin(Pi/4),h,0,lo};
Point(20) = {-R*Sin(Pi/4),h,0,lo};

Delete {
  Point{19};
}

Delete {
  Point{20};
}

Circle(1) = {15,1,11};
Circle(2) = {13,1,2};
Circle(3) = {2,1,14};
Circle(4) = {11,1,16};
Circle(5) = {16,1,10};
Circle(6) = {14,1,4};
Circle(7) = {13,1,3};
Circle(8) = {15,1,9};
Line(9) = {9,3};
Line(10) = {15,13};
Line(11) = {11,2};
Line(12) = {16,14};
Line(13) = {10,4};
Line(14) = {12,11};


Line(15) = {15,17};
Line(16) = {17,8};
Line(17) = {8,12};
Line(18) = {12,7};
Line(19) = {7,5};
Line(20) = {18,16};
//Line(21) = {18,5};
Line(22) = {5,10};
Line(23) = {6,17};
Line(24) = {9,6};


Delete {
  Line{20};
}

Delete {
  Point{18};
}


Line Loop(25) = {14,-1,15,16,17};
Plane Surface(26) = {25};

Line Loop(27) = {-14,18,19,22,-5,-4};
Plane Surface(28) = {27};

Line Loop(29) = {7,-9,-8,10};
Plane Surface(30) = {29};

Line Loop(31) = {-2,-10,1,11};
Plane Surface(32) = {31};

Line Loop(33) = {-3,-11,4,12};
Plane Surface(34) = {33};

Line Loop(35) = {-6,-12,5,13};
Plane Surface(36) = {35};

Line Loop(37) = {8,24,23,-15};
Plane Surface(38) = {37};

Symmetry { 0.0,1.0,0.0,0.0 } {
  Duplicata { Surface{38,26,28,36,34,32,30}; }
}

Translate {0.0,assym,0.0} {
  Point{7};
}

Translate {0.0,assym,0.0} {
  Point{12};
}

Translate {0.0,assym,0.0} {
  Point{8};
}


Transfinite Line {15} = nl2+1 Using Progression 1.01;
Transfinite Line {24} = nl2+1 Using Progression 1.01;
Transfinite Line {-43} = nl2+1 Using Progression 1.01;

// radially away from cylinder
Transfinite Line{-9,-10,-11,-12,-13,59,64,69} = nr+1 Using Progression 1.05;

Transfinite Line {8,5,6,7,4,3,2,1,23} = nc+1 Using Progression 1.0;

Transfinite Line {55,40,73,58,56,63,68,46,42} = nc+1 Using Progression 1.0;

Transfinite Surface {38} = {9,6,17,15};
Transfinite Surface {30} = {3,9,15,13};
Transfinite Surface {32} = {2,13,15,11};
Transfinite Surface {34} = {14,2,11,16};
Transfinite Surface {36} = {4,14,16,10};

//Recombine Surface {38,30,32,34,36,28,26};

Transfinite Surface {39} = {30,20,9,6};
Transfinite Surface {72} = {20,122,3,9};
Transfinite Surface {67} = {38,104,122,20};
Transfinite Surface {62} = {76,86,104,38};
Transfinite Surface {57} = {76,10,4,86};

//Recombine Surface {50,57,62,67,72,39,44};

Delete {
  Point{1};
}

Physical Surface(100) = {38,30,32,34,36,28,26,-50,-57,-62,-67,-72,-39,-44};

Physical Line(1) = {19,53};               // Inlet
Physical Line(2) = {7,2,3,6,58,63,68,73}; // cylinder
Physical Line(3) = {16,23,42,48};         // outlet
Physical Line(4) = {18,17,52,49};         // top and bottom wall

//Geometry.Normals = 100;

