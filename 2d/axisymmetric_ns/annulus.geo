r1 = 1.0; // radius of inner cylinder
r2 = 2.0; // radius of outer cylinder
h  = 2.0; // height of cylinder

n1     = 25; // points along radial
n2     = 25; // points along axial

p      = 0.1;// progression

// First cell height
dh = 0.5*(r2-r1)*p/((1+p)^(n1/2) - 1);
Printf("First cell height = %e\n", dh);

Point(1) = {r1, 0, 0};
Point(2) = {r2, 0, 0};
Point(3) = {r2, h, 0};
Point(4) = {r1, h, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,4};
Transfinite Line{1,3} = n1 Using Bump p;
Transfinite Line{2,4} = n2 Using Bump p;

Physical Line(100001) = {1};
Physical Line(100002) = {2};
Physical Line(100003) = {3};
Physical Line(100004) = {4};

Physical Surface(100005) = {1};
