model = setModel('Schnakenberg');
 
fem = Bilinear2D;
geometry = RectangleR(0,100,0,50,model.h);

p = zeros(geometry.nPoints,1); 
x = geometry.x; 
y = geometry.y;
p(( (x-50).^2+(y-25).^2) <= 4) = -3;

% Initial cond  
u0 = model.u0*ones(geometry.nPoints,1) + p;
v0 = model.v0*ones(geometry.nPoints,1);

 
FEM.solveImplicit(fem,geometry,model,u0,v0,2000,model.dt)