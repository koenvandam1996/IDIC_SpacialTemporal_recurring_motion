function phi = dofbuild_bspline(na,nb,p)
% this functions builds a phi matrix which specifies the "control net" or
% "mesh" for the B-spline functions where;
%    na = number of knots along the x direction
%    nb = number of knots along the y direction
%    p  = the polynomial order of the B-spline functions
% The resulting phi matrix can be directly supplied to globalDIC2D.m

% number of knots
na_knots = na + 2*p;
nb_knots = nb + 2*p;

na_fun   = na_knots-(p+1);
nb_fun   = nb_knots-(p+1);

x = 1:na_fun;
y = 1:nb_fun;

[X Y] = meshgrid(x,y);


P = [X(:) Y(:)];

N = size(P,1);

for k = 1:3:N*3
    kk = ceil(k/3);
    phi(k  ,:) = [P(kk,:) 1 na nb p];
    phi(k+1,:) = [P(kk,:) 2 na nb p];
    phi(k+2,:) = [P(kk,:) 3 na nb p];
end
