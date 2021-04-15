% interp2dRRBF.m - interpolation of a 2d function with step gradients

phi = iqx();           % inverse quadratic RBF
s = 8;                 % shape parameter
N = 5000;              % number of centers            
RAT = true;            %  true -> rational,  false -> standard

 % centers, quasi-random Hamersley points on the unit square
[xc, yc] = rbfCenters.squareCenters(N,0,1,false,0,false);  
x = linspace(0,1,60);
[x,y] = meshgrid(x,x);    % 3600 uniformly space evaluation points 
x = x(:); y = y(:);

% --------------- the function being interpolated -------------------------

% f(x,y) = atan( 125[ sqrt( [x - x0]^2 + [y - y0]^2 ) - 0.92 ] )

x0 = 1.5; y0 = 0.25;    % function parameters
r = rbfx.distanceMatrix2d(xc,yc,x0,y0);
re = rbfx.distanceMatrix2d(x,y,x0,y0);
fun = @(rx) atan(125*(rx - 0.92));
f = fun(r);
fe = fun(re);

% -------------- rational RBF method --------------------------------------

if RAT
   tic
   ra = rbfRational.interpolate2d(phi,s,f,xc,x,yc,y);
   toc

% --------------- standard RBF method -------------------------------------

else
    tic
    [r, rx, ry] = phi.distanceMatrix2d(xc,yc);
    B = phi.rbf(r,s);

    mu = 0;
    safe = true;              % use mldivide rather than Cholesky directly
    a = rbfx.solve(B,f',mu,safe);

    [re, rx, ry] = phi.distanceMatrix2d(xc,yc,x,y);
    H = phi.rbf(re,s);

    ra = H*a;
    toc

end

% -------------------------------------------------------------------------

er = norm( ra - fe', inf)

tri = delaunay(x,y);
trisurf(tri,x,y,ra);
view(150,19);
colormap parula


