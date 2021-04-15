
% centroExtensions.m
%
% Depending on how the centers were extended to be symmetric, RBF
% differentiation matrices will have a (skew) centrosymmetric structure.
% The following reference can be consulted for details:
% "Radial Basis Function Methods - the case of symmetric domains."  
% Under review, Numerical Methods for Partial Differential Equations, 2016.

clear, home, close all

phi = iqx();

s = 4.5;        % shape parameter
N = 4000;      % number of centers on the covering square
mu = 5e-15;    % MDI regularization parameter
safe = false;  % use Cholesky factorization

%  symType     0 - y axis (-x, y)  -> odd order operators wrt x are skew-centro
%                                     odd oder operators wrt y and all even order operators are centro
%                                     mixed odd operators such as the divergence have neither type of symmetry
%                                        as the sum of a skew-centro and centro matrix is in general not skew-centro or centro
%
%              1 - x axis ( x,-y)  -> odd order operators wrt y are skew-centro
%                                     odd oder operators wrt x and all even order operators are centro
%                                     mixed odd operators such as the divergence have neither type of symmetry
%
%              2 - origin (-x,-y)  -> odd order operators skew-centro
%                                     even order operators centro

symType = 2    % set to 0, 1, or 2

[x, y] = rbfCenters.circleCenters(N,true,0,1,false);  
[xc,yc] = rbfCentro.centroCenters(x,y,symType,false);             
[r, rx, ry] = rbfx.distanceMatrix2d(xc,yc);

B = phi.rbf(r,s);
[kappaB, kappaL, kappaM] = rbfCentro.centroConditionNumber(B,mu);

H1x = phi.D1(r,s,rx);
H1y = phi.D1(r,s,ry);
H2 = phi.D2(r,s,rx);
G = phi.G(r, s, rx, ry);
L = phi.L(r, s);

Fn = F2b;
f = Fn.F(xc,yc);
a = rbfCentro.solveCentro(B,f,mu,safe);

disp(' ')
fprintf('system matrix condition number: %4.2e\n',kappaB);
fprintf('B:  '); rbfCentro.hasSymmetry(B);
fprintf('D1x: '); rbfCentro.hasSymmetry(H1x);
fprintf('D1y: '); rbfCentro.hasSymmetry(H1y);
fprintf('D2: '); rbfCentro.hasSymmetry(H2);
fprintf('divergence: '); rbfCentro.hasSymmetry(G);
fprintf('Laplacian: '); rbfCentro.hasSymmetry(L);
disp(' ')

d1xSymType = 1;
d1ySymType = 1;
d2SymType = 2;
gSymType = 1;
LSymType = 2;

if symType == 0
    d1ySymType = 2;
elseif symType == 1
    d1xSymType = 2;
end

[L1x,M1x] = rbfCentro.centroDecomposeMatrix(H1x,d1xSymType);          
fx = rbfCentro.centroMult(a,L1x,M1x,d1xSymType);
fprintf('d1x error = %4.2e\n',norm(fx - Fn.x1(xc,yc)));

[L1y,M1y] = rbfCentro.centroDecomposeMatrix(H1y,d1ySymType);          
fy = rbfCentro.centroMult(a,L1y,M1y,d1ySymType);
fprintf('d1y error = %4.2e\n',norm(fy - Fn.y1(xc,yc)));

[L2,M2] = rbfCentro.centroDecomposeMatrix(H2,d2SymType);
fxx = rbfCentro.centroMult(a,L2,M2,d2SymType);
fprintf('d2 error = %4.2e\n',norm(fxx - Fn.x2(xc,yc)));

% The divergence is only skew-centro with an origin extension.
% However, the centro and skew part can be calculated separately
if symType == 2
    [LG,MG] = rbfCentro.centroDecomposeMatrix(G,gSymType);
    fg = rbfCentro.centroMult(a,LG,MG,gSymType);
    fprintf('divergence error = %4.2e\n',norm(fg - Fn.G(xc,yc)));
elseif symType == 0 || symType == 1
    fprintf('divergence error = %4.2e\n',norm((fx + fy) - Fn.G(xc,yc)));
end

[LL,ML] = rbfCentro.centroDecomposeMatrix(L,LSymType);
fL = rbfCentro.centroMult(a,LL,ML,LSymType);
fprintf('Laplacian error = %4.2e\n',norm(fL - Fn.L(xc,yc)));


