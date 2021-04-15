% diffRRBF.m - 1d interpolation and 1st and 2nd derivatives

clear, home, close all

phi = iqx();   % inverse quadratic RBF
s = 4.0;       % shape parameter
N = 100;       % number of centers
t = 0.53;      % function parameter
nu = 0.0025;   % viscosity coefficient

M = 201;       % number of evaluation points
xc = linspace(-1,1,N)';    % centers
x  = linspace(-1,1,M)';    % evaluation points

% ---- function and exact derivatives -------------------------------------

F = @(x,t,nu) (0.1*exp(-0.05*(x+1.0 - 0.5 + 4.95*t)/nu) + ...
    0.5*exp(-0.25*(x+1.0 - 0.5 + 0.75*t)/nu) + ...
    exp(-0.5*(x+1.0-0.375)/nu))./(exp(-0.05*(x+1.0 - 0.5 + 4.95*t)/nu) + ...
    exp(-0.25*(x+1.0 - 0.5 + 0.75*t)/nu) + exp(-0.5*(x+1.0-0.375)/nu));

FP = @(x,t,n) (-1/200).*exp(1).^((1/400).*n.^(-1).*(75+24.*t+100.*x)).*(25.*exp( ...
  1).^((99/400).*n.^(-1).*t)+81.*exp(1).^((1/80).*n.^(-1).*(8+15.*t+ ...
  16.*x))+16.*exp(1).^((1/80).*n.^(-1).*(23+36.*x))).*(exp(1).^(( ...
  99/400).*n.^(-1).*t)+exp(1).^((1/80).*n.^(-1).*(23+36.*x))+exp(1) ...
  .^((1/400).*n.^(-1).*(75+24.*t+100.*x))).^(-2).*n.^(-1);

FPP = @(x,t,n) (1/4000).*exp(1).^((1/400).*n.^(-1).*(75+24.*t+100.*x)).*(exp(1) ...
 .^((99/400).*n.^(-1).*t)+exp(1).^((1/80).*n.^(-1).*(23+36.*x))+ ...
 exp(1).^((1/400).*n.^(-1).*(75+24.*t+100.*x))).^(-3).*((-125).* ...
 exp(1).^((99/200).*n.^(-1).*t)+64.*exp(1).^((1/40).*n.^(-1).*(23+ ...
  36.*x))+(-729).*exp(1).^((1/200).*n.^(-1).*(20+87.*t+40.*x))+729.* ...
  exp(1).^((1/80).*n.^(-1).*(31+15.*t+52.*x))+125.*exp(1).^((1/400) ...
  .*n.^(-1).*(75+123.*t+100.*x))+(-64).*exp(1).^((1/200).*n.^(-1).*( ...
  95+12.*t+140.*x))+182.*exp(1).^((1/400).*n.^(-1).*(115+99.*t+180.* ...
  x))).*n.^(-2);

f = F(xc,t,nu);
fe = F(x,t,nu);
fe1 = FP(xc,t,nu);
fe2 = FPP(xc,t,nu);
m1 = max(abs(fe1));
m2 = max(abs(fe2));

% ------------------ standard RBF setup -----------------------------------

r = rbfx.distanceMatrix1d(xc);
re = rbfx.distanceMatrix1d(xc,x);
B = phi.rbf(r,s);
kappaB = cond(B)
H = phi.rbf(re,s);
H1 = phi.D1(r,s,r);
H2 = phi.D2(r,s,r);
a = rbfx.solve(B,f,0,false);

% ------------------ interpolation ----------------------------------------

sa = H*a;
ra = rbfRational.interpolate1d(phi,s,f,xc,x);
erInterpRat = norm( ra - fe, inf)
erInterpStd = norm( sa - fe, inf)
plot(x,ra,'b',x,fe,'r--')


% ---------------------- 1st derivative -----------------------------------

sa1 = H1*a;
[ra1, ra2] = rbfRational.diff1dOne(phi,s,f,xc,true);
erD1Rat = norm( ra1 - fe1, inf)/m1   % relative error
erD1Std = norm( sa1 - fe1, inf)/m1
figure()
plot(xc,ra1,'b',xc,fe1,'r--')

% ----------------------- 2nd derivative ----------------------------------

sa2 = H2*a;
erD2Rat = norm( ra2 - fe2, inf)/m2
erD2Std = norm( sa2 - fe2, inf)/m2
figure()
plot(xc,ra2,'b',xc,fe2,'r--')

% -------------------------------------------------------------------------

% relative pointwise errors
% semilogy(x, abs(ra - fe), 'g',x, abs(sa - fe), 'b--')
% semilogy(xc, abs(ra1 - fe1)/m1, 'g',xc, abs(sa1 - fe1)/m1, 'b--')
% semilogy(xc, abs(ra2 - fe2)/m2, 'g',xc, abs(sa2 - fe2)/m2, 'b--')