% stepFunctionRRBF.m - interpolation of a 1d discontinuous function

close all
phi = iqx();                          % inverse quadratic RBF

fun = @(x) -1.0*(x<0) + 1.0*(x>=0);   % discontinuous function

M = 201;            % number of evaluation points
N = 80; s = 3.5;    % number of centers and corresponding shape
%N = 10; s = 0.45;

xc = linspace(-1,1,N)';     % centers
x  = linspace(-1,1,M)';     % evaluation points
f = fun(xc);                
fe = fun(x);

% ------- rational RBF method ---------------------------------------------

ra = rbfRational.interpolate1d(phi,s,f,xc,x);

% -------- standardard RBF method -----------------------------------------

r = rbfx.distanceMatrix1d(xc);
re = rbfx.distanceMatrix1d(xc,x);

B = phi.rbf(r,s);
H = phi.rbf(re,s);
a = rbfx.solve(B,f,0,false);
fs = H*a;

% -------------------------------------------------------------------------

plot(x,ra,'b',x,fs,'g',x,fe,'r--') %,xc,f,'k*')
legend('rational','standard','exact','Location','northwest')
xlabel('x'), ylabel('f(x)')

figure()     % pointwise error
semilogy(x,abs( ra - fe),'b', x, abs( fs - fe ), 'g')
legend('rational','standard')
xlabel('x'), ylabel('|error|')