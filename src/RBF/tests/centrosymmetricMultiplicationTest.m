% centrosymmetricMultiplicationTest

clear, home, format compact

phi = gax();
s = 15;  N = 60;
x = -cos((0:N-1)*pi/(N-1))';
f = x.^2;  f1 = 2*x; f2 = 2.0;
mu = 5e-15;
safe = false;   %

% ----------------- order 0 ---------------------------


r =  phi.distanceMatrix1d(x);
B = phi.rbf(r,s);
fprintf('kappa(B) = %4.1e\n',cond(B));
rc = phi.distanceMatrix1d(x(1:N/2),x); 
Bc = phi.rbf(rc,s);

% ------------- order 1 -------------------------------

rho = 1;
Hc = phi.D1(rc,s,rc);
D1c = rbfCentro.centroDM(Bc,Hc,N,rho,mu,safe);
[L,M] = rbfCentro.centroDecomposeMatrix(D1c,rho);
ac = rbfCentro.centroMult(f,L,M,rho);
fprintf('f1 centro error = %4.2e\n',norm(f1 - ac, 2));

H = phi.D1(r,s,r);
D1 = phi.dm(B,H,mu,safe);
a = D1*f;
fprintf('f1 std error = %4.2e\n',norm(f1 - a, 2));

% ------------ order 2 -------------------------------

rho = 2;
Hc = phi.D2(rc,s,rc);
D2c = rbfCentro.centroDM(Bc,Hc,N,rho,mu,safe);
[L,M] = rbfCentro.centroDecomposeMatrix(D2c,rho);
ac = rbfCentro.centroMult(f,L,M,rho);
fprintf('f2 centro error = %4.2e\n',norm(f2 - ac, 2));

H = phi.D2(r,s,r);
D2 = phi.dm(B,H,mu,safe);
a = D2*f;
fprintf('f2 std error = %4.2e\n',norm(f2 - a, 2));



