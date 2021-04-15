% interpBenchExtended.m

warning off
tic

QUAD = true;   % true for quadruple precision, false for double

K = 1.5*sqrt(2);
[xc,yc,bpi] =  rbfCentersLib(4);   xc = xc/K;    yc = yc/K; 
[x,y,~] =  rbfCentersLib(401);    x = x/K;      y = y/K;  

 
if QUAD  % convert centers/execution points to extended, then all other calcs done in xprec
    mp.Digits(34); x = mp(x); y = mp(y);  xc = mp(xc);  yc = mp(yc);  
end
     
N = length(xc); M = length(x);
     
fn = F2d();          %  Franke function
f = fn.F(xc,yc);
fe = fn.F(x,y);

phi = iqx();

[r, rx, ry] = rbfx.distanceMatrix2d(xc,yc);
[re, rx, ry] = rbfx.distanceMatrix2d(xc,yc,x,y);
    
S = 7.2:-0.5:0.2;     Sn = length(S);
for k = 1:Sn
    s = S(k);
    B = phi.rbf(r,s);  
    a = rbfx.solve(B,f);   
    H = phi.rbf(re,s);
    fa = H*a;
end

if QUAD, tx = toc, else, td = toc, end


% after running the script in both precisions the ratio of times is:
% tx/td     % ratio of execution times
%   51.9    % typical ratio - quadruple takes about 52 times longer

warning on
