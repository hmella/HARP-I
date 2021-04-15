% rational RBF approximation of insviscid Burgers' equation
%   u_t + 0.5(u^2)_x = 0   
% the solution forms a shock at time t = 1/pi

function inviscidBurgersRRBF()
  phi = iqx();             % inverse quadratic RBF   
  s = 6.0;                 % shape parameter
  N = 100;                 % number of centers
  k = 0.005;               % time step size
  finalTime = 0.75;
  x = rbfCenters.xKT(N,-1,1,0.999);  % boudary clustered centers

  r = rbfx.distanceMatrix1d(x);   
  B = phi.rbf(r,s);     % system matrix
  H = phi.D1(r,s,r);    % 1st derivative evaluation matrix
  [L, Bi] = rbfRational.rbfRationalCoreMultipleSetUp(B);

  v = sin(pi*(x+1));
  hh1 = animatedline(x,v,'Color','b','LineWidth',2);

   axis([-1 1 -1.2 1.2])
   n=0;
   while n*k < finalTime   
      v = rbfMisc.rk4(v,n*k,k,@F);
      n = n+1;
      v(1) = 0; v(end) = 0;   % zero Dirichlet BCs
      clearpoints(hh1);   addpoints(hh1,x,v);  drawnow()
   end  % while
   
   exact = zeros(N,1);
   for i=1:N
     exact(i) = exactSolution(x(i),finalTime);
   end
   
   er = norm( v - exact, inf)
   plot(x,v,'b',x,exact,'r--');
   
   figure()
   semilogy( x(2:end-1), abs(v(2:end-1) - exact(2:end-1)), 'b')
  
% -----------------------------------------------------------------------
   
  function fp = F(u,t)
     u(1) = 0; u(end) = 0;
         
     % non-conservative form: works, but slight overshoot at the shock
%       ra1 = rbfRational.diff1dMultiple(u,B,Bi,L,H);
%       fp =  -u.*(ra1);
         
      % conservative form works better
        ra1 = rbfRational.diff1dMultiple(0.5*u.^2,B,Bi,L,H);
        fp =  -ra1;  
  end

    function s = exactSolution(X,T)
        shockTime = 1/pi;
        f = @(x) X - x + T*sin(pi*x);

        % after t=1/pi, for some X there will be more
        % than one zero, find the correct one
           if X>=0                
             x0 = fzero(f,1);
           else
             x0 = fzero(f,-1);
           end
          s = -sin(pi*x0);
    end

% -----------------------------------------------------------------------
 
end