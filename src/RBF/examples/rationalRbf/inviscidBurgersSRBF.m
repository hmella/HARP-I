% standard RBF approximation of insviscid Burgers' equation
%   u_t + 0.5(u^2)_x = 0   
% the solution forms a shock at time t = 1/pi and the standard RBF
%   method solution becomes unstable and blows up

function inviscidBurgersSRBF()
  phi = iqx();             % inverse quadratic RBF   
  s = 6.0;                 % shape parameter
  N = 100;                 % number of centers
  k = 0.005;               % time step size
  finalTime = 0.7;
  x = rbfCenters.xKT(N,-1,1,0.999);  % boudary clustered centers

   r = rbfx.distanceMatrix1d(x);
   B = phi.rbf(r,s);
   H1 = phi.D1(r,s,r);
   D1 = phi.dm(B,H1,0,true); 
                                  
   v = sin(pi*(x+1));
   hh1 = animatedline(x,v,'Color','b','LineWidth',2);
   axis([-1 1 -1.2 1.2])
   n=0;
   while n*k < finalTime   
      v = rbfMisc.rk4(v,n*k,k,@F);
      n = n+1;
      v(1) = 0; v(end) = 0;
      clearpoints(hh1);   addpoints(hh1,x,v);  drawnow()
   end  % while

% ------------------------------------------------------------------------

  function fp = F(u,t)
         u(1) = 0; u(end) = 0;
         fp =  -D1*(0.5*u.^2);   % conservative form
%         fp =  -u.*(D1*u);      % non-conservative form
  end 

% -----------------------------------------------------------------------
 
end