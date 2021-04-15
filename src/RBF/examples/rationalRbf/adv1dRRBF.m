% rational RBF solution of the 1d advection equation
%   u_t + u_x = 0  on the interval [-1,1] with a Dirichlet BC u(-1) = 1

function adv1dRRBF()
  phi = iqx();            % inverse quadratic RBF
  s = 6.0;                % shape parameter
  N = 100;                % number of centers
  k = 0.005;              % time step size
  finalTime = 0.75;
  x = rbfCenters.xKT(N,-1,1,0.999);   % boundary clustered centers
  
  r = rbfx.distanceMatrix1d(x);   
  B = phi.rbf(r,s);     % system matrix
  H = phi.D1(r,s,r);    % 1st derivative evaluation matrix
  [L, Bi] = rbfRational.rbfRationalCoreMultipleSetUp(B);
                                 
  v = exactSol(x,0);                                                                             
  exact = exactSol(x,0);
 
   hh1 = animatedline(x,v,'Color','b','LineWidth',2);
   hh2 = animatedline(x,exact,'Color','r','LineWidth',2,'LineStyle','--');
   axis([-1 1 -1.1 1.1])
 
   n=0;
   while n*k < finalTime   
      v = rbfMisc.rk4(v,n*k,k,@F);
      n = n+1;
      v(1) = 1.0;    % Dirichlet boundary condition

     exact = exactSol(x,n*k);
     clearpoints(hh1);   addpoints(hh1,x,v);
     clearpoints(hh2);   addpoints(hh2,x,exact);
     drawnow();
   end  % while

   exact = exactSol(x,finalTime);
   er = norm(v-exact,inf)
 

% ------------------------------------------------------------------------
   
 function ex = exactSol(x,t)
       x2 = x - t;
       ex = 1.*( x2<=-0.75 )  - 1.*( x2>-0.75 );
 end 
   
   
  function fp = F(u,t)
         u(1) = 1.0;
         ra1 = rbfRational.diff1dMultiple(u,B,Bi,L,H);
         fp =  -(ra1);
  end

% -----------------------------------------------------------------------
 
end