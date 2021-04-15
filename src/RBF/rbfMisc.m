

classdef rbfMisc


    
methods(Static)
        
        
% ------- 4th order Runge-Kutta -------------------------------------------

% rk4  Fourth-order Runge-Kutta ODE method. Advances an ODE system from
%      time t^n to time t^{n+1}
%
% input
%
%   V  solution at time t^n
%   t  time t^n
%   k  time step size delta t
%   F  the F of the ODE system v_t = F[v,t]
%
% output
%
%   v  solution at time t^{n+1}


 function v = rk4(V,t,k,F)
    s1 = feval(F,V,t);                                 
    s2 = feval(F,V + k*s1/2,t+k/2);                     
    s3 = feval(F,V + k*s2/2,t+k/2);                        
    s4 = feval(F,V + k*s3,t+k);                       
             
    v = V + k*(s1 + 2*s2 + 2*s3 + s4)/6;  
 end
 


% -------------------------------------------------------------------------


end % methods


    
% --------------------------------------------------------------------------- 
    
   
end % classdef