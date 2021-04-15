
classdef iqx < rbfx
    methods
        function obj = iqx()  % constructor
          obj@rbfx();  % call constructor of the superclass 
        end
       
        function v = rbf(obj,r,s), v = 1./(1 + (s.*r).^2 ); end
       
        function d = D1(obj,r,s,x), d = -(2*x.*s.^2)./(1.0 + (s.*r).^2 ).^2; end
       
        function d = D2(obj, r, s, x)
            d = 2*s.^2.*(-r.^2.*s.^2 + 4*s.^2*x.^2 - 1)./(r.^2*s.^2 + 1).^3;
        end
        
        function d = D3(obj, r, s, x) 
            d = (-48*x.^3.*s.^6)./(r.^2.*s.^2 + 1).^4 + (24*x.*s.^4)./(r.^2.*s.^2 + 1).^3;
        end
        
        function d = D4(obj, r, s, x) 
            d = (384*x.^4.*s.^8)./(r.^2.*s.^2 + 1).^5 - (288*x.^2.*s.^6)./(r.^2.*s.^2 + 1).^4 + (24*s.^4)./(r.^2.*s.^2 + 1).^3; 
        end
        
        function d = G(obj, r, s, x, y)   % Gradient
           d = -2*s.^2.*(x + y)./(r.^2.*s.^2 + 1).^2;
        end
        
        function d = L(obj, r, s)         % Laplacian
           d = 4*s.^2.*(r.^2.*s.^2 - 1)./(1 + (s.*r).^2 ).^3;
        end
        
         % x and y not used but required by the abstract function definition in the superclass
        function d = B(obj, r, s, x, y)   % Biharmonic operator   
           d = 64*( s.^4 - 4*r.^2.*s.^6 + r.^4.*s.^8  )./(1 + (s.*r).^2 ).^5;
        end

     
% D12       
% mixed partial derivative
%      D_{xyy}  d = D12( r, s, x, y ) 
% or   D_{yxx}  d = D12( r, s, y, x) 
% depending on the order of the x and y arguments
       
       function d = D12(obj, r, s, x, y) 
           d = 8*x.*s.^4.*(1 - 5*y.^2.*s.^2 + x.^2.*s.^2  )./(1 + (s.*r).^2 ).^4;
       end
       
       function d = D22(obj, r, s, x, y) 
           d = -8*s.^4.*(-1 + 4*r.^2.*s.^2 + 5*(x.^4 + y.^4).*s.^4 - 38*x.^2.*y.^2.*s.^4  )./(1 + (s.*r).^2 ).^5;
       end
      
   end % methods
end  % class