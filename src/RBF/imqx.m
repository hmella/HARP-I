
classdef imqx < rbfx
    methods
        function obj = imqx()  % constructor
          obj@rbfx();  % call constructor of the superclass 
        end
       
        function v = rbf(obj,r,s), v = 1./sqrt(1 + (s.*r).^2 ); end
       
        function d = D1(obj,r,s,x), d = -(x.*s.^2)./(1.0 + (s.*r).^2 ).^(1.5); end
       
        function d = D2(obj, r, s, x)
            d = s.^2.*(-r.^2.*s.^2 + 3*s.^2*x.^2 - 1)./(r.^2*s.^2 + 1).^(2.5);
        end
        
        function d = D3(obj, r, s, x) 
            d = 3*s^4*x.*( -5*s.^2.*x.^2 + 3*s.^2.*r.^2 + 3  )./(r.^2.*s.^2 + 1).^(3.5);
        end
        

        function d = D4(obj, r, s, x) 
            d = 3*s.^4.*(35*s.^4*x.^4 - 30*s.^2*x.^2.*(s.^2.*r.^2 + 1) + 3*(s.^2.*r.^2 + 1).^2)./(r.^2.*s.^2 + 1).^(4.5); 
 
        end
        

        function d = G(obj, r, s, x, y)   % Gradient
           d = -s.^2.*(x + y)./(r.^2.*s.^2 + 1).^(1.5);
        end
        
        function d = L(obj, r, s)         % Laplacian
           d = s.^2.*(r.^2.*s.^2 - 2)./(1 + (s.*r).^2 ).^(2.5);
        end
        
        function d = B(obj, r, s, x, y)   % Biharmonic operator   
           d = 3*s.^4.*(35*s.^4.*x.^4 + 70*s.^4.*x.^2.*y.^2 + 35*s.^4.*y.^4 - 30*s.^2.*x.^2.*(s.^2.*r.^2 + 1) - 30*s.^2.*y.^2.*(s.^2.*r.^2 + 1) - 10*s.^2.*r.^2.*(s.^2.*r.^2 + 1) + 8*(s.^2.*r.^2 + 1).^2)./(1 + (s.*r).^2 ).^(4.5);
        end

     
% D12       
% mixed partial derivative
%      D_{xyy}  d = D12( r, s, x, y ) 
% or   D_{yxx}  d = D12( r, s, y, x) 
% depending on the order of the x and y arguments
       
       function d = D12(obj, r, s, x, y) 
           d = 3*x.*s.^4.*(1 - 4*y.^2.*s.^2 + x.^2.*s.^2  )./(1 + (s.*r).^2 ).^(3.5);
       end
       
       function d = D22(obj, r, s, x, y) 
           d = -3*s.^4.*( -1 + 3*y.^2.*s.^2 + 4*x.^4.*s.^4 + 4*y.^4.*s.^4 + x.^2.*(3*s.^2 - 27.*y.^2*s.^4)  )./(1 + (s.*r).^2 ).^(4.5);
       end
      
   end % methods
end  % class