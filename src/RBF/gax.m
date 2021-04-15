

classdef gax < rbfx
   methods
       function obj = gax()  % constructor
          obj@rbfx();  % call constructor of the superclass 
       end
       
       function v = rbf(obj,r,s), v = exp( -(s.*r).^2 ); end
       
       function d = D1(obj,r,s,x), d = -2*x.*s.^2.*exp(-(s.*r).^2); end
        
       function d = D2(obj, r, s, x) 
           d = 2.0*s.^2.*(2*x.*x.*s.^2 - 1.0).*exp(-(s.*r).^2);
       end
       
       function d = D3(obj, r, s, x) 
           d = -4*exp(-(s.*r).^2).*x.*s.^4.*(2*s.^2.*x.^2 - 3);
       end
       
       function d = D4(obj, r, s, x) 
           d = 4*exp(-(s.*r).^2).*s.^4.*(3 - 12*s.^2.*x.^2 + 4*s.^4.*x.^4);
       end
       
       function d = G(obj, r, s, x, y)   % Gradient
           d = -2*exp(-(s.*r).^2).*s.^2.*(x + y);
       end
       
       function d = L(obj, r, s) 
           d = 4*exp(-(s.*r).^2).*s.^2.*(r.^2.*s.^2 - 1.0);
       end
       
       function d = B(obj, r, s, x, y) 
           d = 16*exp(-(s.*r).^2).*s.^4.*( 2 - 4*y.^2.*s.^2 + x.^4.*s.^4 + y.^4.*s.^4 + 2*x.^2.*s.^2.*( y.^2.*s.^2 - 2 ) );
       end
       
% D12       
% mixed partial derivative
%      D_{xyy}  d = D12( r, s, x, y ) 
% or   D_{yxx}  d = D12( r, s, y, x) 
% depending on the order of the x and y arguments
       
       function d = D12(obj, r, s, x, y) 
           d = exp(-(s.*r).^2).*x.*( 4*s.^4 - 8*y.^2.*s.^6 );
       end
       
       function d = D22(obj, r, s, x, y) 
           d = 4*exp(-(s.*r).^2).*s.^4.*( 2*x.^2.*s.^2 - 1 ).*( 2*y.^2.*s.^2 - 1 );
       end
      
    
  end % public methods

end  % class
       
