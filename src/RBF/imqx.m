% Licensing:
%   Copyright (c) 2025 Hernan Mella
%
%   MRBFT is under the GNU General Public License ("GPL").
%   
%   GNU General Public License ("GPL") copyright permissions statement:
%   **************************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%   This file was based on rbfx class by Scott Sarra.
%   Modifications were made by Hernan Mella on 01/2025.
%
% IMQX - Implements the Inverse Multiquadric Radial Basis Function (RBF)
%
% Description:
%   This class defines the Inverse Multiquadric (IMQ) RBF and its associated methods 
%   for interpolation and differentiation. IMQ RBFs are commonly used in scattered 
%   data approximation and meshless methods due to their smoothness and flexibility.
%   The IMQ RBF is represented as: 
%       φ(r) = 1 / sqrt(1 + (a*r)^2),
%   where `a` is the shape parameter controlling the support size of the RBF.
%
% Syntax:
%   obj = imqx();          % Create an instance of the IMQ RBF class
%   v = obj.rbf(r, a);     % Evaluate the IMQ RBF
%
% Inputs:
%   r - Distance matrix or distances between points (scalar or array).
%   a - Shape parameter controlling the support size of the RBF.
%
% Outputs:
%   v - RBF values evaluated at the given distances.
%
% Example:
%   % Create an instance of the IMQ RBF class
%   phi = imqx();
%
%   % Define distances and shape parameter
%   r = linspace(0, 2, 100);  % Distance matrix
%   a = 1.0;                  % Shape parameter
%
%   % Evaluate the IMQ RBF
%   v = phi.rbf(r, a);
%
%   % Plot the RBF
%   plot(r, v);
%   xlabel('Distance (r)');
%   ylabel('RBF Value');
%   title('Inverse Multiquadric Radial Basis Function');
%
% Methods:
%   rbf(r, a)       - Computes the IMQ RBF values for a given distance matrix `r` 
%                     and shape parameter `a`.
%   D1(r, a, x)     - Computes the first derivative of the IMQ RBF with respect to `x`.
%   D2(r, a, x)     - Computes the second derivative of the IMQ RBF with respect to `x`.
%   D3(r, a, x)     - Computes the third derivative of the IMQ RBF with respect to `x`.
%   D4(r, a, x)     - Computes the fourth derivative of the IMQ RBF with respect to `x`.
%   G(r, a, x, y)   - Computes the gradient of the IMQ RBF (partial derivatives w.r.t. x and y).
%   L(r, a)         - Computes the Laplacian of the IMQ RBF.
%   B(r, a, x, y)   - Computes the biharmonic operator for the IMQ RBF.
%   D12(r, a, x, y) - Computes mixed second-order derivatives (partial derivatives w.r.t. x and y).
%   D22(r, a, x, y) - Computes second-order derivatives in the second dimension.
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% Notes:
%   - The IMQ RBF is smooth and infinitely differentiable, making it suitable for 
%     applications requiring high accuracy.
%   - This implementation is part of the `rbfx` framework and can be extended 
%     for advanced applications such as PDE solving.
%   - See related methods in: Mella et al., "HARP-I: A Harmonic Phase 
%     Interpolation Method for the Estimation of Motion From Tagged MR Images," 
%     IEEE Transactions on Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092

classdef imqx < rbfx
    methods
        function obj = imqx()  
            obj@rbfx();  
          end
         
          function v = rbf(obj,r,a), v = 1./sqrt(1 + (a.*r).^2 ); end
         
          function d = D1(obj,r,a,x), d = -(x.*a.^2)./(1.0 + (a.*r).^2 ).^(1.5); end
         
          function d = D2(obj, r, a, x)
              d = a.^2.*(-r.^2.*a.^2 + 3*a.^2*x.^2 - 1)./(r.^2*a.^2 + 1).^(2.5);
          end
          
          function d = D3(obj, r, a, x) 
              d = 3*a^4*x.*( -5*a.^2.*x.^2 + 3*a.^2.*r.^2 + 3  )./(r.^2.*a.^2 + 1).^(3.5);
          end
          
    
          function d = D4(obj, r, a, x) 
              d = 3*a.^4.*(35*a.^4*x.^4 - 30*a.^2*x.^2.*(a.^2.*r.^2 + 1) + 3*(a.^2.*r.^2 + 1).^2)./(r.^2.*a.^2 + 1).^(4.5); 
    
          end
          
    
          function d = G(obj, r, a, x, y)   % Gradient
             d = -a.^2.*(x + y)./(r.^2.*a.^2 + 1).^(1.5);
          end
          
          function d = L(obj, r, a)         % Laplacian
             d = a.^2.*(r.^2.*a.^2 - 2)./(1 + (a.*r).^2 ).^(2.5);
          end
          
          function d = B(obj, r, a, x, y)   % Biharmonic operator   
             d = 3*a.^4.*(35*a.^4.*x.^4 + 70*a.^4.*x.^2.*y.^2 + 35*a.^4.*y.^4 - 30*a.^2.*x.^2.*(a.^2.*r.^2 + 1) - 30*a.^2.*y.^2.*(a.^2.*r.^2 + 1) - 10*a.^2.*r.^2.*(a.^2.*r.^2 + 1) + 8*(a.^2.*r.^2 + 1).^2)./(1 + (a.*r).^2 ).^(4.5);
          end
    
         
         function d = D12(obj, r, a, x, y) 
             d = 3*x.*a.^4.*(1 - 4*y.^2.*a.^2 + x.^2.*a.^2  )./(1 + (a.*r).^2 ).^(3.5);
         end
         
         function d = D22(obj, r, a, x, y) 
             d = -3*a.^4.*( -1 + 3*y.^2.*a.^2 + 4*x.^4.*a.^4 + 4*y.^4.*a.^4 + x.^2.*(3*a.^2 - 27.*y.^2*a.^4)  )./(1 + (a.*r).^2 ).^(4.5);
         end
      
   end 
end  