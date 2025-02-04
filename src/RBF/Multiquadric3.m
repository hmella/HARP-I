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
% MULTIQUADRIC3 - Implements the Multiquadric3 Radial Basis Function (RBF)
%
% Description:
%   This class defines the Multiquadric3 RBF and its associated methods. The 
%   Multiquadric3 RBF is defined as:
%       φ(r) = (r^2 + a^2)^(3/2),
%   where `r` is the distance and `a` is the shape parameter controlling the scaling.
%   The Multiquadric3 RBF introduces an exponent of 3/2 for smoother interpolation.
%
% Syntax:
%   obj = Multiquadric3();     % Create an instance of the Multiquadric3 RBF class
%   v = obj.rbf(r, a);         % Evaluate the Multiquadric3 RBF
%
% Inputs:
%   r - Distance matrix or distances between points (scalar or array).
%   a - Shape parameter controlling the scaling of the RBF.
%   x, c - Coordinates for derivative calculations.
%
% Outputs:
%   v - RBF values evaluated at the given distances.
%
% Example:
%   % Create an instance of the Multiquadric3 RBF class
%   phi = Multiquadric3();
%
%   % Define distances and shape parameter
%   r = linspace(0, 2, 100);  % Distance matrix
%   a = 1.0;                  % Shape parameter
%
%   % Evaluate the Multiquadric3 RBF
%   v = phi.rbf(r, a);
%
%   % Plot the RBF
%   plot(r, v);
%   xlabel('Distance (r)');
%   ylabel('RBF Value');
%   title('Multiquadric3 Radial Basis Function');
%
% Methods:
%   rbf(r, a)       - Computes the Multiquadric3 RBF values for a given distance matrix `r` 
%                     and shape parameter `a`.
%   D1(r, a, x, c)  - Computes the first derivative of the Multiquadric3 RBF with respect 
%                     to `x` and centers `c`.
%   D2(r, a, x)     - Placeholder for the second derivative of the Multiquadric3 RBF 
%                     (returns an empty array).
%   D3(r, a, x)     - Placeholder for the third derivative of the Multiquadric3 RBF 
%                     (returns an empty array).
%   D4(r, a, x)     - Placeholder for the fourth derivative of the Multiquadric3 RBF 
%                     (returns an empty array).
%   G(r, a, x, y)   - Placeholder for gradient-related calculations (returns an empty array).
%   L(r, a)         - Placeholder for Laplacian-related calculations (returns an empty array).
%   B(r, a, x, y)   - Placeholder for biharmonic operator calculations (returns an empty array).
%   D12(r, a, x, y) - Placeholder for mixed second-order derivatives (returns an empty array).
%   D22(r, a, x, y) - Placeholder for second-order derivatives in the second dimension 
%                     (returns an empty array).
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% Notes:
%   - The Multiquadric3 RBF introduces an additional smoothing factor with an 
%     exponent of 3/2, making it suitable for applications requiring higher-order smoothness.
%   - This implementation is part of the `rbfx` framework and serves as an advanced 
%     RBF example.
%   - See related methods in: Mella et al., "HARP-I: A Harmonic Phase Interpolation 
%     Method for the Estimation of Motion From Tagged MR Images," IEEE Transactions 
%     on Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092


classdef Multiquadric3 < rbfx
    methods
      function obj = Multiquadric3()  
        obj@rbfx();  
      end

      function v = rbf(obj,r,a), v = (r.^2 + a.^2).^(3/2); end % Multiquadratic3

      function d = D1(obj,r,a,x,c), d = 3*(x - c').*(r.^2 + a.^2).^(1/2); end
     
      function d = D2(obj, r, a, x), d = []; end
      
      function d = D3(obj, r, a, x), d = []; end
      
      function d = D4(obj, r, a, x), d = []; end
      
      function d = G(obj, r, a, x, y), d = []; end
      
      function d = L(obj, r, a), d = []; end
      
      function d = B(obj, r, a, x, y), d = []; end

     
     function d = D12(obj, r, a, x, y), d = []; end
     
     function d = D22(obj, r, a, x, y), d = []; end

    end 
end  
