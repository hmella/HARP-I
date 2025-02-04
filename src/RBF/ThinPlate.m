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
% THINPLATE - Implements the Thin Plate Radial Basis Function (RBF)
%
% Description:
%   This class defines the Thin Plate RBF, a popular choice for scattered data 
%   interpolation due to its smoothness and ability to model higher-order deformations. 
%   The Thin Plate RBF is defined as:
%       φ(r) = r^2 * log(r), where r is the distance.
%
% Syntax:
%   obj = ThinPlate();        % Create an instance of the Thin Plate RBF class
%   v = obj.rbf(r, a);        % Evaluate the Thin Plate RBF
%
% Inputs:
%   r - Distance matrix or distances between points (scalar or array).
%   a - Unused parameter for compatibility with the abstract base class `rbfx`.
%
% Outputs:
%   v - RBF values evaluated at the given distances.
%
% Example:
%   % Create an instance of the Thin Plate RBF class
%   phi = ThinPlate();
%
%   % Define distances
%   r = linspace(0.1, 2, 100); % Distance matrix (avoiding r = 0 to prevent log(0))
%
%   % Evaluate the Thin Plate RBF
%   v = phi.rbf(r, []);
%
%   % Plot the RBF
%   plot(r, v);
%   xlabel('Distance (r)');
%   ylabel('RBF Value');
%   title('Thin Plate Radial Basis Function');
%
% Methods:
%   rbf(r, a)       - Computes the Thin Plate RBF values for a given distance matrix `r`.
%                     The parameter `a` is unused but retained for compatibility.
%   D1(r, a, x)     - Placeholder for the first derivative of the Thin Plate RBF 
%                     (returns an empty array).
%   D2(r, a, x)     - Placeholder for the second derivative of the Thin Plate RBF 
%                     (returns an empty array).
%   D3(r, a, x)     - Placeholder for the third derivative of the Thin Plate RBF 
%                     (returns an empty array).
%   D4(r, a, x)     - Placeholder for the fourth derivative of the Thin Plate RBF 
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
%   - The Thin Plate RBF is suitable for interpolation tasks requiring smoothness 
%     and is particularly effective in 2D.
%   - The parameter `a` is included for compatibility but is not used in this implementation.
%   - Avoid `r = 0` to prevent undefined behavior due to the logarithm term log(r).
%   - See related methods in: Mella et al., "HARP-I: A Harmonic Phase Interpolation 
%     Method for the Estimation of Motion From Tagged MR Images," IEEE Transactions 
%     on Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092


classdef ThinPlate < rbfx
    methods
        function obj = ThinPlate()  
            obj@rbfx();  
          end
    
          function v = rbf(obj,r,a)
              v = r.^2.*log(r);
              v(isnan(v)) = 0;
          end % Thin plate
    
          function d = D1(obj,r,a,x), d = []; end
         
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
