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
% GAUSSIAN - Implements the Gaussian Radial Basis Function (RBF)
%
% Description:
%   This class defines the Gaussian RBF and its associated methods for 
%   interpolation. The Gaussian RBF is a smooth, infinitely differentiable 
%   function commonly used in scattered data approximation, machine learning, 
%   and meshless methods. It is part of the `rbfx` framework and inherits 
%   from its superclass.
%
% Syntax:
%   obj = Gaussian();        % Create an instance of the Gaussian RBF class
%   v = obj.rbf(r, a);       % Evaluate the Gaussian RBF
%
% Inputs:
%   r - Distance matrix or distances between points (scalar or array).
%   a - Shape parameter controlling the support size of the RBF.
%
% Outputs:
%   v - RBF values evaluated at the given distances.
%
% Example:
%   % Create a Gaussian RBF object
%   g = Gaussian();
%
%   % Define distances and shape parameter
%   r = [0, 1; 1, 0];  % Distance matrix
%   a = 2.0;           % Shape parameter
%
%   % Compute Gaussian RBF values
%   v = g.rbf(r, a);
%
%   % Display results
%   disp(v);
%
% Methods:
%   rbf(r, a)        - Computes the Gaussian RBF values for a given distance matrix `r` 
%                      and shape parameter `a`.
%   D1(r, a, x)      - Placeholder for the first derivative of the Gaussian RBF 
%                      (currently not implemented, returns an empty array).
%   D2(r, a, x)      - Placeholder for the second derivative of the Gaussian RBF 
%                      (currently not implemented, returns an empty array).
%   D3(r, a, x)      - Placeholder for the third derivative of the Gaussian RBF 
%                      (currently not implemented, returns an empty array).
%   D4(r, a, x)      - Placeholder for the fourth derivative of the Gaussian RBF 
%                      (currently not implemented, returns an empty array).
%   G(r, a, x, y)    - Placeholder for gradient-related calculations 
%                      (currently not implemented, returns an empty array).
%   L(r, a)          - Placeholder for Laplacian-related calculations 
%                      (currently not implemented, returns an empty array).
%   B(r, a, x, y)    - Placeholder for mixed boundary conditions 
%                      (currently not implemented, returns an empty array).
%   D12(r, a, x, y)  - Placeholder for mixed second-order derivatives 
%                      (currently not implemented, returns an empty array).
%   D22(r, a, x, y)  - Placeholder for second-order derivatives in the 
%                      second dimension (currently not implemented, returns an empty array).
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% Notes:
%   - This class implements the Gaussian RBF and is compatible with the 
%     `rbfx` framework. It can be extended with specific derivative 
%     implementations as needed.
%   - See related methods in: Mella et al., "HARP-I: A Harmonic Phase 
%     Interpolation Method for the Estimation of Motion From Tagged MR Images," 
%     IEEE Transactions on Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%

classdef Gaussian < rbfx
    methods
      function obj = Gaussian() 
        obj@rbfx(); 
      end

      function v = rbf(obj,r,a), v = exp(-(r/a).^2); end % Gaussian

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