% Licensing:
%   Copyright (c) 2016-2017 Scott A. Sarra
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
%   This file was originally part of MRBFtoolbox by Scott Sarra.
%   Modifications were made by Hernan Mella on 01/2025.
%
% Notes: 
%    *Distance Matrix Functions*:
%    - New: General distanceMatrixND for N-dimensional spaces.
%    - Original: Adds dimension-specific functions (distanceMatrix1d, distanceMatrix2d, distanceMatrix3d), 
%      with "New" versions using implicit expansion for MATLAB R2016b+.

% RBFX - Abstract class defining the interface for Radial Basis Functions (RBFs).
%
% Description:
%   This abstract class serves as a base class for all radial basis functions (RBF) implementations.
%   It defines common methods required for RBF operations, such as computing the RBF, gradients,
%   Laplacian, and distance matrices. Specific RBF types, such as Gaussian, Multiquadric, and 
%   Wendland, extend this class and provide concrete implementations of the abstract methods.
%
%
% Methods:
%   Abstract:
%     v = rbf(obj, r, a)        - Computes the RBF value.
%     d = D1(obj, r, a, x)      - First derivative with respect to x.
%     d = D2(obj, r, a, x)      - Second derivative with respect to x.
%     d = D3(obj, r, a, x)      - Third derivative with respect to x.
%     d = D4(obj, r, a, x)      - Fourth derivative with respect to x.
%     d = G(obj, r, a, x, y)    - Gradient.
%     d = L(obj, r, a)          - Laplacian.
%     d = B(obj, r, a, x, y)    - Biharmonic operator.
%     d = D12(obj, r, a, x, y)  - Mixed partial derivative.
%     d = D22(obj, r, a, x, y)  - Mixed partial derivative.
%
%   Static:
%     [r, rd] = distanceMatrixNDNew(pointsA, pointsB)
%       - Computes pairwise Euclidean distances and differences in N-dimensional space.
%
%     a = solve(B, f, eta, safe)
%       - Solves a linear system B*a = f with optional regularization (eta) and safety (Cholesky).
%
%     D = dm(B, H, eta, safe)
%       - Computes a derivative matrix using RBF evaluations and optional regularization.
%
%     varargout = variableShape(sMin, sMax, opt, N, M)
%       - Generates variable shape parameters for RBFs based on specified range and strategy.
%
%
classdef rbfx

% ----------------------------------------------------------------------------
% ----------------------- Abstract methods------------------------------------
%             used to define a common interface for all subclasses
%                    must be implemented by all subclasses.
% ----------------------------------------------------------------------------

   methods(Abstract = true)    
       v = rbf(obj,r,a);         % RBF definition
       d = D1(obj,r,a,x);        % first derivative wrt x
       d = D2(obj, r, a, x);     % second derivative wrt x
       d = D3(obj, r, a, x);     % third derivative wrt x
       d = D4(obj, r, a, x);     % fourth derivative wrt x
       d = G(obj, r, a, x, y);   % Gradient
       d = L(obj, r, a);         % Laplacian
       d = B(obj, r, a, x, y);   % Biharmonic operator
       d = D12(obj, r, a, x, y); % mixed partial derivative 
       d = D22(obj, r, a, x, y); % mixed partial derivative        
   end
   
% ----------------------------------------------------------------------------
% ---------------- static methods --------------------------------------------
% ----------------------------------------------------------------------------
    
    methods(Static)
         
% ----------------------------------------------------------------------------
% ------------------ distance matrices ---------------------------------------
% ----------------------------------------------------------------------------

% ----------------------------------------------------------------------------
% DISTANCEMATRIX2DNEW - Computes the 2D Euclidean distance matrix.
%
% Description:
%   This function calculates pairwise Euclidean distances between points in 
%   2D space. It can compute distances either:
%     1) Within a single set of points.
%     2) Between two different sets of points.
%
% Syntax:
%   [r, rx, ry] = distanceMatrix2dNew(xc, yc)
%   [r, rx, ry] = distanceMatrix2dNew(xc, yc, x, y)
%
% Inputs:
%   xc, yc  - [vectors] Coordinates of the first set of points.
%   x, y    - [vectors, optional] Coordinates of the second set of points.
%             If not provided, distances are computed within the first set 
%             of points.
%
% Outputs:
%   r       - [matrix] Euclidean distance matrix. 
%             If only (xc, yc) are provided, size(r) = [N x N].
%             If (x, y) are also provided, size(r) = [M x N], where N and M
%             are the number of points in the first and second sets, respectively.
%   rx      - [matrix] Differences in the x-coordinates.
%   ry      - [matrix] Differences in the y-coordinates.
%
% Example:
%   % Compute distances within a single set of points
%   xc = [1; 2; 3];
%   yc = [4; 5; 6];
%   [r, rx, ry] = distanceMatrix2dNew(xc, yc);
%
%   % Compute distances between two sets of points
%   x = [7; 8];
%   y = [9; 10];
%   [r, rx, ry] = distanceMatrix2dNew(xc, yc, x, y);
%
%   % Display the distance matrix
%   disp(r);
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (Benjamin.lopezf@usm.cl)
%
% Notes:
%   - Computes both the Euclidean distance matrix `r` and the coordinate 
%     differences (`rx` for x-coordinates, `ry` for y-coordinates).
%   - Efficient for small to medium-sized datasets in 2D.

    function [r, rx, ry] = distanceMatrix2dNew(xc,yc,x,y)
        xc = xc(:);   yc = yc(:);
        if nargin==2
            rx = xc - xc';
            ry = yc - yc'; 
            r = sqrt( rx.^2 + ry.^2 ); 
        else
            x = x(:);   y = y(:);
            rx = x - xc';
            ry = y - yc'; 
            r = sqrt( rx.^2 + ry.^2 ); 
        end
    end

% DISTANCEMATRIX3DNEW - Computes the 3D Euclidean distance matrix.
%
% Description:
%   This function calculates pairwise Euclidean distances between points in 
%   3D space. It can compute distances either:
%     1) Within a single set of points.
%     2) Between two different sets of points.
%
% Syntax:
%   [r, rx, ry, rz] = distanceMatrix3dNew(xc, yc, zc)
%   [r, rx, ry, rz] = distanceMatrix3dNew(xc, yc, zc, x, y, z)
%
% Inputs:
%   xc, yc, zc - [vectors] Coordinates of the first set of points.
%   x, y, z    - [vectors, optional] Coordinates of the second set of points.
%                If not provided, distances are computed within the first set 
%                of points.
%
% Outputs:
%   r          - [matrix] Euclidean distance matrix. 
%                If only (xc, yc, zc) are provided, size(r) = [N x N].
%                If (x, y, z) are also provided, size(r) = [M x N], where N and M
%                are the number of points in the first and second sets, respectively.
%   rx         - [matrix] Differences in the x-coordinates.
%   ry         - [matrix] Differences in the y-coordinates.
%   rz         - [matrix] Differences in the z-coordinates.
%
% Example:
%   % Compute distances within a single set of points
%   xc = [1; 2; 3];
%   yc = [4; 5; 6];
%   zc = [7; 8; 9];
%   [r, rx, ry, rz] = distanceMatrix3dNew(xc, yc, zc);
%
%   % Compute distances between two sets of points
%   x = [10; 11];
%   y = [12; 13];
%   z = [14; 15];
%   [r, rx, ry, rz] = distanceMatrix3dNew(xc, yc, zc, x, y, z);
%
%   % Display the distance matrix
%   disp(r);
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (Benjamin.lopezf@usm.cl)
%
% Notes:
%   - Computes both the Euclidean distance matrix `r` and the coordinate 
%     differences (`rx` for x-coordinates, `ry` for y-coordinates, `rz` for z-coordinates).
%   - Efficient for small to medium-sized datasets in 3D.

    function [r, rx, ry, rz] = distanceMatrix3dNew(xc,yc,zc,x,y,z)
        xc = xc(:);   yc = yc(:);  zc = zc(:);
        if nargin==3
            rx = xc - xc';
            ry = yc - yc';
            rz = zc - zc';
            r = sqrt( rx.^2 + ry.^2 + rz.^2 );
        else
            x = x(:);   y = y(:);  z = z(:);
            rx = x - xc';
            ry = y - yc';
            rz = z - zc';
            r = sqrt( rx.^2 + ry.^2 + rz.^2 ); 
        end
    end

% ----------------------------------------------------------------------------
% ---------- regularized SPD linear system solvers ---------------------------
% ----------------------------------------------------------------------------
       
% solve - solves the SPD linear system B a = f for a with the option to regularize
%         by the method of diagonal increments (MDI)
% inputs 
%    B    N x N symmetric positive definite (SPD) matrix
%    f    N x 1 vector
%   eta    (optional) MDI regularization parameter.  Use eta = 0 for no regularization
%  safe   true - uses backslash with error checking etc.
%         false - uses a Cholesky factorization.  Faster, but it the matrix is 
%                 severely ill-conditioned and/or the regularization parameter is too
%                 small the matrix may fail to be numerically SPD and the Cholesky
%                 factorization will fail
%
% outputs
%    a    N x 1 solution vector
%
% example usage:
%   1) mdiExample.m, 2) mdiRegularization.m, 3) rbfInterpConvergence.m

    function a = solve(B, f, eta, safe) 
        if ~exist('eta','var'), eta = 5e-15;  end
            if ~exist('safe','var'), safe = true;  end
        
            % Si eta>0, se suma eta en la diagonal de B => B = B + eta*I
            if eta > 0
                N = length(f);          % Número de filas en f (o Nx1)
                % (eta*norm(B,2) a veces se usa para escalar la regularización 
                %  según la magnitud de B, pero aquí se deja comentado)
                B(1:N+1:end) = B(1:N+1:end) + eta;
            end
            
            % Modo "seguro" o directo => B\f
            if safe
                a = B \ f;
            else
                % Modo "no seguro": factoriza B con Cholesky, 
                % asumiendo B ~ simétrica p.d. para mayor eficiencia.
                L = chol(B,'lower');
                a = L'\( L\f );
            end
        end
            
% -------------------------------------------------------------------------       
 
% dm - forms the deriviative matrix D = H*inv(B) by solving the system D B = H for D
% inputs
%    B   N x N SPD system matrix
%    H   N x N derivative evaluation matrix
%   eta    (optional) MDI regularization parameter.  Use eta = 0 for no regularization
%  safe   true - uses backslash with error checking etc.
%         false - uses a Cholesky factorization.  Faster, but if the matrix is 
%                 severely ill-conditioned and/or the regularization parameter is too
%                 small the matrix may fail to be numerically SPD and the Cholesky
%                 factorization will fail
%
% outputs
%    D   N x N differentiation matrix
%
% example usage:
%  1) diffusionReactionCentro.m

    function D = dm(B,H,eta,safe)      
        if ~exist('eta','var'), eta = 5e-15;  end
        if ~exist('safe','var'), safe = true;  end
        
        if eta>0
            s = size(B); N = s(1);
            B(1:N+1:end) = B(1:N+1:end) + eta;  % B = B + eta*eye(N);
        end
          
        if safe
            D = H/B; 
        else 
            L = chol(B,'lower');
            D = (L'\(L\H'))';   
       end
    end
       
% ----------------------------------------------------------------------------   
% -------------- RBF interp---------------------------------------------------
% ----------------------------------------------------------------------------

% RBFINTERP2D - Performs Radial Basis Function (RBF) interpolation in 2D.
%
% Description:
%   This function interpolates 2D data points using Radial Basis Functions (RBF). 
%   Given a set of input coordinates (`x`), their corresponding values (`y`), 
%   and query points (`x0`), the function computes interpolated values (`xe`) 
%   at the query points using the specified RBF kernel.
%
% Syntax:
%   [xe, re] = RBFInterp2D(x, y, p, s, x0, phi)
%
% Inputs:
%   x   - [2 x N] Coordinates of the N control points in 2D space.
%   y   - [C x N] Values associated with each of the N control points 
%                 (C channels, e.g., scalar or vector data).
%   p   - Smoothing factor for regularization (scalar).
%   a   - Scale parameter for the RBF kernel (scalar or vector).
%   x0  - [2 x M] Coordinates of the M query points.
%   phi - Struct containing the RBF kernel function handle (e.g., phi.rbf).
%
% Outputs:
%   xe  - [M x C] Interpolated values at the M query points for each of the C channels.
%   re  - [M x N] Pairwise distance matrix between query points and control points.
%
% Example:
%   % Define control points and their values
%   x = [1 2 3; 10 12 14];  % [2 x N]
%   y = [5 9 10];           % [1 x N]
%   x0 = [1.5 2.5; 11 13];  % [2 x M]
%
%   % Define RBF kernel (Gaussian)
%   phi.rbf = @(r, s) exp(-(r./s).^2);  
%   a = 1; p = 1e-8;
%
%   % Perform RBF interpolation
%   [xe, re] = RBFInterp2D(x, y, p, s, x0, phi);
%   disp('Interpolated values:'), disp(xe);
%
% Author:
%   Hernán Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% Licensing:
%   This Source Code Form is subject to the terms of the Mozilla Public License, 
%   version 2.0. If a copy of the MPL was not distributed with this file, You can 
%   obtain one at http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - This function uses `distanceMatrix2dNew` to compute Euclidean distances.
%   - Regularization via the `p` parameter ensures numerical stability.
%   - The RBF interpolation methodology is related to the techniques described 
%     in Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical 
%     Imaging, 2021.

    function [xe, re] = RBFInterp2D(x,y,p,a,x0,phi)
        
        % RBFinterp performs a RBF interpolation of the data (x,y)
        % at x0, using phi as RBF and p as smoothing factor  
    
        % Distance matrices
        r  = rbfx.distanceMatrix2dNew(x(1,:),x(2,:));
        re = rbfx.distanceMatrix2dNew(x(1,:),x(2,:),...
                x0(1,:),x0(2,:));
    
        % RBF evaluation and weights estimation
        B = phi.rbf(r,a);
        w = rbfx.solve(B,y',p,true);
    
        % RBF evaluation at x0
        Psi = phi.rbf(re,a);
    
        % Interpolated values
        xe = Psi*w;
    
    end

% RBFINTERP3D - Performs Radial Basis Function (RBF) interpolation in 3D.
%
% Description:
%   This function interpolates 3D data points using Radial Basis Functions (RBF). 
%   Given a set of input coordinates (`x`), their corresponding values (`y`), 
%   and query points (`x0`), the function computes interpolated values (`xe`) 
%   at the query points using the specified RBF kernel.
%
% Syntax:
%   [xe, re] = RBFInterp3D(x, y, p, a, x0, phi)
%
% Inputs:
%   x   - [3 x N] Coordinates of the N control points in 3D space.
%   y   - [C x N] Values associated with each of the N control points 
%                 (C channels, e.g., scalar or vector data).
%   p   - Smoothing factor for regularization (scalar).
%   a   - Scale parameter for the RBF kernel (scalar or vector).
%   x0  - [3 x M] Coordinates of the M query points.
%   phi - Struct containing the RBF kernel function handle (e.g., phi.rbf).
%
% Outputs:
%   xe  - [M x C] Interpolated values at the M query points for each of the C channels.
%   re  - [M x N] Pairwise distance matrix between query points and control points.
%
% Example:
%   % Define control points and their values
%   x = [1 2 3; 10 12 14; 5 6 7]; % [3 x N]
%   y = [5 9 10];                 % [1 x N]
%   x0 = [1.5 2.5; 11 13; 6 6.5]; % [3 x M]
%
%   % Define RBF kernel (Gaussian)
%   phi.rbf = @(r, s) exp(-(r./s).^2);  
%   a = 1; p = 1e-8;
%
%   % Perform RBF interpolation
%   [xe, re] = RBFInterp3D(x, y, p, a, x0, phi);
%   disp('Interpolated values:'), disp(xe);
%
% Author:
%   Hernán Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% Licensing:
%   This Source Code Form is subject to the terms of the Mozilla Public License, 
%   version 2.0. If a copy of the MPL was not distributed with this file, You can 
%   obtain one at http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - This function uses `distanceMatrix3dNew` to compute Euclidean distances in 3D space.
%   - Regularization via the `p` parameter ensures numerical stability.
%   - The RBF kernel function (`phi.rbf`) defines the type of interpolation, such as 
%     Gaussian, multiquadric, or thin plate spline.
%   - The RBF interpolation methodology is related to the techniques described 
%     in Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical 
%     Imaging, 2021.
%

    function [xe, re] = RBFInterp3D(x,y,p,a,x0,phi)
        
        % RBFinterp performs a RBF interpolation of the data (x,y)
        % at x0, using phi as RBF and p as smoothing factor  

        % Distance matrices
        r  = rbfx.distanceMatrix3dNew(x(1,:),x(2,:),x(3,:));
        re = rbfx.distanceMatrix3dNew(x(1,:),x(2,:),x(3,:),...
                x0(1,:),x0(2,:),x0(3,:));

        % RBF evaluation and weights estimation
        B = phi.rbf(r,a);
        w = rbfx.solve(B,y',p,true);

        % RBF evaluation at x0
        Psi = phi.rbf(re,a);

        % Interpolated values
        xe = Psi*w;

    end

% ----------------------------------------------------------------------------   
% -------------- variable shape parameters -----------------------------------
% ----------------------------------------------------------------------------

% variableShape 
%
% example usage: variableShapeInterp1d.m
%
%  inputs
%   sMin    minimum value of the shape parameter
%   sMax    maximum value of the shape parameter
%    N      number of columns
%    M      number of rows
%  opt  1, exponentially varying (Kansa)
%          Computers and Mathematics with Applications v. 24, no. 12, 1992. 
%       2, linearly varying
%       3, randonly varying (Sarra and Sturgil)
%          Engineering Analysis with Boundary Elements, v. 33, p. 1239-1245, 2009.
%
%  outputs
%    s1   N x N matrix with constant shapes in each column
%         call as, s1 = rbfx.variableShape(sMin,sMax,N)
%    s2   M x N matrix with constant shapes in each column 
%         (optional, for interpolation evaluation matrix) 
%          call as, [s1, s2] = rbfx.variableShape(sMin,sMax,N,M)
%
% example usage:
%   1) variableShapeInterp1d.m

    function varargout = variableShape(sMin,sMax,opt,N,M)
        if nargin<5, M = []; end
        nOutputs = nargout;
        varargout = cell(1,nOutputs);

        if opt==1          
            sMin = sMin^2;      sMax = sMax^2;
            s = sqrt( sMin*(sMax/sMin).^((0:N-1)./(N-1)) );      
        elseif opt==2   
            s = sMin + ((sMax - sMin)/(N-1)).*(0:N-1);
        else     
            s = rand(1,N);
            s = sMin + (sMax - sMin)*s; 
        end
        
        if nOutputs==2
            varargout{1} = repmat(s,N,1);  
            varargout{2} = repmat(s,M,1);
        else
           varargout{1} = repmat(s,N,1);  
        end
    end 
    



end 
    
   
end 