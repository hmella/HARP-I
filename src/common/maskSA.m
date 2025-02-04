% Copyright (c) 2025 Hernan Mella
%
% MASKSA - Create a logical mask for points inside or on the borders of polygons.
%
% Syntax:
%   tf = maskSA(X, Y, C)
%
% Description:
%   This function computes a logical mask for points in a 2D plane, identifying
%   points that are within the outer polygon but not the inner polygon, or that
%   lie exactly on the borders of either polygon.
%
% Inputs:
%   X, Y - Coordinates of the points to evaluate. Can be vectors or matrices.
%   C    - Cell array containing two polygons:
%          C{1} - Coordinates of the outer polygon (Nx2 matrix).
%          C{2} - Coordinates of the inner polygon (Mx2 matrix).
%
% Output:
%   tf   - Logical array of the same size as X and Y. True for points that:
%          - Lie inside the outer polygon but outside the inner polygon.
%          - Lie on the borders of either polygon.
%
% Example:
%   % Define outer and inner polygons
%   C{1} = [0, 0; 5, 0; 5, 5; 0, 5; 0, 0]; % Outer square
%   C{2} = [2, 2; 4, 2; 4, 4; 2, 4; 2, 2]; % Inner square
%
%   % Define grid points
%   [X, Y] = meshgrid(0:0.5:5, 0:0.5:5);
%
%   % Compute logical mask
%   tf = maskSA(X, Y, C);
%
%   % Plot result
%   scatter(X(tf), Y(tf), 'r', 'filled'); % Points in mask
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% License:
%   This Source Code Form is subject to the terms of the Mozilla Public 
%   License, v. 2.0. If a copy of the MPL was not distributed with this 
%   file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - Inspired by concepts for geometric masking in Mella et al., "HARP-I: A Harmonic
%     Phase Interpolation Method for the Estimation of Motion From Tagged MR Images,"
%     IEEE Transactions on Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%
% Notes:
%   - This implementation aligns with methods described in:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on 
%     Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%

function tf = maskSA(X,Y,C)

    [inep,onep] = inpolygon(X,Y,C{1}(:,1),C{1}(:,2));
    [inen,onen] = inpolygon(X,Y,C{2}(:,1),C{2}(:,2));
    tf = (inep & ~inen) | onep | onen;
end