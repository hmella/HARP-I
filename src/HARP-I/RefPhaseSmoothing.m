% Copyright (c) 2025 Hernan Mella
%
% REFPHASESMOOTHING - Smooths two phase images (Xpha, Ypha) using RBF interpolation
%
% Description:
%   This function applies radial basis function (RBF) interpolation to smooth
%   two phase channels (Xpha and Ypha). It constructs a coordinate grid and
%   performs an RBF-based reconstruction of the phase values, preserving
%   spatial consistency while reducing noise or small discontinuities.
%
% Syntax:
%   [Xphas, Yphas] = RefPhaseSmoothing(Xpha, Ypha, method, a_constant, eta)
%
% Inputs:
%   Xpha       - 2D array representing the first phase channel.
%   Ypha       - 2D array (same size as Xpha) representing the second phase channel.
%   method     - Method or handle specifying the RBF kernel (e.g., 'ThinPlate').
%   a_constant  - Scale factor (scalar) for the RBF; can influence the spread or stiffness. Equation (16) and (17) on the paper.
%   eta     - Smoothing/regularization parameter passed to the RBF solver.
%
% Outputs:
%   Xphas      - Smoothed version of Xpha (same size as Xpha).
%   Yphas      - Smoothed version of Ypha (same size as Ypha).
%
% Example:
%   % Suppose Xpha and Ypha are two phase images of size [256 x 256]
%   Xpha = rand(256, 256);
%   Ypha = rand(256, 256);
%   method = 'ThinPlate';     % Or a function handle that returns phi.rbf
%   a_constant = 150;
%   eta = 1e-08;
%
%   % Smooth the reference phases
%   [Xphas, Yphas] = RefPhaseSmoothing(Xpha, Ypha, method, a_constant, eta);
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (Benjamin.lopezf@usm.cl)
%
% License:
%   This Source Code Form is subject to the terms of the Mozilla Public License,
%   version 2.0. If a copy of the MPL was not distributed with this file,
%   you can obtain one at http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - The mask of valid points is determined where both Xpha and Ypha are non-NaN.
%   - After interpolation, smoothed values are written back into arrays of
%     the same size as the input phases.
%   - Implements the techniques described in:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical 
%     Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%   - Note: Temporal fitting is currently a placeholder in this implementation 
%     and has not been fully developed for 3D data.

function [Xphas, Yphas] = RefPhaseSmoothing(Xpha,Ypha,method,a,eta)

    % Create grid
    Isz = size(Xpha);
    [X, Y] = meshgrid(1:Isz(2),1:Isz(1));

    % Get radial basis functions
    phi = feval(method);

    % Smooth phase using RBF interpolations
    tf = ~isnan(Xpha) & ~isnan(Ypha);
    fint = rbfx.RBFInterp2D([X(tf), Y(tf)]',[Xpha(tf), Ypha(tf)]',...
            eta,a(1),[X(tf), Y(tf)]',phi);
    Xphas = NaN(Isz(1:2));
    Yphas = NaN(Isz(1:2));
    Xphas(tf) = fint(:,1);
    Yphas(tf) = fint(:,2);

end