% Copyright (c) 2025 Hernan Mella
%
% REFPHASESMOOTHING3D - Smooths three 3D phase images (Xpha, Ypha, Zpha) using RBF interpolation
%
% Description:
%   This function applies radial basis function (RBF) interpolation to smooth
%   three phase channels in a 3D grid (Xpha, Ypha, Zpha). It constructs a 3D
%   coordinate mesh and performs RBF-based reconstruction of the phase values,
%   allowing for noise reduction or continuity improvements in multi-channel,
%   volumetric phase data.
%
% Syntax:
%   [Xphas, Yphas, Zphas] = RefPhaseSmoothing3D(Xpha, Ypha, Zpha, method, RBFFactor, eta)
%
% Inputs:
%   Xpha       - 3D array (size [dimY x dimX x dimZ]) representing the first phase channel.
%   Ypha       - 3D array (same size as Xpha) representing the second phase channel.
%   Zpha       - 3D array (same size as Xpha) representing the third phase channel.
%   method     - Method or handle specifying the RBF kernel (e.g., 'ThinPlate').
%   RBFFactor  - Scale factor (scalar) for the RBF, influencing the radial spread or stiffness.
%   eta     - Smoothing/regularization parameter passed to the RBF solver.
%
% Outputs:
%   Xphas      - Smoothed version of Xpha (same 3D size as Xpha).
%   Yphas      - Smoothed version of Ypha (same 3D size as Ypha).
%   Zphas      - Smoothed version of Zpha (same 3D size as Zpha).
%
% Example:
%   % Suppose Xpha, Ypha, and Zpha are three 3D phase volumes of size [128 x 128 x 32]
%   Xpha = rand(128, 128, 32);
%   Ypha = rand(128, 128, 32);
%   Zpha = rand(128, 128, 32);
%   method    = 'ThinPlate'; 
%   RBFFactor = 150;
%   eta    = 1e-08;
%
%   % Smooth the 3D reference phases
%   [Xphas, Yphas, Zphas] = RefPhaseSmoothing3D(Xpha, Ypha, Zpha, method, RBFFactor, eta);
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% License:
%   This Source Code Form is subject to the terms of the Mozilla Public License,
%   version 2.0. If a copy of the MPL was not distributed with this file,
%   you can obtain one at http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - The 3D mask of valid points is determined where all three volumes 
%     (Xpha, Ypha, Zpha) are non-NaN.
%   - The smoothed output preserves the original 3D dimensions, but replaces 
%     valid voxels with the interpolated phase values.
%   - Implements the techniques described in:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical 
%     Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%   - Note: Temporal fitting is currently a placeholder in this implementation 
%     and has not been fully developed for 3D data.
%

function [Xphas, Yphas, Zphas] = RefPhaseSmoothing3D(Xpha,Ypha,Zpha,method,RBFFactor,eta)

    % Create grid
    Isz = size(Xpha,1:3);
    [X, Y, Z] = meshgrid(1:Isz(2),1:Isz(1),1:Isz(3));

    % Get radial basis functions
    phi = feval(method);
    s = RBFFactor;

    % Smooth phase using RBF interpolations
    tf = ~isnan(Xpha) & ~isnan(Ypha) & ~isnan(Zpha);
    fint = rbfx.RBFInterp3D([X(tf), Y(tf), Z(tf)]',[Xpha(tf), Ypha(tf), Zpha(tf)]',...
            eta,s(1),[X(tf), Y(tf), Z(tf)]',phi);
    Xphas = NaN(Isz);
    Yphas = NaN(Isz);
    Zphas = NaN(Isz);
    Xphas(tf) = fint(:,1);
    Yphas(tf) = fint(:,2);
    Zphas(tf) = fint(:,3);

end