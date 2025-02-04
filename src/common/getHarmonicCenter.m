% Copyright (c) 2025 Hernan Mella
%
% GETHARMONICCENTER - Determines the harmonic center in k-space.
%
% Description:
%   This function identifies the harmonic center in k-space by analyzing
%   a windowed region of the frequency domain. The harmonic center is 
%   the point of maximum intensity within the specified region.
%
% Syntax:
%   c = getHarmonicCenter(K, omega, window, FOV)
%   c = getHarmonicCenter(K, omega, window, FOV, center)
%
% Inputs:
%   K       - Complex k-space data (2D matrix).
%   omega   - Frequency vector [omega_x, omega_y] specifying the harmonic 
%             frequency.
%   window  - Half-dimensions of the windowed region in k-space [rows, cols].
%   FOV     - Field of view of the image in the spatial domain [FOV_x, FOV_y].
%   center  - (Optional) Initial guess for the harmonic center [row, col]. 
%             If not provided, it is calculated based on omega and FOV.
%
% Outputs:
%   c       - Coordinates of the harmonic center [row, col].
%
% Example:
%   % Simulated k-space data
%   K = rand(128, 128) + 1i * rand(128, 128); % Complex k-space
%
%   % Frequency and window parameters
%   omega = [0.3, 0.2];
%   window = [5, 5];
%   FOV = [256, 256];
%
%   % Find the harmonic center
%   harmonic_center = getHarmonicCenter(K, omega, window, FOV);
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
%   - This implementation aligns with methods described in:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on 
%     Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%

function [c] = getHarmonicCenter(K, omega, window, FOV, center)
    % Window center
    if nargin < 5
        % Kspace center
        s = round(0.5*size(K));
        c = round(s + omega.*FOV/(2*pi));
    else
        c = center;
    end
        
    % Indices for windowed K-space
    row0 = c(2)-window(2); row1 = c(2)+window(2);
    col0 = c(1)-window(1); col1 = c(1)+window(1);
    
    % Extract windowed K-space
    k = K(row0:row1,col0:col1);
    
    % Look for maximun
    k = abs(k);
    v = 0.0;
    for row = 1:size(k,1)
        for col = 1:size(k,2)
            if k(row,col) > v
                v = k(row,col);
                ii = row;
                jj = col;
            end
        end
    end

    % Local to global coefficients
    I = row0 + ii - 1;
    J = col0 + jj - 1;
    c = [I,J];

end
