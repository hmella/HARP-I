% Copyright (c) 2025 Hernan Mella
%
% DOWNSAMPLE - Reduces the resolution of an N-dimensional array.
%
% Description:
%   This function downsamples an input array `A` along specified dimensions 
%   by a given factor. The user can optionally specify a shift to control 
%   the starting point of the downsampling process.
%
% Syntax:
%   B = Downsample(A, dims, factor)
%   B = Downsample(A, dims, factor, shift)
%
% Inputs:
%   A        - Input array (N-dimensional).
%   dims     - Dimensions along which the array will be downsampled (vector).
%   factor   - Downsampling factor (integer).
%   shift    - (Optional) Shift to control the starting index for downsampling. 
%              Default: 0.
%
% Outputs:
%   B        - Downsampled array.
%
% Example:
%   % Create a 3D array
%   A = rand(100, 100, 50);
%
%   % Downsample along the first and second dimensions by a factor of 2
%   B = Downsample(A, [1, 2], 2);
%
%   % Downsample along the first dimension with a factor of 3 and a shift of 1
%   B = Downsample(A, 1, 3, 1);
%
%   % Visualize the size of the downsampled array
%   disp(size(B));
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopex (Benjamin.lopezf@usm.cl)
%
% License:
%   This Source Code Form is subject to the terms of the Mozilla Public License, 
%   version 2.0. If a copy of the MPL was not distributed with this file, 
%   you can obtain one at http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - This function implements spatial downsampling based on subsampling indices.
%   - It allows for flexible selection of dimensions and customizable shifts, 
%     making it suitable for applications in image processing, data compression, 
%     and spatial resolution adjustment.
%   - Related Section in Paper:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical 
%     Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.

function B = Downsample(A,dims,factor,shift)
    % Default shift is zero
    if nargin < 4
        shift = 0;
    end

    % Number of dimensions of the array
    Ndims = numel(size(A));

    % Array dimensions
    a_dims = 1:Ndims;
    a_dims(dims) = [];

    % Indices
    idx = cell([1 numel(dims)]);
    for i=dims
        idx{i} = (1+shift):factor:size(A,i);
    end
    for i=a_dims
        idx{i} = ':';
    end

    % Donwsample array
    S.subs = idx;
    S.type = '()';
    B = subsref(A,S);

end
