% Copyright (c) 2025 Hernan Mella
%
% WRAP - Wraps a phase to its principal range [-pi, pi].
%
% Syntax:
%   W = Wrap(phase)
%   W = Wrap(phase, row, col)
%
% Description:
%   This function adjusts phase values to ensure they fall within the 
%   principal range of [-pi, pi]. If row and column indices are provided,
%   it returns the value at the specified location.
%
% Inputs:
%   phase - Matrix of phase values (numeric).
%   row   - (Optional) Row index (positive integer).
%   col   - (Optional) Column index (positive integer).
%
% Outputs:
%   W - Matrix of phase values wrapped to the range [-pi, pi], or a scalar 
%       value corresponding to the specified indices.
%
% Examples:
%   % Basic example:
%   phase = [-3*pi, pi; 0, 3*pi];
%   W = Wrap(phase);
%
%   % Example with indices:
%   W = Wrap(phase, 1, 2);
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
%   Adapted for illustrative purposes in motion estimation.
%   Refer to Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for
%   the Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical Imaging, 
%   vol. 40, no. 4, pp. 1240-1251, April 2021, for methods inspiring this functionality.
%

function [W] = Wrap(phase, row, col)
    W = mod(phase + pi, 2 * pi) - pi;
    if nargin > 1
        W = W(row, col);
    end
end
