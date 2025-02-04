% Copyright (c) 2025 Hernan Mella
%
% FLATTEN - Reshapes a matrix into a single row or column vector.
%
% Description:
%   This function takes a matrix `A` and converts it into either a row or
%   column vector, based on the specified `format`. By default, the output 
%   is a column vector.
%
% Syntax:
%   B = flatten(A)
%   B = flatten(A, format)
%
% Inputs:
%   A       - Input matrix of any size.
%   format  - (Optional) Output format:
%             - 'row': Converts the matrix to a row vector.
%             - 'col': Converts the matrix to a column vector (default).
%
% Outputs:
%   B       - Reshaped vector (row or column).
%
% Examples:
%   % Convert a matrix to a column vector (default)
%   A = [1, 2; 3, 4];
%   B = flatten(A);
%   % Output: B = [1; 3; 2; 4]
%
%   % Convert a matrix to a row vector
%   A = [1, 2; 3, 4];
%   B = flatten(A, 'row');
%   % Output: B = [1, 3, 2, 4]
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

function B = flatten(A,format)
  if nargin < 2
    format = 'col';
  end

  if strcmp(format,'row')
    B = A(:)';
  elseif strcmp(format,'col')
    B = A(:);
  end

end