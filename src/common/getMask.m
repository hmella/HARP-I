% Copyright (c) 2025 Hernan Mella
%
% GETMASK - Generates a binary mask based on provided contours.
%
% Description:
%   This function creates a binary mask that delimits the area defined by
%   the input contours. The mask can be used for applications such as image
%   segmentation or spatial analysis.
%
% Syntax:
%   mask = getMask('ParameterName', ParameterValue, ...)
%
% Inputs:
%   'MaskSize'   - Dimensions of the binary mask [rows, columns].
%   'Contour'    - Structure containing the contour information:
%                  - Position: Coordinates of the contour (Nx2 array).
%   'ContourRes' - (Optional) Resolution of the contour segments. 
%                  Default: 0.5.
%
% Outputs:
%   mask         - Binary mask (matrix) with the same size as specified 
%                  by MaskSize. The mask contains `1` inside the contours
%                  and `0` outside.
%
% Example:
%   % Define mask dimensions
%   maskSize = [256, 256];
%
%   % Example contour
%   contour.Position = {[50, 50; 200, 50; 200, 200; 50, 200; 50, 50]};
%
%   % Generate binary mask
%   mask = getMask('MaskSize', maskSize, 'Contour', contour, 'ContourRes', 0.1);
%
%   % Visualize the mask
%   imshow(mask);
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
%   - This implementation is based on the methods described in the paper:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical 
%     Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%   - The function is a computational realization of the contour-based segmentation 
%     methods detailed in the aforementioned publication.
%
% Notes:
%   - This implementation aligns with methods described in:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on 
%     Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%

function [mask] = getMask(varargin)
% Default arguments
defapi = struct(...
        'MaskSize',           [],...
        'Contour',            [],...
        'ContourRes',         0.5);

% check inputs
api = parseinputs(defapi, [], varargin{:});

% contours loop
linesegments = {};
masksegments = {};
for i = 1:numel(api.Contour.Position)

    % contours segment
    linesegments{i} = clinesegments(api.Contour.Position{1,i},...
        true,true(size(api.Contour.Position{1,i})),...
        false,api.ContourRes);

    % mask segments
    tmp = linesegments{i};
    masksegments{i} = vertcat(tmp{:,1});

end

% mask
[X, Y] = meshgrid(1:api.MaskSize(2),1:api.MaskSize(1));
mask = maskSA(X,Y,masksegments);

return;

end

