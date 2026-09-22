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
[api, extra_args] = parseinputs(defapi, [], varargin{:});
assertNoUnknownArgs(mfilename, extra_args);

% validate
if isempty(api.Contour) || ~isstruct(api.Contour) && ~isobject(api.Contour)
    error([mfilename ':invalidContour'], ...
        ['''Contour'' must be a cline object (or struct) with a Position ' ...
         'field; got %s. Note the parameter is singular: ''Contour''.'], ...
        class(api.Contour));
end
if isempty(api.MaskSize) || numel(api.MaskSize) < 2
    error([mfilename ':invalidMaskSize'], '''MaskSize'' must be [rows cols].');
end

% contours loop
Ncontour     = numel(api.Contour.Position);
linesegments = cell(1, Ncontour);
masksegments = cell(1, Ncontour);
for i = 1:Ncontour

    pos = api.Contour.Position{1,i};

    % Per-point flags carried by the contour, describing what was drawn.
    isclosed = contour_flag(api.Contour, 'IsClosed', i, true,  1);
    iscurved = contour_flag(api.Contour, 'IsCurved', i, true,  size(pos,1));
    iscorner = contour_flag(api.Contour, 'IsCorner', i, false, size(pos,1));

    % contours segment
    linesegments{i} = clinesegments(pos, isclosed, iscurved, iscorner, ...
                                    api.ContourRes);

    % mask segments
    tmp = linesegments{i};
    masksegments{i} = vertcat(tmp{:,1});

end

% mask
[X, Y] = meshgrid(1:api.MaskSize(2),1:api.MaskSize(1));
mask = maskSA(X,Y,masksegments);

end

% =========================================================================
function val = contour_flag(C, name, idx, default, n)
%CONTOUR_FLAG Fetch a per-contour flag, falling back to a scalar default.
    val = repmat(default, [n 1]);
    if ~isprop(C, name) && ~isfield(C, name)
        return
    end
    raw = C.(name);
    if iscell(raw)
        if idx > numel(raw); return; end
        raw = raw{idx};
    end
    if isempty(raw); return; end
    raw = logical(raw(:));
    if isscalar(raw)
        val = repmat(raw, [n 1]);
    elseif numel(raw) == n
        val = raw;
    end
    if n == 1
        val = val(1);
    end
end

