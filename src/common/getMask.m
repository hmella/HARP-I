function [mask] = getMask(varargin)
%GETMASK Summary of this function goes here
%   Detailed explanation goes here

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

