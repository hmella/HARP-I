% Copyright (c) 2025 Hernan Mella
%
% GETSTRAINBYSEGMENTS - Computes strain values for specific radial and circumferential segments of a region.
%
% Description:
%   This function divides a region of interest into angular segments based on a defined origin and orientation.
%   It calculates mean strain values for each segment in both circumferential (CC) and radial (RR) directions
%   across multiple frames.
%
% Syntax:
%   strain = getStrainBySegments('CC', CC, 'RR', RR, 'Mask', mask, 'Origin', origin, 
%                                'Orientation', theta, 'Nseg', Nseg, 'ClockWise', ClockWise)
%
% Inputs:
%   CC            - (Optional) Circumferential strain matrix (size: [rows, cols, frames]).
%   RR            - (Optional) Radial strain matrix (size: [rows, cols, frames]).
%   Mask          - Binary mask defining the region of interest (ROI) within the image. Default: entire image.
%   Origin        - (Optional) Coordinates [x, y] for the center of segmentation. Default: geometric center of the mask.
%   Orientation   - (Optional) Angular orientation of the region in polar coordinates. Default: calculated using `cart2pol`.
%   Nseg          - (Optional) Number of angular segments to divide the region into. Default: 6.
%   ClockWise     - (Optional) Boolean flag to define clockwise or counterclockwise segmentation. Default: false (counterclockwise).
%   PositionA     - (Optional) X-coordinate defining the primary spoke angle. Default: maximum X-coordinate.
%   PositionB     - (Optional) X-coordinate of the origin. Default: origin x-coordinate.
%
% Outputs:
%   strain - A structure containing:
%       - `segments`       : 2D matrix where each pixel is assigned to a segment (NaN outside the mask).
%       - `segments_CC`    : Mean circumferential strain for each segment across all frames (size: [Nseg, frames]).
%       - `segments_RR`    : Mean radial strain for each segment across all frames (size: [Nseg, frames]).
%
% Example:
%   % Example: Divide a circular mask into 6 segments and compute strain
%   CC = randn(128, 128, 10); % Circumferential strain (10 frames)
%   RR = randn(128, 128, 10); % Radial strain (10 frames)
%   mask = createCircularMask(128, 128); % Binary mask
%   origin = [64, 64]; % Center of the mask
%   strain = getStrainBySegments('CC', CC, 'RR', RR, 'Mask', mask, 'Origin', origin, 'Nseg', 6);
%
% Key Features:
%   - Automatically calculates the angular orientation of the segments based on the defined origin and geometry.
%   - Divides the region of interest into angular segments for detailed strain analysis.
%   - Computes mean strain values for both circumferential (CC) and radial (RR) components.
%   - Handles NaN values in the input strain matrices.
%   - Allows clockwise or counterclockwise segmentation.
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% Notes:
%   - The polynomial fitting ensures smooth trajectories while preserving temporal consistency.
%   - Boundary conditions are enforced using first and second derivatives at the initial and final frames.
%   Refer to Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for
%   the Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical Imaging, 
%   vol. 40, no. 4, pp. 1240-1251, April 2021, for methods inspiring this functionality.
%  - This based on the code of DENSEanalysis: segmentmodel.m (https://github.com/denseanalysis/denseanalysis/blob/main/DENSE_utilities/DENSE_toolbox/segmentmodel.m)
%
% Notes:
%   - This implementation aligns with methods described in:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on 
%     Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%



function [strain] = getStrainBySegments(varargin)
    defapi = struct(...
            'CC',               [],...
            'RR',               [],...
            'X',                [],...
            'Y',                [],...
            'Mask',             [],...
            'Origin',           [],...
            'Orientation',      [],...
            'Nseg',             6,...
            'ClockWise',        false,...
            'PositionA',        [],...
            'PositionB',        [],...
            'Frames',           [],...
            'OutlierThreshold', 0.5);

    % Check input
    [api, extra_args] = parseinputs(defapi, [], varargin{:});
    assertNoUnknownArgs(mfilename, extra_args);

    errid = sprintf('%s:invalidInput', mfilename);

    % image size
    Isz = size(api.Mask);

    % number of frames (CC may be omitted when only RR is wanted)
    if ~isempty(api.CC)
        Nfr = size(api.CC,3);
    elseif ~isempty(api.RR)
        Nfr = size(api.RR,3);
    else
        error(errid, 'At least one of CC or RR must be supplied.');
    end
    if ~isempty(api.Frames)
        if ~isempty(api.CC); api.CC = api.CC(:,:,api.Frames); end
        if ~isempty(api.RR); api.RR = api.RR(:,:,api.Frames); end
        Nfr = numel(api.Frames);
    end

    % parse coordinates
    if isempty(api.X) || isempty(api.Y)
        [X,Y] = meshgrid(1:Isz(2),1:Isz(1));
    else
        X = api.X;
        Y = api.Y;
    end
    if ~isequal(size(X), Isz) || ~isequal(size(Y), Isz)
        error(errid, 'X and Y must be the same size as Mask (%s).', mat2str(Isz));
    end

    % parse origin
    origin = api.Origin;
    if isempty(origin)
        origin = [mean(X(api.Mask)), mean(Y(api.Mask))];
    end

    if ~isnumeric(origin) || size(origin,2)~=2
        error(errid,'Invalid Origin parameter.');
    end

    % parse orientation
    theta = api.Orientation;
    if isempty(theta)
        theta = cart2pol(X-origin(1),Y-origin(2));
    end

    if ~isnumeric(theta) || ~all(size(theta)==Isz)
        error(errid,'%s','Invalid Orientation.');
    end

    % throw warning
    if ~isempty(api.Orientation) && ~isempty(api.Origin)
        warning(sprintf('%s:ignoringParameter',mfilename),'%s',...
            'It is unneccesary to enter an Origin value ',...
            'when externally defining all Orientation values.');
    end

    % parse positions
    if isempty(api.PositionA)
        api.PositionA = max(X(:));
    end

    if isempty(api.PositionB)
        api.PositionB = origin(2);   % y of the primary spoke -> angle0 == 0
    end

    % primary spoke angle
    angle0 = atan2(api.PositionB-origin(2),api.PositionA-origin(1));

    % minor spoke angles
    if api.ClockWise
        angles = linspace(0,2*pi,api.Nseg+1);
    else
        angles = linspace(2*pi,0,api.Nseg+1);
    end
    angles = angles - pi;

    % Work in the rotated frame, so every pixel lands in exactly one bin.
    theta = mod(theta - angle0 + pi, 2*pi) - pi;
    
    % initialize output
    strain = struct(...
        'segments',      NaN(Isz),...
        'segments_CC',   NaN([api.Nseg,Nfr]),...
        'segments_RR',   NaN([api.Nseg,Nfr]));

    % Slice each frame once rather than once per segment.
    CCfr = cell(1,Nfr); RRfr = cell(1,Nfr);
    for fr = 1:Nfr
        if ~isempty(api.CC); CCfr{fr} = api.CC(:,:,fr); end
        if ~isempty(api.RR); RRfr{fr} = api.RR(:,:,fr); end
    end

    % prepare output
    segments = NaN(Isz);
    for seg=1:api.Nseg

        % positions
        if api.ClockWise
            pos = and(theta > angles(seg), theta <= angles(seg+1));
        else
            pos = and(theta <= angles(seg), theta >= angles(seg+1));
        end

        % segments
        segments(pos) = seg;

        % strain, restricted to the mask
        inseg = pos & api.Mask;
        for fr=1:Nfr
            if ~isempty(api.CC)
                sel = inseg & ~(CCfr{fr} > api.OutlierThreshold);
                strain.segments_CC(seg,fr) = mean(CCfr{fr}(sel),'omitnan');
            end
            if ~isempty(api.RR)
                sel = inseg & ~(RRfr{fr} > api.OutlierThreshold);
                strain.segments_RR(seg,fr) = mean(RRfr{fr}(sel),'omitnan');
            end
        end

    end

    % segments
    segments(~api.Mask) = NaN;

    % output
    strain.segments = segments;

    return;

end