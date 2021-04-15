function [strain] = getStrainBySegments(varargin)
%SegmentStrain estimates the strain for each segment of 
% a SA view of the left ventricle.
%
%

    % Default arguments
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
            'PositionB',        []);

    % Check input
    api = parseinputs(defapi, [], varargin{:});

    % image size
    Isz = size(api.Mask);

    % number of frames
    Nfr = size(api.CC,3);

    % parse coordinates
    if or(isempty(api.X),isempty(api.Y))
        [X,Y] = meshgrid(1:Isz(1),1:Isz(2));
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
        api.PositionB = origin(1);
    end

    % primary spoke angle
    angle0 = atan2(api.PositionB-origin(2),api.PositionA-origin(1));

    % minor spoke angles
    if api.ClockWise
        angles = linspace(0,2*pi,api.Nseg+1);
    else
        angles = linspace(2*pi,0,api.Nseg+1);
    end
    angles = angles + angle0 - pi;
    
    % initialize output
    strain = struct(...
        'segments',      NaN(Isz),...
        'segments_CC',   NaN([api.Nseg,Nfr]),...
        'segments_RR',   NaN([api.Nseg,Nfr]));

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

        % strain
        for fr=1:Nfr
            if ~isempty(api.CC)
                CCfr = api.CC(:,:,fr);
                strain.segments_CC(seg,fr) = mean(CCfr(and(~(CCfr>0.5),pos)),'omitnan');
            end
            if ~isempty(api.RR)
                RRfr = api.RR(:,:,fr);
                strain.segments_RR(seg,fr) = mean(RRfr(and(~(RRfr>0.5),pos)),'omitnan');
            end
        end

    end

    % segments
    segments(~api.Mask) = NaN;

    % output
    strain.segments = segments;

    return;

end