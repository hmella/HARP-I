function a = rbfShapeSchedule(a_constant, dist, Nfr)
%RBFSHAPESCHEDULE Per-frame RBF shape parameter schedule.
%
%   a = rbfShapeSchedule(a_constant, dist, Nfr)
%
%   A_CONSTANT is either a scalar (constant schedule) or a vector whose
%   min/max bound a schedule of length NFR. DIST selects the shape:
%
%     'Constant'         - a = a_constant (scalar form)
%     'DecreasingLinear' - max -> min
%     'Linear'           - min -> max
%     'Quadratic'        - min + |(max-min)^0.4 * t|^2.5, t in [-1,1]
%
%   The paper (eq. 16-17) uses a single value for every frame; the schedule
%   is an undocumented extension. Matching is case-insensitive and an
%   unrecognised name is an error rather than a silent no-op -- previously
%   the default 'quadratic' matched none of the case-sensitive branches, so
%   `a` kept its input length and a(i) threw for i >= 3.
%
%   See also HARPI, HARPI_3D.

    if isempty(a_constant)
        error([mfilename ':missingShapeParameter'], ...
            ['''a_constant'' is required. For the cubic multiquadric the ' ...
             'paper uses a = 150/(pixels per tag wavelength).']);
    end

    if ~isnumeric(a_constant) || ~all(isfinite(a_constant(:)))
        error([mfilename ':invalidShapeParameter'], ...
            '''a_constant'' must be finite and numeric.');
    end

    if isscalar(a_constant)
        a = a_constant*ones([1 Nfr]);
        return
    end

    valid = {'Constant','DecreasingLinear','Linear','Quadratic'};
    dist  = validatestring(dist, valid, mfilename, 'RBFFacDist');

    lo = min(a_constant(:));
    hi = max(a_constant(:));

    switch dist
        case 'Constant'
            a = lo*ones([1 Nfr]);
        case 'DecreasingLinear'
            a = hi - (hi-lo)*linspace(0,1,Nfr);
        case 'Linear'
            % Ramp lo -> hi. The old expression was hi + (hi-lo)*t, which
            % ran from hi to 2*hi-lo and never visited lo at all.
            a = lo + (hi-lo)*linspace(0,1,Nfr);
        case 'Quadratic'
            a = lo + abs(power(hi-lo,0.4)*linspace(-1,1,Nfr)).^2.5;
    end
end
