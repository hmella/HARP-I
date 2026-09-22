function out = capture_reference(varargin)
%CAPTURE_REFERENCE Snapshot every stage of the HARP-I 2D pipeline.
%
%   capture_reference()
%   capture_reference('Output', '/path/ref.mat', 'Data', '/path/data_2D.mat')
%   out = capture_reference(...)
%
%   Runs the demo_2D pipeline with plotting disabled and saves every stage
%   boundary -- not just the final result -- so that a later difference can be
%   attributed to one specific change. The same file is consumed by the Python
%   port's test suite, hence the -v7.3 format (readable by h5py).
%
%   Always run single-threaded for reproducibility:
%       matlab -batch "maxNumCompThreads(1); capture_reference"
%
%   See also COMPARE_REFERENCE.

    root = fileparts(fileparts(mfilename('fullpath')));
    root = fileparts(root);   % .../HARP-I

    p = struct('Output', '', 'Data', fullfile(root, 'data', 'data_2D.mat'), ...
               'Save', true);
    for k = 1:2:numel(varargin)
        p.(validatestring(varargin{k}, fieldnames(p))) = varargin{k+1};
    end

    if isempty(p.Output)
        p.Output = fullfile(root, 'tools', 'regression', 'refs', ...
                            sprintf('ref_%s.mat', git_sha(root)));
    end

    %% ---- data ----------------------------------------------------------
    S    = load(p.Data);
    data = S.data;
    I    = data.image;
    Isz  = size(I);
    Nfr  = Isz(4);
    mask = logical(data.mask);

    pxsz = data.pixelsize/1000;
    s    = data.tagspacing/1000;
    ke   = 2*pi/s*pxsz;

    %% ---- stage 1: harmonic filter --------------------------------------
    filt = HARPFilter(struct('Image', I, 'CentralFreq', ke, ...
        'Direction', deg2rad([0 90]), 'FilterType', 'Transmission', ...
        'Butterworth_cuttoff', 17, 'Butterworth_order', 10));
    Ih = filt.filter(I);

    %% ---- stage 2: wrapped harmonic phase --------------------------------
    Xpha = squeeze(angle(Ih(:,:,1,:)));
    Ypha = squeeze(angle(Ih(:,:,2,:)));

    %% ---- stage 3: phase quality + unwrapping ----------------------------
    qual_x = zeros(size(Xpha));
    qual_y = zeros(size(Ypha));
    uXpha  = NaN(size(Xpha));
    uYpha  = NaN(size(Ypha));
    for i = 1:Nfr
        m = mask(:,:,i);
        qual_x(:,:,i) = phasequality2_mask(Xpha(:,:,i), [1 1], m, 4);
        qual_y(:,:,i) = phasequality2_mask(Ypha(:,:,i), [1 1], m, 4);
        uXpha(:,:,i)  = unwrap2(Xpha(:,:,i), 'Mask', m, 'PixelSize', [1 1], ...
                                'Seed', 'auto', 'Connectivity', 4);
        uYpha(:,:,i)  = unwrap2(Ypha(:,:,i), 'Mask', m, 'PixelSize', [1 1], ...
                                'Seed', 'auto', 'Connectivity', 4);
    end

    %% ---- stage 4: temporal phase consistency ----------------------------
    % Captured at BOTH jump budgets: njumps=1 is what HARPI effectively uses
    % today (the nargin off-by-one silently clamps it), njumps=8 is the
    % intended behaviour. Keeping both makes that fix verifiable.
    phase_jumps_n1 = TemporalPhaseConsistency(uXpha, uYpha, 1, 1, []);
    phase_jumps_n8 = TemporalPhaseConsistency(uXpha, uYpha, 1, 8, [1 Isz(1) 1 Isz(2)]);

    %% ---- stage 5: reference phase smoothing ------------------------------
    a = 150/(s/pxsz(1));
    [Xphas, Yphas] = RefPhaseSmoothing(uXpha(:,:,1), uYpha(:,:,1), ...
                                       'Multiquadric3', a/5, 15);

    %% ---- stage 6: motion estimation --------------------------------------
    args = struct( ...
        'Mask',                 mask, ...
        'Frames',               1:Nfr, ...
        'Show',                 false, ...
        'Method',               'Multiquadric3', ...
        'a_constant',           [1 1]*a, ...
        'RBFFacDist',           'Linear', ...
        'eta_constant',         15, ...
        'Connectivity',         4, ...
        'RefPhaseSmoothing',    true, ...
        'TemporalFitting',      true, ...
        'TemporalFittingOrder', 10);
    t_harpi = tic;
    harpi   = HARPI(Ih, args);
    t_harpi = toc(t_harpi);

    dxr = squeeze(harpi.FittedMotion(:,:,1,:));
    dyr = squeeze(harpi.FittedMotion(:,:,2,:));

    %% ---- stage 7: strain ---------------------------------------------------
    [X, Y] = meshgrid(1:Isz(2), 1:Isz(1));
    st = pixelstrain(struct('X', X, 'Y', Y, 'mask', mask(:,:,1), ...
        'times', 1:Nfr, 'dx', dxr, 'dy', dyr, 'Origin', [], 'Orientation', []));

    %% ---- stage 8: global curves --------------------------------------------
    GCS = mean(st.CC, 1, 'omitnan')';
    GRS = mean(st.RR, 1, 'omitnan')';

    %% ---- assemble ----------------------------------------------------------
    out = struct();
    out.Ih             = Ih;
    out.Xpha           = Xpha;
    out.Ypha           = Ypha;
    out.qual_x         = qual_x;
    out.qual_y         = qual_y;
    out.uXpha          = uXpha;
    out.uYpha          = uYpha;
    out.phase_jumps_n1 = phase_jumps_n1;
    out.phase_jumps_n8 = phase_jumps_n8;
    out.Xphas          = Xphas;
    out.Yphas          = Yphas;
    out.u_raw          = harpi.RawMotion;
    out.u_fit          = harpi.FittedMotion;
    out.traj_raw       = harpi.RawTrajectories;
    out.traj_fit       = harpi.FittedTrajectories;
    out.unwrapped      = harpi.UnwrappedPhases;
    out.st             = st;
    out.GCS            = GCS;
    out.GRS            = GRS;
    out.meta           = struct( ...
        'git_sha',     git_sha(root), ...
        'matlab',      version, ...
        'threads',     maxNumCompThreads, ...
        'timestamp',   datestr(now, 'yyyy-mm-ddTHH:MM:SS'), ...
        'data',        p.Data, ...
        't_harpi',     t_harpi);

    %% ---- save ---------------------------------------------------------------
    if p.Save
        outdir = fileparts(p.Output);
        if ~isfolder(outdir); mkdir(outdir); end
        save(p.Output, '-struct', 'out', '-v7.3');
        fprintf('reference written: %s\n', p.Output);
        fprintf('  HARPI time %.2f s, threads %d, sha %s\n', ...
                t_harpi, maxNumCompThreads, out.meta.git_sha);
    end
end

% =========================================================================
function sha = git_sha(root)
    sha = 'unknown';
    [st, o] = system(sprintf('git -C "%s" rev-parse --short HEAD', root));
    if st == 0
        sha = strtrim(o);
        [~, d] = system(sprintf('git -C "%s" status --porcelain -- src HARPIsetup.m', root));
        if ~isempty(strtrim(d)); sha = [sha '-dirty']; end
    end
end
