function info = harpi_paths(varargin)
%HARPI_PATHS Put HARP-I and its DENSEanalysis dependency on the MATLAB path.
%
%   harpi_paths()
%   harpi_paths('DENSEanalysis', '/path/to/denseanalysis')
%   info = harpi_paths(...)
%
%   Adds DENSEanalysis first and HARP-I's src/ last, so HARP-I's own
%   functions take precedence. DENSEanalysis ships an incompatible
%   PIXELSTRAIN and its own PARSEINPUTS; if those win the search path,
%   demo_2D.m fails with "Invalid spline data".
%
%   Any remaining shadowed function names are printed, and returned in
%   INFO.collisions when an output is requested.
%
%   See also HARPIsetup.

    % ---- parse -----------------------------------------------------------
    p = struct('DENSEanalysis', '', 'Verbose', true);
    for k = 1:2:numel(varargin)
        name = validatestring(varargin{k}, fieldnames(p));
        p.(name) = varargin{k+1};
    end

    root = fileparts(fileparts(mfilename('fullpath')));

    % ---- locate DENSEanalysis -------------------------------------------
    dense = p.DENSEanalysis;
    if isempty(dense)
        dense = locate_denseanalysis(root);
    end

    % ---- add DENSEanalysis first, HARP-I last ---------------------------
    if ~isempty(dense) && isfolder(dense)
        utils = fullfile(dense, 'DENSE_utilities');
        if isfolder(utils)
            addpath(genpath(utils));
        else
            addpath(genpath(dense));
        end
    elseif p.Verbose
        warning('harpi_paths:noDENSEanalysis', ...
            ['DENSEanalysis was not found. Phase unwrapping (unwrap2/unwrap3) and\n' ...
             'the segmentation GUI (cline/imcline/getcline) will be unavailable.\n' ...
             'Clone it from https://github.com/denseanalysis/denseanalysis and pass\n' ...
             'its location: harpi_paths(''DENSEanalysis'', ''/path/to/denseanalysis'')']);
    end

    addpath(genpath(fullfile(root, 'src')));   % HARP-I wins any tie

    % ---- report collisions ----------------------------------------------
    collisions = check_collisions(root, p.Verbose);

    if nargout > 0
        info = struct('root', root, 'denseanalysis', dense, 'collisions', {collisions});
    end
end

% =========================================================================
function dense = locate_denseanalysis(root)
% Look for DENSEanalysis on the existing path, then in sibling directories.
    dense = '';

    w = which('unwrap2');
    if ~isempty(w)
        % .../denseanalysis/DENSE_utilities/unwrap_toolbox/unwrap2.m
        dense = fileparts(fileparts(fileparts(w)));
        return
    end

    parent = fileparts(root);
    candidates = {fullfile(parent, 'denseanalysis'), ...
                  fullfile(parent, 'DENSEanalysis'), ...
                  fullfile(root, 'denseanalysis')};
    for k = 1:numel(candidates)
        if isfolder(fullfile(candidates{k}, 'DENSE_utilities'))
            dense = candidates{k};
            return
        end
    end
end

% =========================================================================
function collisions = check_collisions(root, verbose)
% Report HARP-I functions that are shadowed by another copy on the path.
%
% Builds one name -> files map from every directory on the path, which is much
% faster than calling WHICH('-all') per function.

    srcdir = fullfile(root, 'src');

    dirs = regexp(path, pathsep, 'split');
    dirs = dirs(~cellfun(@isempty, dirs));

    names  = {};
    owners = {};
    for d = 1:numel(dirs)
        listing = dir(fullfile(dirs{d}, '*.m'));
        for k = 1:numel(listing)
            [~, fname] = fileparts(listing(k).name);
            names{end+1}  = fname;                               %#ok<AGROW>
            owners{end+1} = fullfile(dirs{d}, listing(k).name);  %#ok<AGROW>
        end
    end

    mine = dir(fullfile(srcdir, '**', '*.m'));
    mynames = cell(1, numel(mine));
    for k = 1:numel(mine)
        [~, mynames{k}] = fileparts(mine(k).name);
    end

    collisions = {};
    for k = 1:numel(mynames)
        hits = owners(strcmp(names, mynames{k}));
        if numel(hits) > 1 && ~strncmp(hits{1}, srcdir, numel(srcdir))
            collisions{end+1} = struct('name', mynames{k}, ...
                'winner', hits{1}, 'all', {hits}); %#ok<AGROW>
        end
    end

    if verbose && ~isempty(collisions)
        fprintf(2, '\n*** HARP-I path collision: %d function(s) shadowed ***\n', ...
                numel(collisions));
        for k = 1:numel(collisions)
            fprintf(2, '  %s -> %s\n', collisions{k}.name, collisions{k}.winner);
        end
        fprintf(2, ['Re-run harpi_paths() to restore the correct ordering.\n' ...
                    'HARP-I must come after DENSEanalysis on the path.\n\n']);
    end
end
