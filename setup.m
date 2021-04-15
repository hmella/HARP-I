function setup(varargin)
% SETUP initialize the Cardiac-Motion environment.
%
% Copyright (c) 2016 DENSEanalysis Contributors

    %% EXISITING UTILITY FOLDERS ON SEARCH PATH
    % If a folder on the search path contains certain words, or if certain
    % files are found in the search path, we likely have an existing
    % installation of this software. We may need to remove some folders
    % from the MATLAB search path to avoid potential conflicts.

    % current MATLAB search path
    curdir = regexp(path, ['\' pathsep], 'split')';

    % for comparison, absolute search path
    curdirabs = cell(size(curdir));
    for k = 1:numel(curdir)
        if isdir(curdir{k})
            [stat,info] = fileattrib(curdir{k});
            curdirabs{k} = info.Name;
        end
    end
    path

    % strings & files to locate
    strings = {'common_utilities', 'HARPI_utilities',...
               'HARP_utilities', 'PCSPAMM_utilities',...
               'RBF_utilities', 'Segmentation_utilities',...
               'SinMod_utilities','Test','external_utilities',...
               'HARP_SPHR','PhaseUnwrapping_utilities'};
    files   = {};

    % locate files in search path, translate into strings
    folders = cell(size(files));
    for k = 1:numel(files)
        if exist(files{k},'file')~=0
            F = which(files{k},'-all');
            [P,F,E] = cellfun(@fileparts,F,'uniformoutput',0);
            [P,F,E] = cellfun(@fileparts,P,'uniformoutput',0);
            folders{k} = P(:);
        end
    end
    folders = cat(1,folders{:});

    % all strings to check
    if isempty(strings)
        tmp1 = {};
    else
        tmp1 = strcat('*',strings(:),'*');
    end
    if isempty(folders)
        tmp2 = {};
    else
        tmp2 = strcat(folders(:),'*');
    end

    allstrings = [tmp1;tmp2];
    allstrings = unique(allstrings);

    % search paths similar to candidate toolbox folders
    tf = false(size(curdir));
    for k = 1:numel(allstrings)
        wild = regexptranslate('wildcard',allstrings{k});
        wild = ['^' wild '$'];
        idx  = regexp(curdirabs,wild);
        tf   = tf | cellfun(@(x)~isempty(x),idx);
    end

    % candidate directories to remove
    curdir    = curdir(tf);
    curdirabs = curdirabs(tf);



    %% UPDATE MATLAB SEARCH PATH

    % utility directory
    basedir = fileparts(mfilename('fullpath'));
    utildir = fullfile(basedir, 'utilities');

    % absolute folder path
    [stat,info] = fileattrib(utildir);
    if stat==0
        error(sprintf('%s:directoryNotFound',mfilename),...
            'Utility directory could not be located.');
    end
    utildir = info.Name;

    % locate toolbox directories
    % (one folder down from the utility directory)
    d = dir(utildir);
    d = d(3:end);
    d = d([d.isdir]);
    tooldir = cellfun(@(x)fullfile(utildir,x),...
        {d.name},'uniformoutput',0);


    % remove current directories
    rmdir = setdiff(curdirabs,tooldir);
    if ~isempty(rmdir)
        prompt = {'The following folders on the MATLAB search path';...
            'may be due to previous installations of the DENSEanalysis';...
            'program.  Select paths to remove...'};
        [sel,ok] = listdlg(...
            'PromptString', prompt,...
            'ListString',   rmdir,...
            'ListSize',     [500 300],...
            'Name',         'Remove Search Paths',...
            'InitialValue', 1:numel(rmdir));
        drawnow
        if ~ok
            error(sprintf('%s:setupIncomplete',mfilename),...
                'setup was cancelled.');
        end
        rmpath(rmdir{sel});
    end

    % add utility directories
    if ~isempty(tooldir)
        adddir = setdiff(tooldir,curdir);
        if ~isempty(adddir)
            fprintf(['The following folders will be added to ',...
            'the MATLAB search path:\n'])
            fprintf('%s\n',adddir{:});
        end

        toolbox_dirs = {tooldir{:}};
        for i = 1:numel(toolbox_dirs)
          addpath(genpath(toolbox_dirs{i}));
        end
%         addpath(tooldir{:});
        if ~isempty(adddir)
            fprintf('Search path updated\n')
        end
    end




    %% PARSE ADDITIONAL INPUT OPTIONS

    % permanently save the modified path
    savepath([userpath '/pathdef.m']);
    fprintf('Search path saved\n')
 
    % add startup file to ensure the addition of paths
    % every time that matlab start
    fid = fopen([userpath '/startup.m'],'wt');
    fprintf(fid, 'addpath(pathdef);');
    fclose(fid);

    %% COMPILE MEX FILES

    % compilation options
    if regexp(computer, '64$')
        opts = {'-largeArrayDims'};
    end


    %% CHECK FOR LICENSES


end

