% Copyright (c) 2016 DENSEanalysis Contributors
%
% PARSEINPUTS - Parse and validate input arguments with defaults.
%
% Syntax:
%   [valid_args, other_args] = parseinputs(valid_param, valid_default, varargin)
%
% Description:
%   This function processes variable input arguments (varargin), identifying
%   and validating predefined parameters while maintaining defaults. Parameters
%   not in the predefined list are returned as additional arguments.
%
% Inputs:
%   valid_param   - Cell array or structure of valid parameter names.
%   valid_default - Default values corresponding to valid_param.
%   varargin      - Variable input arguments (parameter/value pairs or structures).
%
% Outputs:
%   valid_args    - Structure containing valid parameters with user-specified or default values.
%   other_args    - Additional arguments not in valid_param, returned as parameter/value pairs.
%
% Example:
%   % Define valid parameters
%   valid_param = {'DisplayRange', 'Colormap'};
%   valid_default = {[0 1], 'gray'};
%
%   % Call parseinputs
%   [valid_args, other_args] = parseinputs(valid_param, valid_default, ...
%       'DisplayRange', [0 255], 'SomeOtherParam', 42);
%
%   % Output:
%   % valid_args.DisplayRange -> [0 255]
%   % valid_args.Colormap -> 'gray'
%   % other_args -> {'SomeOtherParam', 42}
%
% Author:
%   Drew Gilliam
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% License:
%   This Source Code Form is subject to the terms of the Mozilla Public
%   License, v. 2.0. If a copy of the MPL was not distributed with this
%   file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%
% Notes:
%   - Inspired by practices in DENSEanalysis software for handling flexible function inputs.
%   - Referenced: Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for
%     the Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical Imaging,
%     vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092.
%

function [valid_args, other_args] = parseinputs...
    (valid_param, valid_default, varargin)
    % check for proper number of inputs
    narginchk(2, Inf);

    % default argument structure
    if isstruct(valid_param)
        valid_args = valid_param;
        valid_param = fieldnames(valid_args);
    else
        valid_args = cell2struct(valid_default(:)', valid_param(:)', 2);
    end

    % return on empty
    if nargin == 2
        other_args = {};
        return;
    end

    % create empty "other_args" structure
    other_args = struct;

    % parse input arguments
    Nvar = numel(varargin);
    k = 1;
    while k <= Nvar

        % parse structure
        if isstruct(varargin{k})

            tags = fieldnames(varargin{k});
            for ti = 1:numel(tags)
                tf = strcmpi(tags{ti},valid_param);
                if any(tf)
                    valid_args.(valid_param{tf}) = varargin{k}.(tags{ti});
                else
                    tag = tags{ti};
                    other_args.(tag) = varargin{k}.(tags{ti});
                end
            end

            k = k+1;

        % parse param/value pair
        elseif ischar(varargin{k}) && k+1 <= Nvar

            tf = strcmpi(varargin{k},valid_param);
            if any(tf)
                valid_args.(valid_param{tf}) = varargin{k+1};
            else
                other_args.(varargin{k}) = varargin{k+1};
            end
            k = k+2;

        % unrecognized input
        else
           error(sprintf('%s:unrecognizedInput',mfilename),...
               'Param/Value pairs not as expected.');
        end

    end

    % separate "other_args" into param/value pairs
    C = [fieldnames(other_args), struct2cell(other_args)]';
    other_args = C(:)';

end

