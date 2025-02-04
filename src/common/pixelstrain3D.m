% Copyright (c) 2016 DENSEanalysis Contributors
%
% PIXELSTRAIN3D Calculate 3D strain tensors based on displacement fields.
%
% Description:
%   This function computes strain tensors (e.g., radial, circumferential, and
%   principal strains) at each voxel location within a 3D segmented volume mask. 
%   It extends the functionality of the original `pixelstrain` implementation 
%   to three dimensions by introducing an additional spatial dimension (Z).
%
%   The function calculates strain values by determining the deformation
%   gradient tensor at each voxel, applying coordinate system rotations,
%   and computing principal strains and orientations.
%
% Inputs:
%   'X'                - [3D matrix] Grid of X-coordinates.
%   'Y'                - [3D matrix] Grid of Y-coordinates.
%   'Z'                - [3D matrix] Grid of Z-coordinates.
%   'mask'             - [3D matrix] Logical mask indicating regions of interest.
%   'times'            - [vector] Time points for displacement measurements.
%   'dx'               - [4D matrix] Precomputed X-displacements over time.
%   'dy'               - [4D matrix] Precomputed Y-displacements over time.
%   'dz'               - [4D matrix] Precomputed Z-displacements over time.
%   'Origin'           - [vector] Origin for strain orientation calculations.
%                         Default: mean(X(mask)), mean(Y(mask)), mean(Z(mask)).
%   'Orientation'      - [3D matrix] Voxel-level orientation angles (radians).
%                         Default: calculated based on `Origin`.
%
% Outputs:
%   strain             - [struct] Structure containing the following fields:
%                         - maskimage: Processed mask for visualization.
%                         - Other strain components can be extended similarly.
%
%
% Autor:
%   DENSEanalysis: - Andrew Gilliam <http://www.adgilliam.com>
%                  - Andrew Scott <https://github.com/andydscott>
%                  - David vanMaanen <https://github.com/dpvanmaan>
%                  - Jonathan Suever <https://suever.com>
%
%   Copyright (c) 2016 DENSEanalysis Contributors
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
%
% Differences from `pixelstrain`:
%   1. Includes an additional Z-dimension for inputs (`Z`, `dz`) and computations.
%   2. Neighbor calculations and mask validations are extended to 3D.
%   3. Orientation and strain computations consider voxel-level 3D coordinates.
%
% Licensing:
%   This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
%   If a copy of the MPL was not distributed with this file, You can obtain one at
%   http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - This implementation is a 3D extension of the 2D `pixelstrain` function.
%   - The added Z-dimension allows for the analysis of volumetric datasets.
%   - The function requires displacement fields (`dx`, `dy`, `dz`) to be 
%     precomputed and supplied as inputs.

function strain = pixelstrain3D(varargin)

    
    %% PARSE APPLICATION DATA

    errid = sprintf('%s:invalidInput',mfilename);
    
    % parse input data
    defapi = struct(...
        'X',                [],...
        'Y',                [],...
        'Z',                [],...
        'mask',             [],...
        'times',            [],...
        'dx',               [],...
        'dy',               [],...
        'dz',               [],...
        'Origin',           [],...
        'Orientation',      []);
    api = parseinputs(defapi,[],varargin{:});
    
    % check X/Y/Z/mask
    X = api.X;
    Y = api.Y;
    Z = api.Z;    
    mask = api.mask;

    % check times
    time = api.times;
    if ~isnumeric(time) || ~ismatrix(time)
        error(errid,'Invalid times.');
    end

    % additional parameters
    Ntime = numel(time);
    Isz   = size(mask);


    % pixel trajectories
    % note this effectively checks the spline data as well
    xtrj = NaN([Isz,Ntime]);
    ytrj = xtrj;
    ztrj = xtrj;    

    if ~isempty(api.dx) && ~isempty(api.dy) && ~isempty(api.dz)
        for k = 1:Ntime
            xtrj(:,:,:,k) = X + api.dx(:,:,:,k);
            ytrj(:,:,:,k) = Y + api.dy(:,:,:,k);
            ztrj(:,:,:,k) = Z + api.dz(:,:,:,k);
        end
    end

    % parse origin
    origin = api.Origin;
    if isempty(origin)
        origin = [mean(X(mask)), mean(Y(mask))];
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

    % cos/sin calculation (saving computational effort)
    ct = cos(theta);
    st = sin(theta);

    % eliminate any invalid mask locations
    h = [0 1 0; 1 0 0; 0 0 0];
    tmp_SA = false([Isz,4]);
    tmp_LA = false([Isz,4]);    
    for i=1:Isz(3)
        for k = 1:4
            tmp_SA(:,:,i,k) = conv2(double(mask(:,:,i)),h,'same')==2;
            h = rot90(h);
        end
    end
    for i=1:Isz(2)
        for k = 1:4
            tmp_LA(:,i,:,k) = conv2(squeeze(double(mask(:,i,:))),h,'same')==2;
            h = rot90(h);
        end
    end
    mask = any(tmp_SA,4) & any(tmp_LA,4) & mask;
    clear tmp_SA tmp_LA
    

    %% STRAIN CALCULATION

    % determine number of neighbors
    h = [0 1 0; 1 0 1; 0 1 0];
    Nneighbor_SA = zeros(Isz);
    for i=1:Isz(3)
        Nneighbor_SA(:,:,i) = mask(:,:,i).*conv2(double(mask(:,:,i)),h,'same');
    end
    Nneighbor_LA = zeros(Isz);
    for i=1:Isz(2)
        Nneighbor_LA(:,i,:) = squeeze(mask(:,i,:)).*conv2(squeeze(double(mask(:,i,:))),h,'same');
    end    

    % initialize output strain structure
    tmp = NaN([Isz Ntime]);

    strain = struct(...
        'maskimage',    [],...
        'RR',       tmp,...
        'RC',       tmp,...
        'RL',       tmp,...
        'CR',       tmp,...
        'CC',       tmp,...
        'CL',       tmp,...
        'LR',       tmp,...
        'LC',       tmp,...
        'LL',       tmp);
    
    % strain calculation at each point
    dx = zeros(3,6);
    dX = zeros(3,6);
    tf = false(3,6);

    for fr = 1:Ntime
        fprintf('\n Estimating strain on frame %d',fr)
        for k = 1:Isz(3)
            for j = 1:Isz(2)
                for i = 1:Isz(1)

                    if mask(i,j,k)

                        dx(:) = 0;
                        dX(:) = 0;
                        tf(:) = false;

                        % X direction
                        if (i-1 >= 1) && mask(i-1,j,k) && Nneighbor_SA(i,j,k) > 1
                            tf(1) = true;
                            dx(:,1) = [xtrj(i-1,j,k,fr) - xtrj(i,j,k,fr);...
                                       ytrj(i-1,j,k,fr) - ytrj(i,j,k,fr);...
                                       ztrj(i-1,j,k,fr) - ztrj(i,j,k,fr)];
                            dX(:,1) = [X(i-1,j,k) - X(i,j,k); ...
                                       Y(i-1,j,k) - Y(i,j,k); ...
                                       1e-10];
                        end

                        if (i+1 <= Isz(1)) && mask(i+1,j,k) && Nneighbor_SA(i,j,k) > 1
                            tf(2) = true;
                            dx(:,2) = [xtrj(i+1,j,k,fr) - xtrj(i,j,k,fr);...
                                       ytrj(i+1,j,k,fr) - ytrj(i,j,k,fr);...
                                       ztrj(i+1,j,k,fr) - ztrj(i,j,k,fr)];
                            dX(:,2) = [X(i+1,j,k) - X(i,j,k); ...
                                       Y(i+1,j,k) - Y(i,j,k); ...
                                       1e-10];
                        end

                        % Y direction
                        if (j-1 >= 1) && mask(i,j-1,k) && Nneighbor_SA(i,j,k) > 1
                            tf(3) = true;
                            dx(:,3) = [xtrj(i,j-1,k,fr) - xtrj(i,j,k,fr);...
                                       ytrj(i,j-1,k,fr) - ytrj(i,j,k,fr);...
                                       ztrj(i,j-1,k,fr) - ztrj(i,j,k,fr)];
                            dX(:,3) = [X(i,j-1,k) - X(i,j,k);...
                                       Y(i,j-1,k) - Y(i,j,k);...
                                       1e-10];
                        end

                        if (j+1 <= Isz(2)) && mask(i,j+1,k) && Nneighbor_SA(i,j,k) > 1
                            tf(4) = true;
                            dx(:,4) = [xtrj(i,j+1,k,fr) - xtrj(i,j,k,fr);...
                                       ytrj(i,j+1,k,fr) - ytrj(i,j,k,fr);...
                                       ztrj(i,j+1,k,fr) - ztrj(i,j,k,fr)];
                            dX(:,4) = [X(i,j+1,k) - X(i,j,k);...
                                       Y(i,j+1,k) - Y(i,j,k);...
                                       1e-10];
                        end

                        % Z direction
                        if (k-1 >= 1) && mask(i,j,k-1) && Nneighbor_LA(i,j,k) > 1
                            tf(5) = true;
                            dx(:,5) = [xtrj(i,j,k-1,fr) - xtrj(i,j,k,fr);...
                                       ytrj(i,j,k-1,fr) - ytrj(i,j,k,fr);...
                                       ztrj(i,j,k-1,fr) - ztrj(i,j,k,fr)];
                            dX(:,5) = [1e-10;...
                                       1e-10;...
                                       Z(i,j,k-1) - Z(i,j,k)];
                        end

                        if (k+1 <= Isz(3)) && mask(i,j,k+1) && Nneighbor_LA(i,j,k) > 1
                            tf(6) = true;
                            dx(:,6) = [xtrj(i,j,k+1,fr) - xtrj(i,j,k,fr);...
                                       ytrj(i,j,k+1,fr) - ytrj(i,j,k,fr);...
                                       ztrj(i,j,k+1,fr) - ztrj(i,j,k,fr)];
                            dX(:,6) = [1e-10;...
                                       1e-10;...
                                       Z(i,j,k+1) - Z(i,j,k)];
                        end
                        

                        % average deformation gradient tensor
                        Fave = dx(:,tf)/dX(:,tf);

                        % x/y strain tensor
                        E = 0.5*(Fave'*Fave - eye(3));

                        % coordinate system rotation matrix
                        % (Note this is the transpose of the vector rotation matrix)
                        Rot = [ct(i,j) st(i,j) 0; -st(i,j) ct(i,j) 0; 0 0 1];

                        % radial/circumferential strain tensor
                        Erot = Rot*E*Rot';

                        strain.RR(i,j,k,fr) = Erot(1); % 11
                        strain.RC(i,j,k,fr) = Erot(2); % 12
                        strain.RL(i,j,k,fr) = Erot(3); % 13
                        strain.CR(i,j,k,fr) = Erot(4); % 21
                        strain.CC(i,j,k,fr) = Erot(5); % 22
                        strain.CL(i,j,k,fr) = Erot(6); % 23
                        strain.LR(i,j,k,fr) = Erot(7); % 31
                        strain.LC(i,j,k,fr) = Erot(8); % 32
                        strain.LL(i,j,k,fr) = Erot(9); % 33
                    end

                end

            end

        end

    end

    %% FACE/VERTEX OUTPUT

    % unique faces within mask
    tf = mask & Nneighbor_SA>1 & Nneighbor_LA>1;
    strain.maskimage = tf;
    idx0 = find(tf);
    idx  = idx0(:,:,:,ones(Ntime,1));
    for fr = 1:Ntime
        idx(:,:,:,fr) = idx(:,:,:,fr) + Isz(3)*Isz(2)*Isz(1)*(fr-1);
    end

    tags = setdiff(fieldnames(strain),...
        {'orientation','maskimage'});
    for ti = 1:numel(tags)
        tag = tags{ti};
        strain.(tag) = reshape(strain.(tag)(idx),[numel(idx0),Ntime]);
    end

end
    
