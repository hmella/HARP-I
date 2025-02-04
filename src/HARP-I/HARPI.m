% Copyright (c) 2025 Hernan Mella
%
% HARPI - Harmonic Phase Interpolation
%
% Description:
%   This function performs harmonic phase interpolation on a given input
%   image sequence. It includes phase unwrapping, displacement estimation, 
%   and optional temporal fitting, and is particularly designed for motion 
%   estimation in tagged MR images.
%
% Syntax:
%   metadata = HARPI(Ih, 'ParameterName', ParameterValue, ...)
%
% Inputs:
%   Ih          - Input image sequence (4D array, size: [rows, cols, channels, frames]).
%   varargin    - Optional parameters specified as name-value pairs:
%       'Mask'              - Binary mask (default: true(size(Ih))).
%       'FOV'               - Field of view (default: []).
%       'PixelSize'         - Pixel size (default: []).
%       'Frames'            - Frames to be analyzed (default: all frames).
%       'Show'              - Flag to display intermediate results (default: false).
%       'Method'            - Interpolation method (default: 'ThinPlate').
%       'a_constant'         - Radial basis function factor (default: []). Equation (16) and (17) on the paper.
%       'RBFFacDist'        - Distribution of RBF factors ('quadratic', default).
%       'eta_constant'      - Spatial smoothing factor (default: 1e-08). Equation (8) on the paper.
%       'Seed'              - Unwrapping seed (default: 'auto').
%       'Connectivity'      - Connectivity for unwrapping (default: 4).
%       'RefPhaseSmoothing' - Flag for reference phase smoothing (default: false).
%       'ROI'               - Region of interest (default: [1, rows, 1, cols]).
%       'TemporalFitting'   - Flag for temporal fitting (default: false).
%       'TemporalFittingOrder' - Temporal fitting order (default: 10).
%       'DownsamplingFac'   - Downsampling factor (default: 1).
%
% Outputs:
%   metadata - Struct containing:
%       RawMotion           - Raw displacement field.
%       FittedMotion        - Temporally fitted displacement field.
%       RawTrajectories     - Raw tissue trajectories.
%       FittedTrajectories  - Temporally fitted tissue trajectories.
%       UnwrappedPhases     - Unwrapped harmonic phase images.
%
% Example:
%   % Load image sequence and mask
%   Ih = rand(256, 256, 2, 10);  % Example 4D image sequence
%   mask = true(size(Ih, 1), size(Ih, 2));
%
%   % Perform harmonic phase interpolation
%   metadata = HARPI(Ih, 'Mask', mask, 'Frames', 1:10, 'Show', true);
%
%   % Access fitted trajectories
%   fitted_trajectories = metadata.FittedTrajectories;
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (Benjamin.lopezf@usm.cl)
%
% License:
%   This Source Code Form is subject to the terms of the Mozilla Public License, 
%   version 2.0. If a copy of the MPL was not distributed with this file, 
%   you can obtain one at http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - This is the **main function** of the repository and encapsulates the 
%     core topic addressed in the accompanying paper.
%   - The function implements the methods described in:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical 
%     Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - This function contains the primary algorithm described in the paper 
%     and is essential for the estimation of motion in tagged MR images.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%
%   See also: parseinputs, Downsample, unwrap2, TemporalPhaseConsistency, 
%             RefPhaseSmoothing, RBFInterp, TemporalFitting

function [metadata] = HARPI(Ih,varargin)

    %HARPI Summary of this function goes here
    %   Detailed explanation goes here
  
    %% PARSE INPUT
    % Default arguments
    defapi = struct(...
            'Mask',              true(size(Ih)),...
            'FOV',               [],...
            'PixelSize',         [],...
            'Frames',            [],...
            'Show',              false,...
            'Method',            'ThinPlate',...
            'a_constant',         [],... % Equation (16) and (17) on the paper.
            'RBFFacDist',        'quadratic',...
            'eta_constant',  1e-08,...
            'UnwrappingSeed',    'auto',...
            'Connectivity',      4,...
            'RefPhaseSmoothing', false,...
            'ROI',               [1 size(Ih,1) 1 size(Ih,2)],...
            'TemporalFitting',   false,...
            'TemporalFittingOrder', 10,...
            'DownsamplingFac',      1);
  
    % Check input
    api = parseinputs(defapi, [], varargin{:});
  
    % Image size
    Isz = size(Ih(:,:,1,1));
  
    % input parameters
    mask = api.Mask;                    % Mask
    fr   = api.Frames;                  % frames to be analyzed
    Nfr  = numel(fr);                   % number of frames
    method = api.Method;                % interpolation method
    eta = api.eta_constant;             % smoothing factor	
    seed = api.UnwrappingSeed;          % unwrapping seed
    dfac = api.DownsamplingFac;
  
    % Resampling mask and images
    mask = mask(:,:,fr);
    fr = 1:Nfr;
  
    % Harmonic phases
    Xpha = squeeze(angle(Ih(:,:,1,:)));
    Ypha = squeeze(angle(Ih(:,:,2,:)));
  
    % Down-sample phase images to avoid large arrays
    dXpha = Downsample(Xpha,1:2,dfac);
    dYpha = Downsample(Ypha,1:2,dfac);
    dmask = Downsample(mask,1:2,dfac);
  
    if api.Show
        figure;
        for i=1:Nfr
            subplot 121
            imagesc(dXpha(:,:,i),'AlphaData',dmask(:,:,i))
            set(gca, 'Ydir', 'normal')
            colormap(gray)
            subplot 122
            imagesc(dYpha(:,:,i),'AlphaData',dmask(:,:,i))
            set(gca, 'Ydir', 'normal')
            colormap(gray)
            pause(0.1)
        end
    end
  
    [X,Y] = meshgrid(1:Isz(2),1:Isz(1));
  
  
    %% PHASE UNWRAPPING AND NORMALIZATION
    % Phase unwrapping loop
    for i=1:Nfr
        % Unwrap phases
        dXpha(:,:,i) = unwrap2(dXpha(:,:,i),'Mask',dmask(:,:,i),'PixelSize',[1 1],...
                            'Seed',seed,'Connectivity',api.Connectivity);
        dYpha(:,:,i) = unwrap2(dYpha(:,:,i),'Mask',dmask(:,:,i),'PixelSize',[1 1],...
                            'Seed',seed,'Connectivity',api.Connectivity);
  
        % Unwrap fully sampled reference if needed
        if and(i==1,dfac~=1)
            UV = NaN(Isz); UV(1:dfac:end,1:dfac:end) = dXpha(:,:,i);
            Xphar = unwrap2(Xpha(:,:,i),'Mask',mask(:,:,i),'PixelSize',[1 1],...
                                  'Seed',seed,'Connectivity',api.Connectivity,...
                                  'UnwrappedImage',UV);
            UV(1:dfac:end,1:dfac:end) = dYpha(:,:,i);
            Yphar = unwrap2(Ypha(:,:,i),'Mask',mask(:,:,i),'PixelSize',[1 1],...
                                'Seed',seed,'Connectivity',api.Connectivity,...
                                'UnwrappedImage',UV);
        elseif and(i==1,dfac==1)
            Xphar = dXpha(:,:,1);
            Yphar = dYpha(:,:,1);
        end
        
    end
  
    % Show unwrapped phases
    if api.Show
        figure;
        for i=1:Nfr
            subplot 121;
            imagesc(dXpha(:,:,i),'AlphaData',dmask(:,:,i));
            colormap(jet);
            axis equal;
            set(gca, 'Ydir', 'normal')
            subplot 122;
            imagesc(dYpha(:,:,i),'AlphaData',dmask(:,:,i));
            colormap(jet);
            axis equal;
            set(gca, 'Ydir', 'normal')
            pause(0.1)
        end
    end
  
  
    %% Displacement estimation
    % Get radial basis functions
    phi = feval(method);
    a = api.a_constant;
  
    % TODO: add documentation for this part of the code
    if numel(a) > 1
        if strcmp(api.RBFFacDist,'DecreasingLinear')
            a = max(a) - (max(a)-min(a))*linspace(0,1,Nfr);
        elseif strcmp(api.RBFFacDist,'Linear')
            a = max(a) + (max(a)-min(a))*linspace(0,1,Nfr);
        elseif strcmp(api.RBFFacDist,'Quadratic')
            a = min(a) + abs(power(max(a)-min(a),0.4)*linspace(-1,1,Nfr)).^2.5;
        end
    else
        a = api.a_constant*ones([1 Nfr]);
    end
  
    % Grids
    [dX, dY] = meshgrid(1:dfac:Isz(2),1:dfac:Isz(1));
    [X, Y] = meshgrid(1:Isz(2),1:Isz(1));
  
    % Tissue phase in the reference domain
    Pha_ref = [flatten(Xphar),flatten(Yphar)]';
    tf_ref = ~isnan(Xphar) & ~isnan(Yphar);
  
    % Tissue positions on the deformed domain
    G  = [dX(:),dY(:)]';
  
    % Tissue trajectories
    traj = NaN([Isz 2 Nfr]);
    traj(:,:,:,1) = cat(3,X,Y);
    traj(~repmat(tf_ref,[1 1 2 Nfr])) = NaN;
  
    % Tissue position in the deformed domain mapped to the reference domain
    x_traj = NaN(Isz); % x(X)
    y_traj = NaN(Isz); % y(X)
  
    % Temporal phase consistency
    ref    = 1;
    njumps = 8;
    phase_jumps = TemporalPhaseConsistency(dXpha,dYpha,ref,njumps,api.ROI);
    dXpha = dXpha + reshape(phase_jumps(1,:),[1 1 numel(phase_jumps(1,:))]);
    dYpha = dYpha + reshape(phase_jumps(2,:),[1 1 numel(phase_jumps(2,:))]);
  
    % Reference phase
    if api.RefPhaseSmoothing
        smooth_f = a(i)/5; % This parameter was chosen empirically
        [Xphar,Yphar] = RefPhaseSmoothing(Xphar,Yphar,method,smooth_f,eta);    
        Pha_ref = [flatten(Xphar),flatten(Yphar)]';
    end
  
    % Interpolation loop
    for i = 2:Nfr
  
        % Tissue phase on the deformed domain
        tf = ~isnan(dXpha(:,:,i)) & ~isnan(dYpha(:,:,i));
        Pha  = [flatten(dXpha(:,:,i)),flatten(dYpha(:,:,i))]';
  
        % RBF interpolation
        [xx_vec, Psi] = rbfx.RBFInterp2D(Pha(:,tf),G(:,tf),eta,a(i),Pha_ref(:,tf_ref),phi);
        x_traj(tf_ref) = xx_vec(:,1);
        y_traj(tf_ref) = xx_vec(:,2);
        
        % Update trajectories
        traj(:,:,1,i) = x_traj;
        traj(:,:,2,i) = y_traj;
  
        if api.Show
            fprintf('\nShowing processed data at frame %.0d\n',i)
            ss = 5;
            A=Ih(:,:,1,i);
            B=Ih(:,:,2,i);
  
            subplot 121
            imagesc(X(:), Y(:), abs(A).*abs(B)); 
            caxis([0 0.4*max(abs(A).*abs(B),[],'all')]); colormap gray; hold on
            q = quiver(X(1:ss:end,1:ss:end), Y(1:ss:end,1:ss:end),...
                traj(1:ss:end,1:ss:end,1,i)-X(1:ss:end,1:ss:end),...
                traj(1:ss:end,1:ss:end,2,i)-Y(1:ss:end,1:ss:end)); hold off
            q.AlignVertexCenters = 'on';
            q.AutoScale = 'off';
            q.Color = 'r';
            q.Marker = 'none';
            q.MarkerEdgeColor = 'r';
            q.MarkerFaceColor = 'r';
            q.MarkerSize = 2;
            q.MaxHeadSize = 0.1;
            axis off equal
            set(gca, 'Ydir', 'normal')
  
            subplot 122
            imagesc(Psi); colorbar; colormap jet
            set(gca, 'Ydir', 'normal')
  
            pause(0.05)
        end
    end
  
    % Output displacements
    u = traj - cat(3,X,Y);
  
    % Temporal fitting
    u_fitted = u;
    traj_fitted = traj;
    if api.TemporalFitting
        args = struct(...
            'Mask',   mask,...
            'Frames', api.Frames,...
            'TemporalOrder', api.TemporalFittingOrder,...
            'Show', api.Show);
        [dx, dy] = TemporalFitting(u, args);
        u_fitted = permute(cat(4,dx,dy),[1 2 4 3]);
        traj_fitted = cat(3,X,Y) + u_fitted;
    end
  
    % Ouputs
    unwrapped_phases = permute(cat(4,dXpha,dYpha),[1 2 4 3]);
    metadata = struct(...
        'RawMotion',          u,...
        'FittedMotion',       u_fitted,...
        'RawTrajectories',    traj,...
        'FittedTrajectories', traj_fitted,...
        'UnwrappedPhases',    unwrapped_phases);
  
    end
  