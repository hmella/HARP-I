% Copyright (c) 2025 Hernan Mella
%
% HARPI_3D - Harmonic Phase Interpolation for 3D Data
%
% Description:
%   This function performs harmonic phase interpolation on 3D image sequences, 
%   enabling phase unwrapping, displacement estimation, and optional temporal 
%   fitting. It is designed for motion estimation in 3D tagged MR images.
%
% Syntax:
%   metadata = HARPI_3D(Ih, 'ParameterName', ParameterValue, ...)
%
% Inputs:
%   Ih          - Input image sequence (5D array, size: [rows, cols, slices, channels, frames]).
%   varargin    - Optional parameters specified as name-value pairs:
%       'Mask'              - Binary mask (default: true(size(Ih))).
%       'FOV'               - Field of view (default: []).
%       'PixelSize'         - Pixel size (default: []).
%       'Frames'            - Frames to be analyzed (default: all frames).
%       'Show'              - Flag to display intermediate results (default: false).
%       'Method'            - Interpolation method (default: 'ThinPlate').
%       'a_constant'         - Radial basis function factor (default: []). Equation (16) and (17) on the paper.
%       'RBFFacDist'        - Distribution of RBF factors ('quadratic', default).
%       'eta_constant'  - Spatial smoothing factor (default: 1e-08). Equation (8) on the paper.
%       'Seed'              - Unwrapping seed (default: 'auto').
%       'Connectivity'      - Connectivity for unwrapping (default: 4).
%       'RefPhaseSmoothing' - Flag for reference phase smoothing (default: false).
%       'ROI'               - Region of interest (default: [1, rows, 1, cols, 1, slices]).
%       'TemporalFitting'   - Flag for temporal fitting (default: false).
%       'TemporalFittingOrder' - Temporal fitting order (default: 10).
%       'DownsamplingFac'   - Downsampling factor (default: 2).
%
% Outputs:
%   metadata - Struct containing:
%       RawMotion           - Raw displacement field.
%       FittedMotion        - Temporally fitted displacement field (if implemented).
%       RawTrajectories     - Raw tissue trajectories.
%       FittedTrajectories  - Temporally fitted tissue trajectories (if implemented).
%       UnwrappedPhases     - Unwrapped harmonic phase images (empty in 3D implementation).
%
% Example:
%   % Load 3D image sequence and mask
%   Ih = rand(128, 128, 64, 3, 10);  % Example 5D image sequence
%   mask = true(size(Ih, 1), size(Ih, 2), size(Ih, 3));
%
%   % Perform harmonic phase interpolation
%   metadata = HARPI_3D(Ih, 'Mask', mask, 'Frames', 1:10, 'Show', true);
%
%   % Access raw trajectories
%   raw_trajectories = metadata.RawTrajectories;
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopex (Benjamin.lopezf@usm.cl)
%
% License:
%   This Source Code Form is subject to the terms of the Mozilla Public License, 
%   version 2.0. If a copy of the MPL was not distributed with this file, 
%   you can obtain one at http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - This function extends the capabilities of HARPI to 3D datasets and is a 
%     **core algorithm** for motion estimation in 3D tagged MR images.
%   - Implements the techniques described in:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical 
%     Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%   - Note: Temporal fitting is currently a placeholder in this implementation 
%     and has not been fully developed for 3D data.
%
%
%   See also: parseinputs, Downsample, unwrap2, TemporalPhaseConsistency, RefPhaseSmoothing, RBFInterp2D, TemporalFitting

function [metadata] = HARPI_3D(Ih,varargin)
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
            'eta_constant',  1e-08,... % Equation (8) on the paper.
            'Seed',              'auto',...
            'Connectivity',      4,...
            'RefPhaseSmoothing', false,...
            'ROI',               [1 size(Ih,1) 1 size(Ih,2)],...
            'TemporalFitting',   false,...
            'TemporalFittingOrder', 10,...
            'DownsamplingFac',      2);
    
    % Check input
    api = parseinputs(defapi, [], varargin{:});
    
    % input parameters
    mask = api.Mask;                    % Mask
    Show = api.Show;                    % Show process
    fr   = api.Frames;                  % frames to be analyzed
    Nfr  = numel(fr);                   % number of frames
    method = api.Method;                % interpolation method
    eta = api.eta_constant;             % smoothing factor	
    seed = api.Seed;                    % unwrapping seed
    dfac = api.DownsamplingFac;
    
    % Resampling mask and images
    mask = mask(:,:,:,fr);
    fr = 1:Nfr;
    
    % Harmonic phases
    Xpha = squeeze(angle(Ih(:,:,:,1,:)));
    Ypha = squeeze(angle(Ih(:,:,:,2,:)));
    Zpha = squeeze(angle(Ih(:,:,:,3,:)));
    
    % Down-sample phase images to avoid large arrays
    dXpha = Downsample(Xpha,1:3,dfac);
    dYpha = Downsample(Ypha,1:3,dfac);
    dZpha = Downsample(Zpha,1:3,dfac);
    dmask = Downsample(mask,1:3,dfac);
    
    % Image sizes
    Isz = size(mask,1:3);
    dIsz = size(dmask,1:3);
    
    
    %% PHASE UNWRAPPING AND NORMALIZATION
    % Phase unwrapping
    try
        load([pwd,'/tmp_unwrapped_phase.mat'],'Xphar','Yphar','Zphar','dXpha','dYpha','dZpha')
    catch
        for i=1:Nfr
            fprintf('\n Unwrapping frame %d of %d',i,Nfr);
            dXpha(:,:,:,i) = unwrap3(dXpha(:,:,:,i),'Mask',dmask(:,:,:,i),'PixelSize',[1 1 1],...
                                'Seed',seed);
            dYpha(:,:,:,i) = unwrap3(dYpha(:,:,:,i),'Mask',dmask(:,:,:,i),'PixelSize',[1 1 1],...
                                'Seed',seed);
            dZpha(:,:,:,i) = unwrap3(dZpha(:,:,:,i),'Mask',dmask(:,:,:,i),'PixelSize',[1 1 1],...
                                'Seed',seed);
            if i==1
                fprintf('\n     Unwrapping reference frame');
                UV = NaN(Isz); UV(1:dfac:end,1:dfac:end,1:dfac:end) = dXpha(:,:,:,i);
                Xphar = unwrap3(Xpha(:,:,:,i),'Mask',mask(:,:,:,i),'PixelSize',[1 1 1],...
                                'Seed',seed,'UnwrappedVolume',UV);
                UV = NaN(Isz); UV(1:dfac:end,1:dfac:end,1:dfac:end) = dYpha(:,:,:,i);
                Yphar = unwrap3(Ypha(:,:,:,i),'Mask',mask(:,:,:,i),'PixelSize',[1 1 1],...
                                'Seed',seed,'UnwrappedVolume',UV);
                UV = NaN(Isz); UV(1:dfac:end,1:dfac:end,1:dfac:end) = dZpha(:,:,:,i);
                Zphar = unwrap3(Zpha(:,:,:,i),'Mask',mask(:,:,:,i),'PixelSize',[1 1 1],...
                                'Seed',seed,'UnwrappedVolume',UV);
            end                        
        end
        save([pwd,'/tmp_unwrapped_phase.mat'],'Xphar','Yphar','Zphar','dXpha','dYpha','dZpha')
    end
    
    if Show
        sl = round(dIsz/2);
        figure;
        for i=1:Nfr
            subplot 131;
            imagesc(dXpha(:,:,sl(3),i),'AlphaData',dmask(:,:,sl(3),i));
            axis equal;
            set(gca, 'Ydir', 'normal')
            subplot 132;
            imagesc(dYpha(:,:,sl(3),i),'AlphaData',dmask(:,:,sl(3),i));
            axis equal;
            subplot 133;
            imagesc(squeeze(dZpha(:,sl(2),:,i)),'AlphaData',squeeze(dmask(:,sl(2),:,i)));
            axis equal;
            set(gca, 'Ydir', 'normal')
            drawnow
        end
    end
    
    
    %% Displacement estimation
    % Get radial basis functions
    phi = feval(method);
    a = api.a_constant; % Equation (16) and (17) on the paper.
    
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
    [X, Y, Z] = meshgrid(1:Isz(2),1:Isz(1),1:Isz(3));
    dX = Downsample(X,1:3,dfac);
    dY = Downsample(Y,1:3,dfac);
    dZ = Downsample(Z,1:3,dfac);
        
    % Tissue positions on the deformed domain
    G  = [dX(:),dY(:),dZ(:)]';
    
    % Tissue trajectories
    tf_ref = ~isnan(Xphar) & ~isnan(Yphar) & ~isnan(Zphar);
    traj = NaN([Isz 3 Nfr]);
    traj(:,:,:,:,1) = cat(4,X,Y,Z);
    traj(~repmat(tf_ref,[1 1 1 3 Nfr])) = NaN;
    
    % Temporal phase consistency
    try
        load([pwd,'/tmp_phase_jumps.mat'],'phase_jumps')
    catch
        ref = 1;
        njumps = 2;
        phase_jumps = TemporalPhaseConsistency3D(dXpha,dYpha,dZpha,ref,njumps);
        save([pwd,'/tmp_phase_jumps.mat'],'phase_jumps');
    end
    dXpha = dXpha + reshape(phase_jumps(1,:),[1 1 1 numel(phase_jumps(1,:))]);
    dYpha = dYpha + reshape(phase_jumps(2,:),[1 1 1 numel(phase_jumps(2,:))]);
    dZpha = dZpha + reshape(phase_jumps(3,:),[1 1 1 numel(phase_jumps(3,:))]);
    
    % Reference phase
    if api.RefPhaseSmoothing
        try
            load([pwd,'/tmp_smoothed_unwrapped_phase.mat'],'dXpha','dYpha','dZpha')
        catch
            for i=1:Nfr
                fprintf('\n Smoothing phase on frame %d of %d',i,Nfr);
                [dXpha(:,:,:,i),dYpha(:,:,:,i),dZpha(:,:,:,i)] = RefPhaseSmoothing3D(dXpha(:,:,:,i),dYpha(:,:,:,i),dZpha(:,:,:,i),method,api.a_constant/5,eta);
            end
            save([pwd,'/tmp_smoothed_unwrapped_phase.mat'],'dXpha','dYpha','dZpha')        
        end
        Pha_ref = [flatten(dXpha(:,:,:,1)),flatten(dYpha(:,:,:,1)),flatten(dZpha(:,:,:,1))]';
    end
    
    % Interpolation loop
    for i = 1:Nfr
    
        fprintf('\n Estimating motion on frame %d of %d',i,Nfr);
    
        % Tissue phase on the deformed domain
        tf = ~isnan(dXpha(:,:,:,i)) & ~isnan(dYpha(:,:,:,i)) & ~isnan(dZpha(:,:,:,i));
        Pha  = [flatten(dXpha(:,:,:,i)),flatten(dYpha(:,:,:,i)),flatten(dZpha(:,:,:,i))]';    
        
        % Estimate motion on each frame in a iterative way to avoid
        % memory issues
        fac = 4;
        for j=1:fac
            
            % Downsample reference phases
            dXphar = Downsample(Xphar,1,fac,j-1);
            dYphar = Downsample(Yphar,1,fac,j-1);
            dZphar = Downsample(Zphar,1,fac,j-1);
    
            % Tissue position in the deformed domain mapped to the reference domain
            x_traj = NaN(size(dXphar)); % x(X)
            y_traj = NaN(size(dYphar)); % y(X)
            z_traj = NaN(size(dZphar)); % z(X)        
            
            % Tissue phase in the reference domain
            Pha_ref = [flatten(dXphar),flatten(dYphar),flatten(dZphar)]';
    
            % Reference tissue positions
            tf_ref = ~isnan(dXphar) & ~isnan(dYphar) & ~isnan(dZphar);
    
            % RBF interpolation
            if j==1
                r  = rbfx.distanceMatrix3dNew(Pha(1,tf),Pha(2,tf),Pha(3,tf));
                B = phi.rbf(r,a(i));
                w = rbfx.solve(B,G(:,tf)',eta,true);
            end
            re = rbfx.distanceMatrix3dNew(Pha(1,tf),Pha(2,tf),Pha(3,tf),...
                         Pha_ref(1,tf_ref),Pha_ref(2,tf_ref),Pha_ref(3,tf_ref));
            Psi = phi.rbf(re,a(i));
            xx_vec = Psi*w;
            x_traj(tf_ref) = xx_vec(:,1);
            y_traj(tf_ref) = xx_vec(:,2);
            z_traj(tf_ref) = xx_vec(:,3);          
    
            % Update trajectories
            r = j:fac:Isz(1);
            traj(r,:,:,1,i) = x_traj;
            traj(r,:,:,2,i) = y_traj;
            traj(r,:,:,3,i) = z_traj;
            
        end
    
        if api.Show
            fprintf('\nShowing processed data at frame %.0d\n',i)
            slices = [20 40 60 80];
            r = 1:dfac:Isz(1);
            figure(1),
            for sl=slices
                mag = sqrt((traj(:,:,sl,1,i)-X(:,:,sl)).^2 + ...
                           (traj(:,:,sl,2,i)-Y(:,:,sl)).^2 + ...
                           (traj(:,:,sl,3,i)-Z(:,:,sl)).^2);
                surf(X(:,:,sl),Y(:,:,sl),Z(:,:,sl),mag,'LineStyle','none'); hold on
                quiver3(X(r,r,sl),Y(r,r,sl),Z(r,r,sl),...
                        squeeze(traj(r,r,sl,1,i)-X(r,r,sl)),...
                        squeeze(traj(r,r,sl,2,i)-Y(r,r,sl)),...
                        squeeze(traj(r,r,sl,3,i)-Z(r,r,sl)),...
                        'AutoScale','off','Color','k'); hold on
            end
            hold off
            colorbar
            caxis([0 14])
            colormap(jet)
            axis off equal    
            view([0 10])
            drawnow
            
        end
        
    end
    
    % Output displacements
    u = traj - cat(4,X,Y,Z);
    
    % Temporal fitting
    u_fitted = [];
    traj_fitted = [];
    try
        if api.TemporalFitting
            args = struct(...
                'Mask',   mask,...
                'Frames', api.Frames,...
                'TemporalOrder', api.TemporalFittingOrder,...
                'Show', api.Show);
            [dx, dy] = TemporalFitting(u, args);
            u_fitted = permute(cat(5,dx,dy,dz),[1 2 3 5 4]);
            traj_fitted = cat(4,X,Y,Z) + u_fitted;
        end
    catch
        fprintf('\nTemporalFitting has not been yet implemented for 3D data!\n')
    end
    
    % Ouputs
    unwrapped_phases = [];
    metadata = struct(...
        'RawMotion',          u,...
        'FittedMotion',       u_fitted,...
        'RawTrajectories',    traj,...
        'FittedTrajectories', traj_fitted,...
        'UnwrappedPhases',    unwrapped_phases);
    
end