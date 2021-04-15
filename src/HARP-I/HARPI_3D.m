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
            'RBFFactor',         [],...
            'RBFFacDist',        'quadratic',...
            'SpatialSmoothing',  1e-08,...
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
    Show = api.Show;                    % Show process?
    fr   = api.Frames;                  % frames to be analyzed
    Nfr  = numel(fr);                   % number of frames
    method = api.Method;                % interpolation method
    pspace = api.SpatialSmoothing;      % smoothing factor	
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
    s = api.RBFFactor;
    
    % TODO: add documentation for this part of the code
    if numel(s) > 1
        if strcmp(api.RBFFacDist,'DecreasingLinear')
            s = max(s) - (max(s)-min(s))*linspace(0,1,Nfr);
        elseif strcmp(api.RBFFacDist,'Linear')
            s = max(s) + (max(s)-min(s))*linspace(0,1,Nfr);
        elseif strcmp(api.RBFFacDist,'Quadratic')
            s = min(s) + abs(power(max(s)-min(s),0.4)*linspace(-1,1,Nfr)).^2.5;
        end
    else
        s = api.RBFFactor*ones([1 Nfr]);
    end
    
    % Grids
    [X, Y, Z] = meshgrid(1:Isz(2),1:Isz(1),1:Isz(3));
    dX = Downsample(X,1:3,dfac);
    dY = Downsample(Y,1:3,dfac);
    dZ = Downsample(Z,1:3,dfac);
        
    % Tissue positions on the deformed domain
    gy  = [dX(:),dY(:),dZ(:)]';
    
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
                [dXpha(:,:,:,i),dYpha(:,:,:,i),dZpha(:,:,:,i)] = RefPhaseSmoothing3D(dXpha(:,:,:,i),dYpha(:,:,:,i),dZpha(:,:,:,i),method,api.RBFFactor/5,pspace);
            end
            save([pwd,'/tmp_smoothed_unwrapped_phase.mat'],'dXpha','dYpha','dZpha')        
        end
        Pha_ref = [flatten(dXpha(:,:,:,1)),flatten(dYpha(:,:,:,1)),flatten(dZpha(:,:,:,1))]';
    end
    
    % Interpolation loop
    for i = 2:Nfr
    
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
    %         [xx_vec, Psi] = RBFInterp3D(Pha(:,tf),gy(:,tf),pspace,s(i),Pha_ref(:,tf_ref),phi);
    %         x_traj(tf_ref) = xx_vec(:,1);
    %         y_traj(tf_ref) = xx_vec(:,2);
    %         z_traj(tf_ref) = xx_vec(:,3);
    
            % RBF interpolation
            if j==1
                r  = rbfx.distanceMatrix3dNew(Pha(1,tf),Pha(2,tf),Pha(3,tf));
                B = phi.rbf(r,s(i));
                w = rbfx.solve(B,gy(:,tf)',pspace,true);
            end
            re = rbfx.distanceMatrix3dNew(Pha(1,tf),Pha(2,tf),Pha(3,tf),...
                         Pha_ref(1,tf_ref),Pha_ref(2,tf_ref),Pha_ref(3,tf_ref));
            Psi = phi.rbf(re,s(i));
            xx_vec = Psi*w;
            x_traj(tf_ref) = xx_vec(:,1);
            y_traj(tf_ref) = xx_vec(:,2);
            z_traj(tf_ref) = xx_vec(:,3);          
    
            % Update trajectories
            r = j:fac:Isz(1);
            traj(r,:,:,1,i) = x_traj;
            traj(r,:,:,2,i) = y_traj;
            traj(r,:,:,3,i) = z_traj;
            
    %         figure(1)
    %         sl=56;
    %         imagesc(squeeze(traj(:,:,sl,1,i)),'AlphaData',squeeze(~isnan(traj(:,:,sl,1,i))));
    %         axis equal;
    %         set(gca, 'Ydir', 'normal')
    %         drawnow
    %         pause
            
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
                % mag = sqrt((traj(:,:,sl,3,i)-Z(:,:,sl)).^2);
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
    % unwrapped_phases = permute(cat(5,dXpha,dYpha,dZpha),[1 2 3 5 4]);
    unwrapped_phases = [];
    metadata = struct(...
        'RawMotion',          u,...
        'FittedMotion',       u_fitted,...
        'RawTrajectories',    traj,...
        'FittedTrajectories', traj_fitted,...
        'UnwrappedPhases',    unwrapped_phases);
    
end