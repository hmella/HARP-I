% Copyright (c) 2025 Hernan Mella
%
% TEMPORALPHASECONSISTENCY3D - Corrects temporal inconsistencies in 3D phase data.
%
% Description:
%   This function addresses temporal inconsistencies in 3D phase data obtained 
%   from imaging techniques, such as tagged MRI. It calculates optimal phase jumps 
%   for each frame to maintain temporal consistency, leveraging distance matrix 
%   optimization in 3D.
%
% Syntax:
%   phase_jumps = TemporalPhaseConsistency3D(Xpha, Ypha, Zpha, ref, njumps)
%
% Inputs:
%   Xpha     - [4D matrix] Unwrapped phase data along the X-axis.
%   Ypha     - [4D matrix] Unwrapped phase data along the Y-axis.
%   Zpha     - [4D matrix] Unwrapped phase data along the Z-axis.
%   ref      - [integer] Index of the reference frame.
%   njumps   - [integer, optional] Number of phase jump levels to test. Default: 1.
%
% Outputs:
%   phase_jumps - [3xN matrix] Optimal phase jumps for each frame (X, Y, Z).
%
% Example:
%   % Correct phase jumps for a 4D dataset
%   phase_jumps = TemporalPhaseConsistency3D(Xpha, Ypha, Zpha, 1, 2);
%
% Author:
%   Hernán Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% Licensing:
%   This Source Code Form is subject to the terms of the Mozilla Public License, 
%   version 2.0. If a copy of the MPL was not distributed with this file, You can 
%   obtain one at http://mozilla.org/MPL/2.0/.
% Notes:
%   - The approach leverages the relationship between distance matrices of
%     reference and current frames to estimate phase jumps, as detailed in the
%     PDF  by Mella, H., et al., "HARP-I: A Harmonic Phase Interpolation Method for
%     the Estimation of Motion From Tagged MR Images," IEEE Transactions on
%     Medical Imaging, vol. 40, no. 4, 2021. section "Strategies for Correcting
%     Temporal Phase Inconsistencies" 
%   - This method ensures phase consistency without imposing restrictions on
%     temporal resolution or segmentations 
%

function phase_jumps = TemporalPhaseConsistency3D(Xpha, Ypha, Zpha, ref, njumps)

        if nargin<5
            njumps = 1;
        end
    
        % Downsample phases to make the process faster
        Xpha = Downsample(Xpha,1:3,2);
        Ypha = Downsample(Ypha,1:3,2);
        Zpha = Downsample(Zpha,1:3,2);
        
        % Debugging
        show = false;
        
        % Number of frames
        Nfr = size(Xpha,4);
        
        % Reference phase
        ref_pha = [flatten(Xpha(:,:,:,ref)),flatten(Ypha(:,:,:,ref)),flatten(Zpha(:,:,:,ref))]';
        tf_ref = ~isnan(Xpha(:,:,:,ref)) & ~isnan(Ypha(:,:,:,ref)) & ~isnan(Zpha(:,:,:,ref));
        
        % Maximum value on reference phase
        ref_pha_max = max(ref_pha,[],2);
    
        % Reference distance matrix
        r_ref = rbfx.distanceMatrix3dNew(ref_pha(1,tf_ref),ref_pha(2,tf_ref),ref_pha(3,tf_ref),...
                                    ref_pha(1,tf_ref),ref_pha(2,tf_ref),ref_pha(3,tf_ref));
        r_ref_max = max(r_ref(:));
    
        % Output phase jumps
        phase_jumps = zeros([3 Nfr]);
    
        % Interpolation loop
        for i = 1:Nfr
            
            fprintf('\n Correcting phase inconsistencies on frame %d of %d',i,Nfr);        
    
            if i ~= 1
    
                % Phase further cardiac phases
                pha = [flatten(Xpha(:,:,:,i)),flatten(Ypha(:,:,:,i)),flatten(Zpha(:,:,:,i))]';
                tf = ~isnan(Xpha(:,:,:,i)) & ~isnan(Ypha(:,:,:,i)) & ~isnan(Zpha(:,:,:,i));
    
                % Estimate jumps to correct temporal inconsistencies
                pha_max = max(pha,[],2);
                delta = floor((ref_pha_max-pha_max)/(2*pi));
                jumps_x = [0, 2*pi*((delta(1)-njumps):(delta(1)+njumps))];
                jumps_y = [0, 2*pi*((delta(2)-njumps):(delta(2)+njumps))];    
                jumps_z = [0, 2*pi*((delta(3)-njumps):(delta(3)+njumps))];    
                init_diff = 1e+100;
                jump_x = 0;
                jump_y = 0;
                jump_z = 0;
                for jx=1:numel(jumps_x)
                    for jy=1:numel(jumps_y)
                        for jz=1:numel(jumps_z)
                            % Current phase
                            pha  = [flatten(Xpha(:,:,:,i))+jumps_x(jx),...
                                    flatten(Ypha(:,:,:,i))+jumps_y(jy),...
                                    flatten(Zpha(:,:,:,i))+jumps_z(jz)]';
    
                            % Distance matrix
                            r = rbfx.distanceMatrix3dNew(pha(1,tf),pha(2,tf),pha(3,tf),...
                            ref_pha(1,tf_ref),ref_pha(2,tf_ref),ref_pha(3,tf_ref));
                            r_max = max(r(:));
    
                            % Check optimality condition
                            if abs(r_max-r_ref_max) < init_diff                
                                init_diff = abs(r_max-r_ref_max);
                                jump_x = jumps_x(jx);
                                jump_y = jumps_y(jy);
                                jump_z = jumps_z(jz);
                            end
                        end
                    end
                end
    
                % Store optimum phase jumps
                phase_jumps(:,i) = [jump_x; jump_y; jump_z];
    
            end
    
        end
    
        % Debugging loop
        if show
          for i = 1:size(mask,4)
    
              % Current phase
              pha  = [flatten(Xpha(:,:,:,i)),flatten(Ypha(:,:,:,i)),flatten(Zpha(:,:,:,i))]';
              tf = mask(:,:,:,i);
    
              % Uncorrected distance matrix
              r = rbfx.Matrix3dNew(pha(1,tf),pha(2,tf),pha(3,tf),...
              ref_pha(1,tf_ref),ref_pha(2,tf_ref),ref_pha(3,tf_ref));
    
              % Corrected phase
              cpha  = [flatten(Xpha(:,:,:,i)),flatten(Ypha(:,:,:,i)),flatten(Zpha(:,:,:,i))]';
              cpha(1,:) = cpha(1,:) + phase_jumps(1,i);
              cpha(2,:) = cpha(2,:) + phase_jumps(2,i);
              cpha(3,:) = cpha(3,:) + phase_jumps(3,i);
    
              % Uncorrected distance matrix
              r = rbfx.distanceMatrix2dNew(cpha(1,tf),cpha(2,tf),...
              ref_pha(1,tf_ref),ref_pha(2,tf_ref));
    
              % plot images
              figure(1)
              t = tiledlayout(2,2);
              t.Padding = 'compact';
              t.TileSpacing = 'compact';
    
              nexttile
              imagesc(r_ref);
              caxis([min(r_ref(:)) max(r_ref(:))])
              box on; axis off;
              set(gca, 'Ydir', 'normal')            
    
              nexttile
              imagesc(r);
              caxis([min(r_ref(:)) max(r_ref(:))])            
              box on; axis off;
              set(gca, 'Ydir', 'normal')            
    
              nexttile
              imagesc(Xpha(:,:,1),'AlphaData',tf_ref); axis(ROI);
              caxis([min(ref_pha(1,tf_ref)) max(ref_pha(1,tf_ref))])
              box on; axis off;
              set(gca, 'Ydir', 'normal')            
    
              nexttile
              imagesc(Xpha(:,:,:,i)+phase_jumps(1,i),'AlphaData',tf); axis(ROI);
              caxis([min(ref_pha(1,tf_ref)) max(ref_pha(1,tf_ref))])
              box on; axis off;
              set(gca, 'Ydir', 'normal')  
    
              colormap(gray)
              pause
    
          end
      end
    
end
