% Copyright (c) 2025 Hernan Mella
%
% TEMPORALPHASECONSISTENCY - Corrects temporal inconsistencies in phase data.
%
% Description:
%   This function addresses temporal inconsistencies in phase data obtained from
%   imaging techniques like tagged Magnetic Resonance Imaging (MRI). The method
%   calculates phase jumps between frames and adjusts them to maintain temporal
%   consistency, leveraging a distance matrix optimization approach.
%
% Syntax:
%   phase_jumps = TemporalPhaseConsistency(Xpha, Ypha, ref, njumps, ROI)
%
% Inputs:
%   Xpha     - [3D matrix] Unwrapped phase data along the X-axis.
%   Ypha     - [3D matrix] Unwrapped phase data along the Y-axis.
%   ref      - [integer] Index of the reference frame.
%   njumps   - [integer, optional] Number of phase jump levels to test. Default: 1.
%   ROI      - [array, optional] Region of interest for debugging visualization.
%
% Outputs:
%   phase_jumps - [2xN matrix] Optimal phase jumps for each frame (X and Y).
%
% Example:
%   % Correct phase jumps for a 3D MRI dataset
%   phase_jumps = TemporalPhaseConsistency(Xpha, Ypha, 1, 2);
%
% Author:
%   Hernán Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% License:
%   This Source Code Form is subject to the terms of the Mozilla Public License, 
%   version 2.0. If a copy of the MPL was not distributed with this file, 
%   you can obtain one at http://mozilla.org/MPL/2.0/.
%
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

function phase_jumps = TemporalPhaseConsistency(Xpha,Ypha,ref,njumps,ROI)

    % njumps is the fourth positional argument and ROI the fifth.
    if nargin < 4 || isempty(njumps)
        njumps = 1;
    end
    if nargin < 5
        ROI = []; %#ok<NASGU>  (used only by the debug plots below)
    end

    % Downsample phases to make the process faster
    Xpha = Downsample(Xpha,1:2,2);
    Ypha = Downsample(Ypha,1:2,2);

    % Debugging
    show = false;

    % Number of frames
    Nfr = size(Xpha,3);

    % Reference phase
    ref_pha = [flatten(Xpha(:,:,ref)),flatten(Ypha(:,:,ref))]';
    tf_ref = ~isnan(Xpha(:,:,ref)) & ~isnan(Ypha(:,:,ref));
    
    % Maximum value on reference phase
    ref_pha_max = max(ref_pha,[],2);

    % Largest distance within the reference point set. Only the maximum of
    % the distance matrix is used, and it is reached at a convex hull
    % vertex, so the full N-by-N matrix is not needed (see extremePoints).
    ref_hull  = extremePoints(ref_pha(:,tf_ref));
    r_ref_max = maxPairDistance(ref_hull, ref_hull);

    % Output phase jumps
    phase_jumps = zeros([2 Nfr]);

    % Interpolation loop
    for i = 1:Nfr

        if i ~= ref

            % Phase further cardiac phases
            pha = [flatten(Xpha(:,:,i)),flatten(Ypha(:,:,i))]';
            tf = ~isnan(Xpha(:,:,i)) & ~isnan(Ypha(:,:,i));

            % Estimate jumps to correct temporal inconsistencies
            pha_max = max(pha,[],2);
            delta = floor((ref_pha_max-pha_max)/(2*pi));
            jumps_x = [0, 2*pi*((delta(1)-njumps):(delta(1)+njumps))];
            jumps_y = [0, 2*pi*((delta(2)-njumps):(delta(2)+njumps))];
            init_diff = Inf;
            jump_x = 0;
            jump_y = 0;

            % A candidate jump only translates this frame's point set, so
            % the hull is computed once per frame rather than once per
            % candidate.
            pha_hull = extremePoints(pha(:,tf));

            for jx=1:numel(jumps_x)
                for jy=1:numel(jumps_y)
                    % Shifted hull of the current phase
                    shifted = [pha_hull(1,:) + jumps_x(jx); ...
                               pha_hull(2,:) + jumps_y(jy)];

                    r_max = maxPairDistance(shifted, ref_hull);

                    % Check optimality condition
                    if abs(r_max-r_ref_max) < init_diff
                        init_diff = abs(r_max-r_ref_max);
                        jump_x = jumps_x(jx);
                        jump_y = jumps_y(jy);
                    end

                end
            end

            % Store optimum phase jumps
            phase_jumps(:,i) = [jump_x; jump_y];

        end

    end

    % Debugging loop
    if show
      r_ref = rbfx.distanceMatrix2dNew(ref_pha(1,tf_ref),ref_pha(2,tf_ref),...
                  ref_pha(1,tf_ref),ref_pha(2,tf_ref));
      for i = 1:Nfr

          % Current phase
          pha  = [flatten(Xpha(:,:,i)),flatten(Ypha(:,:,i))]';
          tf = ~isnan(Xpha(:,:,i)) & ~isnan(Ypha(:,:,i));

          % Uncorrected distance matrix
          r = rbfx.distanceMatrix2dNew(pha(1,tf),pha(2,tf),...
          ref_pha(1,tf_ref),ref_pha(2,tf_ref));

          % Corrected phase
          cpha  = [flatten(Xpha(:,:,i)),flatten(Ypha(:,:,i))]';
          cpha(1,:) = cpha(1,:) + phase_jumps(1,i);
          cpha(2,:) = cpha(2,:) + phase_jumps(2,i);

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
          box on; axis off equal;
          set(gca, 'Ydir', 'normal')            

          nexttile
          imagesc(r);
          caxis([min(r_ref(:)) max(r_ref(:))])            
          box on; axis off equal;
          set(gca, 'Ydir', 'normal')            

          nexttile
          imagesc(Xpha(:,:,1),'AlphaData',tf_ref);
          caxis([min(ref_pha(1,tf_ref)) max(ref_pha(1,tf_ref))])
          box on; axis off equal;
          set(gca, 'Ydir', 'normal')  
          axis(ROI);

          nexttile
          imagesc(Xpha(:,:,i) + phase_jumps(1,i),'AlphaData',tf);
          caxis([min(ref_pha(1,tf_ref)) max(ref_pha(1,tf_ref))])
          box on; axis off equal;
          set(gca, 'Ydir', 'normal')
          axis(ROI);

          colormap(jet)
          sgtitle(sprintf('frame %02d',i))
          pause

      end
  end

end
