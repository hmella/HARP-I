function phase_jumps = TemporalPhaseConsistency(Xpha,Ypha,ref,njumps,ROI)

    if nargin<5
        njumps = 1;
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

    % Reference distance matrix
    r_ref = distanceMatrix2dNew(ref_pha(1,tf_ref),ref_pha(2,tf_ref),...
                ref_pha(1,tf_ref),ref_pha(2,tf_ref));
    r_ref_max = max(r_ref(:));

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
            init_diff = 1e+100;
            jump_x = 0;
            jump_y = 0;
            for jx=1:numel(jumps_x)
                for jy=1:numel(jumps_y)
                    % Current phase
                    pha  = [flatten(Xpha(:,:,i))+jumps_x(jx),...
                            flatten(Ypha(:,:,i))+jumps_y(jy)]';

                    % Distance matrix
                    r = distanceMatrix2dNew(pha(1,tf),pha(2,tf),...
                    ref_pha(1,tf_ref),ref_pha(2,tf_ref));
                    r_max = max(r(:));

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
      for i = 1:Nfr

          % Current phase
          pha  = [flatten(Xpha(:,:,i)),flatten(Ypha(:,:,i))]';
          tf = ~isnan(Xpha(:,:,i)) & ~isnan(Ypha(:,:,i));

          % Uncorrected distance matrix
          r = distanceMatrix2dNew(pha(1,tf),pha(2,tf),...
          ref_pha(1,tf_ref),ref_pha(2,tf_ref));

          % Corrected phase
          cpha  = [flatten(Xpha(:,:,i)),flatten(Ypha(:,:,i))]';
          cpha(1,:) = cpha(1,:) + phase_jumps(1,i);
          cpha(2,:) = cpha(2,:) + phase_jumps(2,i);

          % Uncorrected distance matrix
          r = distanceMatrix2dNew(cpha(1,tf),cpha(2,tf),...
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