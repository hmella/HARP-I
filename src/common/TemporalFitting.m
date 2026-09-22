% Copyright (c) 2025 Hernan Mella
%
% TEMPORALFITTING - Performs temporal fitting of displacement fields using polynomial interpolation.
%
% Description:
%   This function fits displacement fields (u) over time using polynomial interpolation. It ensures temporal smoothness 
%   by fitting a polynomial trajectory to the x and y displacement components for each pixel in a 2D grid. Temporal fitting 
%   is particularly useful for resampling displacement data with constraints on derivatives at boundary frames.
%
% Syntax:
%   [dxr, dyr] = TemporalFitting(u, 'Mask', mask, 'Frames', fr, 'TemporalOrder', order, 'Show', show)
%
% Inputs:
%   u              - Displacement field (size: [rows, cols, components, frames]).
%                    `u(:,:,1,:)` represents x-displacements.
%                    `u(:,:,2,:)` represents y-displacements.
%   Mask           - (Optional) Binary mask defining active pixels to fit. Default: true(rows, cols).
%   Frames         - (Optional) Indices of frames to resample. Default: [] (all frames).
%   TemporalOrder  - (Optional) Temporal order of the polynomial fit. Default: 10.
%   Show           - (Optional) Boolean flag to visualize the results. Default: false.
%
% Outputs:
%   dxr - Fitted/resampled x-displacements (size: [rows, cols, frames]).
%   dyr - Fitted/resampled y-displacements (size: [rows, cols, frames]).
%
% Example:
%   % Displacement field input (example size: 128x128 grid, 2 components, 30 frames)
%   u = randn(128, 128, 2, 30);
%   mask = true(128, 128);
%   frames = 1:30;
%   [dxr, dyr] = TemporalFitting(u, 'Mask', mask, 'Frames', frames, 'TemporalOrder', 5, 'Show', true);
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (benjamin.lopezf@usm.cl)
%
% License:
%   This Source Code Form is subject to the terms of the Mozilla Public 
%   License, v. 2.0. If a copy of the MPL was not distributed with this 
%   file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - The polynomial fitting ensures smooth trajectories while preserving temporal consistency.
%   - Boundary conditions are enforced using first and second derivatives at the initial and final frames.
%   - Refer to Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for
%   the Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical Imaging, 
%   vol. 40, no. 4, pp. 1240-1251, April 2021, for methods inspiring this functionality.
%   - This based on references 48,50 and 51 thats appear in HARP-I
%

function [dxr,dyr] = TemporalFitting(u,varargin)

    defapi = struct(...
        'Mask',             true([size(u,1) size(u,2)]),...
        'Frames',           [],...
        'TemporalOrder',    10,...
        'Show',             false);

    % parse inputs
    [api, extra_args] = parseinputs(defapi,[],varargin{:});
    assertNoUnknownArgs(mfilename, extra_args);    
    
    % Mask
    mask = logical(api.Mask);
    
    % Number of frames and frames to be resampled
    fr = api.Frames;
    Nfr = numel(fr);

    % Image size
    Isz = size(u, [1 2]);

    % Grids
    [X,Y] = meshgrid(1:Isz(1),1:Isz(2));    
    
    % Estimated displacements
    add_zero = false;
    % Reshape rather than squeeze: squeeze would also collapse a singleton
    % image dimension (a 1-by-1 ROI, or a single row/column), silently
    % turning the [rows cols frames] stack into something else.
    Nt = size(u,4);
    if add_zero
        dx = zeros([Isz Nt+1]); dx(:,:,2:end) = reshape(u(:,:,1,:), [Isz Nt]);
        dy = zeros([Isz Nt+1]); dy(:,:,2:end) = reshape(u(:,:,2,:), [Isz Nt]);
    else
        dx = reshape(u(:,:,1,:), [Isz Nt]);
        dy = reshape(u(:,:,2,:), [Isz Nt]);
    end

    % Resampled displacements
    dxr = NaN([Isz Nfr]);
    dyr = NaN([Isz Nfr]);
    
    % temporal resampling range
    if add_zero
        time = 0:Nfr;
        tmp = [true fr];
    else
        time = 1:Nfr;
        tmp = fr;
    end
    time = (time-time(1))/(time(end)-time(1));
    tf = time == 0;
    
    % Vandermonde matrices
    p = api.TemporalOrder:-1:0;
    M  = bsxfun(@power,time(:),p);
    M4 = bsxfun(@power,time(:),3:-1:0);

    % Equality constraints: value at t=0 and slope at both ends.
    Aeq = [time(1).^p;
           p.*time(1).^max(0,p-1);
           p.*time(end).^max(0,p-1)];

    % ------------------------------------------------------------------
    % Batched, equality-constrained polynomial fit.
    %
    % This used to call lsqlin six times per masked pixel. Two of those
    % solved a square 4x4 system (an exactly determined problem handed to an
    % interior-point solver), and the loop-invariant Aeq was rebuilt twice
    % per pixel. The problem is
    %
    %     min ||A c - f||^2   s.t.   Aeq c = beq
    %
    % with A = M(~tf,:) and Aeq both independent of the pixel. Solve it by
    % the null-space method: with Z spanning null(Aeq), any feasible c is
    % c0 + Z*y, and y falls out of an unconstrained least-squares problem in
    % the fixed matrix A*Z. Everything that does not depend on the pixel --
    % Z, pinv(Aeq) and the QR factorisation of A*Z -- is computed once, so
    % each pixel costs one triangular solve and all pixels go through as a
    % single matrix of right-hand sides.
    %
    % The null-space route is used rather than the KKT normal equations
    % because forming A'A would square the condition number of an already
    % ill-conditioned order-10 Vandermonde matrix.
    % ------------------------------------------------------------------

    if Nfr < 4
        error([mfilename ':tooFewFrames'], ...
            'Temporal fitting needs at least 4 frames; got %d.', Nfr);
    end
    if api.TemporalOrder + 1 > Nfr
        error([mfilename ':orderTooHigh'], ...
            'TemporalOrder %d needs at least %d frames; got %d.', ...
            api.TemporalOrder, api.TemporalOrder+1, Nfr);
    end

    A  = M(~tf,:);
    Z  = null(Aeq);            % nc x (nc-3), constant
    P  = pinv(Aeq);            % nc x 3,      constant
    AZ = A*Z;

    if ~isempty(AZ)
        c = cond(AZ);
        if ~isfinite(c) || c > 1e12
            warning([mfilename ':illConditioned'], ...
                ['Temporal fit is ill-conditioned (cond = %.2e); consider a ' ...
                 'TemporalOrder below %d for %d frames.'], c, api.TemporalOrder, Nfr);
        end
    end
    dAZ = decomposition(AZ, 'qr');

    % The historical code indexed a 3-D mask with two subscripts, which
    % selects the first frame. Keep that meaning explicitly.
    if ~ismatrix(mask)
        mask = mask(:,:,1);
    end

    idx = find(mask);
    if ~isempty(idx)
        FX = reshape(dx, prod(Isz), []);  FX = FX(idx, tmp).';   % Nfr x Npix
        FY = reshape(dy, prod(Isz), []);  FY = FY(idx, tmp).';

        M4a = M4(1:4,:);
        M4b = M4(end-3:end,:);

        DXR = fit_component(FX, A, M, M4a, M4b, time, tf, Z, P, dAZ);
        DYR = fit_component(FY, A, M, M4a, M4b, time, tf, Z, P, dAZ);

        dxr = scatter_back(DXR, idx, Isz, Nfr);
        dyr = scatter_back(DYR, idx, Isz, Nfr);
    end

    % PLots
    if api.Show
        figure,
        for i=1:3:Isz(1)
            for j=1:3:Isz(2)
                if mask(i,j)
                    plot(X(i,j), Y(i,j), 'sk', 'MarkerFaceColor', 'k'); hold on
                    plot(squeeze(X(i,j) + dxr(i,j,:)),squeeze(Y(i,j) + dyr(i,j,:)),'Color','k','LineWidth',2); hold on
                    plot(squeeze(X(i,j) + dx(i,j,:)),squeeze(Y(i,j) + dy(i,j,:)),'-.','Color',[0,0,0]+0.5,'LineWidth',2); hold on
                end
            end
        end
        legend('Origins','Fitted','Original')
        hold off
        pause
    end

end

% =========================================================================
function VALS = fit_component(F, A, M, M4a, M4b, time, tf, Z, P, dAZ)
%FIT_COMPONENT Constrained polynomial fit for every pixel of one component.
%   F is Nfr x Npix; returns Nfr x Npix of fitted values.

    % End slopes from local cubic fits. These are square systems, so the
    % backslash gives the exact solution lsqlin was approximating.
    C4a = M4a \ F(1:4,:);
    C4b = M4b \ F(end-3:end,:);

    t1 = time(1);
    t2 = time(end);
    d1 = 3*C4a(1,:)*t1^2 + 2*C4a(2,:)*t1 + C4a(3,:);
    d2 = 3*C4b(1,:)*t2^2 + 2*C4b(2,:)*t2 + C4b(3,:);

    % Constraint values: fix f(0) and both end slopes
    BEQ = [F(1,:); d1; d2];

    C0 = P*BEQ;                       % minimum-norm feasible point
    Y  = dAZ \ (F(~tf,:) - A*C0);     % least squares within the null space
    C  = C0 + Z*Y;                    % coefficients, highest power first

    VALS = M*C;                       % == polyval(c, time) for every column
end

% =========================================================================
function out = scatter_back(VALS, idx, Isz, Nfr)
%SCATTER_BACK Place per-pixel time courses back on the image grid.
    out  = NaN([Isz Nfr]);
    flat = reshape(out, [], Nfr);
    flat(idx,:) = VALS.';
    out  = reshape(flat, [Isz Nfr]);
end
