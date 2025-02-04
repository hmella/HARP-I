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
%   Mask           - (Optional) Binary mask defining active pixels to fit. Default: true(size(u)).
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
        'Mask',             true(size(u)),...
        'Frames',           [],...
        'TemporalOrder',    10,...
        'Show',             false);

    % parse inputs
    api = parseinputs(defapi,[],varargin{:});    
    
    % Mask
    mask = api.Mask;
    
    % Number of frames and frames to be resampled
    fr = api.Frames;
    Nfr = numel(fr);

    % Image size
    Isz = size(u, [1 2]);

    % Grids
    [X,Y] = meshgrid(1:Isz(1),1:Isz(2));    
    
    % Estimated displacements
    add_zero = false;
    if add_zero
        dx = zeros([Isz Nfr+1]); dx(:,:,2:end) = squeeze(u(:,:,1,:));
        dy = zeros([Isz Nfr+1]); dy(:,:,2:end) = squeeze(u(:,:,2,:));
    else
        dx = squeeze(u(:,:,1,:));
        dy = squeeze(u(:,:,2,:));
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

    % Equality constraints
    Aeq_fun = @(T,p) [T(1).^p;
                      p.*T(1).^max(0,p-1);
                      p.*T(end).^max(0,p-1)];
    Beq_fun = @(f,d1,d2) [f(1);
                          d1;
                          d2];
    % Polynomial fitting
    opt = optimoptions('lsqlin','Algorithm','interior-point','Display','off');
    for i=1:Isz(1)
        for j=1:Isz(2)
            if mask(i,j)
                % Fitting xtraj
                f = squeeze(dx(i,j,tmp));
                d1 = polyval(polyder(lsqlin(M4(1:4,:),f(1:4),[],[],[],[],[],[],[],opt)),time(1));
                d2 = polyval(polyder(lsqlin(M4(end-3:end,:),f(end-3:end),[],[],[],[],[],[],[],opt)),time(end));
                Aeq = Aeq_fun(time,p);
                Beq = Beq_fun(f,d1,d2);
                coef = lsqlin(M(~tf,:),f(~tf),[],[],Aeq,Beq,[],[],[],opt);
                if add_zero
                    dxr(i,j,:) = polyval(coef,time(2:end));
                else
                    dxr(i,j,:) = polyval(coef,time);
                end

                % Fitting ytraj
                f = squeeze(dy(i,j,tmp));
                d1 = polyval(polyder(lsqlin(M4(1:4,:),f(1:4),[],[],[],[],[],[],[],opt)),time(1));
                d2 = polyval(polyder(lsqlin(M4(end-3:end,:),f(end-3:end),[],[],[],[],[],[],[],opt)),time(end));
                Aeq = Aeq_fun(time,p);
                Beq = Beq_fun(f,d1,d2);
                coef = lsqlin(M(~tf,:),f(~tf),[],[],Aeq,Beq,[],[],[],opt);
                if add_zero
                    dyr(i,j,:) = polyval(coef,time(2:end));
                else
                    dyr(i,j,:) = polyval(coef,time);
                end
            end
        end
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