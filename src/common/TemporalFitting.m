function [dxr,dyr] = TemporalFitting(u,varargin)

    %% TEMPORAL FITTING
    % we use the function lsqlin and polyval of the optimization
    % toolbox to determine the fitted coefficents at all spatial 
    % locations within the tissue segmentation

    %% PARSE INPUT
    % Default arguments
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
    % Aeq_fun = @(T,p) [T(1).^p;
    %                   T(end).^p;
    %                   p.*T(1).^max(0,p-1);
    %                   p.*T(end).^max(0,p-1)];
    % Beq_fun = @(f,d1,d2) [f(1);
    %                       f(end);
    %                       d1;
    %                       d2];
    Aeq_fun = @(T,p) [T(1).^p;
                      p.*T(1).^max(0,p-1);
                      p.*T(end).^max(0,p-1)];
    Beq_fun = @(f,d1,d2) [f(1);
                          d1;
                          d2];
    % Aeq_fun = @(T,p) [p.*T(1).^max(0,p-1);
    %                   p.*T(end).^max(0,p-1)];
    % Beq_fun = @(f,d1,d2) [d1;
    %                       d2];
    
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

                % if add_zero
                %     figure(1)
                %     subplot 121
                %     plot(squeeze(dx(i,j,2:end))); hold on
                %     plot(squeeze(dxr(i,j,:))); hold off
                %     legend('Original','Fitted')
                %     subplot 122
                %     plot(squeeze(dy(i,j,2:end))); hold on
                %     plot(squeeze(dyr(i,j,:))); hold off
                %     legend('Original','Fitted')                    
                % else
                %     figure(1)
                %     subplot 121
                %     plot(time,squeeze(dx(i,j,:))); hold on
                %     plot(time,squeeze(dxr(i,j,:))); hold off
                %     legend('Original','Fitted')
                %     subplot 122
                %     plot(time,squeeze(dy(i,j,:))); hold on
                %     plot(time,squeeze(dyr(i,j,:))); hold off
                %     legend('Original','Fitted')
                % end
                
            end
        end
    end


    %% PLots
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