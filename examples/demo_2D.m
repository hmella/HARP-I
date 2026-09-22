clear; clc; close all;

%% Load image
load('../data/data_2D.mat')

% Image
I = data.image;

% Image size and number of frames
Isz = size(I);
Nfr = Isz(4);


%% IMAGE FILTERING
% Pixel size, tag spacing, and encoding frequency
pxsz = data.pixelsize/1000;
s    = data.tagspacing/1000;
ke   = 2*pi/s*pxsz;

% Create filter
direction = deg2rad([0 90]);
filter = HARPFilter(struct('Image',I,'CentralFreq',ke,'Direction',direction,...
                           'FilterType','Transmission','Butterworth_cuttoff',17,...
                           'Butterworth_order',10));

% Get harmonic image
If = filter.filter(I);

figure(1)
tiledlayout(2,2,'Padding','compact','TileSpacing','compact')
nexttile
imagesc(abs(itok(I(:,:,1,1)))); set(gca,'YDir','Normal'); axis off
nexttile
imagesc(abs(itok(If(:,:,1,1)))); set(gca,'YDir','Normal'); axis off
nexttile
imagesc(abs(itok(I(:,:,2,1)))); set(gca,'YDir','Normal'); axis off
nexttile
imagesc(abs(itok(If(:,:,2,1)))); set(gca,'YDir','Normal'); axis off
sgtitle('Filtered k-spaces')


%% HARPI displacement and strain
% Segmentation
mask = data.mask;

% Motion estimation
tic
args = struct(...
        'Mask',                 mask,...
        'Frames',               1:Nfr,...
        'Show',                 false,...
        'Method',               'Multiquadric3',...
        'a_constant',           [1 1]*150/(s/pxsz(1)),... % equation (16) and (17) on the paper.
        'RBFFacDist',           'Linear',...
        'eta_constant',         15,... % equation (8) on the paper.
        'Connectivity',         4,...
        'RefPhaseSmoothing',    true,...
        'TemporalFitting',      true,...
        'TemporalFittingOrder', 10);
harpi = HARPI(If, args);
t = toc;
dxr = squeeze(harpi.FittedMotion(:,:,1,:));
dyr = squeeze(harpi.FittedMotion(:,:,2,:));
fprintf('\n Elapsed time HARPI: %d (s)', t)

%% Strain
[X, Y] = meshgrid(1:size(dxr,2), 1:size(dxr,1));
options = struct(...
    'X', X,...
    'Y', Y,...
    'mask',mask(:,:,1),...
    'times',1:Nfr,...
    'dx', dxr,...
    'dy', dyr,...
    'Origin', [],...
    'Orientation', []);
st = pixelstrain(options);

% Which strain is which:
%   st.RR / st.CC  radial and circumferential components of the Lagrangian
%                  strain tensor, as reported in the paper and plotted by
%                  strain_curves2D.m
%   st.p1 / st.p2  1st and 2nd principal strains, the eigenvalues of E
% All four are computed below. The switch selects what the animation shows.
STRAIN_FOR_DISPLAY = 'p2';   % 'CC' circumferential (paper) | 'p2' 2nd principal

expand = @(v) reshape_masked(v, st.maskimage, Isz, Nfr);
RR = expand(st.RR);   % radial (Lagrangian)
CC = expand(st.CC);   % circumferential (Lagrangian)
P1 = expand(st.p1);   % 1st principal
P2 = expand(st.p2);   % 2nd principal

switch STRAIN_FOR_DISPLAY
    case 'CC', StrainMap = CC; strain_label = 'E_{cc} (circumferential)';
    case 'p2', StrainMap = P2; strain_label = '2nd principal strain';
    otherwise, error('Unknown STRAIN_FOR_DISPLAY: %s', STRAIN_FOR_DISPLAY);
end

%% Save motion and strain estimation
save('motion2D_HARPI.mat','harpi','st','RR','CC','P1','P2')

%% Show results
% gif path
gifFileName = 'StrainMap.gif';

% displacements
ui = permute(cat(4,dxr,dyr),[1 2 4 3]);

% colorbar axis
a = -0.5;
b = 0.25;

% quiver subsampling
f = 2;

% quiver data
m = st.maskimage;
x = X(1:f:end,1:f:end);
y = Y(1:f:end,1:f:end);

%
figure,
for i=1:size(ui,4)

    % Create tiled layout
    t = tiledlayout(1,2);
    t.TileSpacing = 'compact';

    % displacements
    uxi = ui(1:f:end,1:f:end,1,i);
    uyi = ui(1:f:end,1:f:end,2,i);

    % Quiver data
    ax1 = nexttile;
    imagesc(abs(I(:,:,1,i).*I(:,:,2,i))); hold on
    quiver(x,y,uxi,uyi,...
           'Color','r','AutoScale','off'); hold off
    set(gca,'YDir','normal') 
    colormap(ax1,'gray')
    axis([50 120 50 120])
    daspect([1 1 1])

    % Strain map
    ax2 = nexttile;
    imagesc(StrainMap(:,:,i),'AlphaData',~isnan(StrainMap(:,:,i))); clim([a b])
    title(strain_label)
    set(gca,'YDir','normal') 
    cc = flipud(jet);
    colormap(ax2,cc)
    axis([50 120 50 120])
    daspect([1 1 1])

    frame = getframe(gcf);
    img = frame2im(frame); 
    [A, map] = rgb2ind(img, 256);

    if i == 1
        imwrite(A, map, gifFileName, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
    else
        imwrite(A, map, gifFileName, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end

    pause(0.1);
end

function out = reshape_masked(vals, maskimage, Isz, Nfr)
%RESHAPE_MASKED Scatter per-pixel strain values back onto the image grid.
    out = NaN([Isz(1) Isz(2) Nfr]);
    out(repmat(maskimage,[1 1 Nfr])) = vals(:);
end
