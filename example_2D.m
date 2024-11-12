clear; clc; close all;

%% Load image
load('data/data.mat')

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
% filter = HARPFilter(struct('Image',I,'CentralFreq',ke,'Direction',direction,...
%                            'FilterType','Transmission'));

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
        'RBFFactor',            [1 1]*150/(s/pxsz(1)),...
        'RBFFacDist',           'Linear',...
        'SpatialSmoothing',     15,...
        'Connectivity',         4,...
        'RefPhaseSmoothing',    true,...
        'ROI',                  [86 132 68 114],...
        'TemporalFitting',      true,...
        'TemporalFittingOrder', 10);
harpi = HARPI(If, args);
t = toc;
dxr = squeeze(harpi.FittedMotion(:,:,1,:));
dyr = squeeze(harpi.FittedMotion(:,:,2,:));
fprintf('\n Elapsed time HARPI: %d (s)', t)

% Strain
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
RR = NaN([Isz(1) Isz(2) Nfr]);
CC = NaN([Isz(1) Isz(2) Nfr]);
RR(repmat(st.maskimage,[1 1 Nfr])) = st.p1(:);
CC(repmat(st.maskimage,[1 1 Nfr])) = st.p2(:);

% Show estimated strain maps
figure(2)
tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
nexttile
imagesc(CC(:,:,8),'AlphaData',st.maskimage); set(gca,'YDir','Normal');
colormap(jet); axis off; colorbar; caxis([-0.5 0.2])
title('CC')
nexttile
imagesc(RR(:,:,8),'AlphaData',st.maskimage); set(gca,'YDir','Normal'); axis off
colormap(flipud(jet)); axis off; colorbar; caxis([-0.2 0.5])
title('RR')


%% Strain by segments
% Number of segments
Nseg = 6;

% HARPI
api = struct(...
    'CC',               CC,...
    'RR',               RR,...
    'Mask',             st.maskimage,...
    'Nseg',             Nseg,...
    'ClockWise',        false,...
    'Frames',           1:Nfr);
harpi_seg = getStrainBySegments(api);

% Plot strain by segments
figure,
errorbar(1:Nfr,mean(harpi_seg.segments_CC),...
         std(harpi_seg.segments_CC),'-','LineWidth',2); hold on
errorbar(1:Nfr,mean(harpi_seg.segments_RR),...
         std(harpi_seg.segments_RR),'-','LineWidth',2); hold on

%% Show results
% Plot results
close all;

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
    imagesc(CC(:,:,i),'AlphaData',~isnan(CC(:,:,i))); caxis([a b])
    set(gca,'YDir','normal') 
    cc = flipud(jet);
    colormap(ax2,cc)
    axis([50 120 50 120])
    daspect([1 1 1])
    pause(0.1)
end
