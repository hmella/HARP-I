clear; clc; close all

% Load data for processing
path = '../data/data_3D.mat';
load(path,'I','mask','info')

% Image size
Isz = size(I,1:3);
Nenc = size(I,4);
Nfr = size(I,5);

% Permute dimensions to achieve taglines ordered in (x,y,z) directions
Itmp = I;
I(:,:,:,1,:) = Itmp(:,:,:,3,:);
I(:,:,:,2,:) = Itmp(:,:,:,1,:);
I(:,:,:,3,:) = Itmp(:,:,:,2,:);
clear Itmp

% Plot LA tagged image to see actual basal displacement
figure,
imagesc(abs(squeeze(I(:,:,80,1,2)).*squeeze(I(:,:,80,2,2)))); set(gca,'YDir','Normal')
axis off equal
box on
colormap gray

figure,
imagesc(permute(abs(squeeze(I(50,:,:,1,2)).*squeeze(I(50,:,:,3,2))),[2 1])); set(gca,'YDir','Normal')
axis off equal
box on
colormap gray
print('-dpng','-r600','LA_tagged_1')

figure,
imagesc(permute(abs(squeeze(I(56,:,:,1,10)).*squeeze(I(56,:,:,3,10))),[2 1])); set(gca,'YDir','Normal')
axis off equal
box on
colormap gray
print('-dpng','-r600','LA_tagged_2')


%% IMAGE FILTERING
% Encoding frequency
spacing = info.TagSpacing;
vxsz    = [info.PixelSpacing' info.PixelSpacing(1)];
ke      = 2*pi./spacing;

% Get filter
filter_type = {'Butterworth','Transmission','Gabor'};
wave_vecs = 2*(spacing/vxsz(1))*[1,0,0;
                                 0,1,0;
                                 0,0,1];
filter = HARPFilter3D(struct('Image',I,'WaveVec',wave_vecs,'SinMod',false,...
                            'FilterType','Butterworth','Gabor_lambda',2,...
                            'Butterworth_order',5));
            
% Get harmonic image
Ih = filter.filter(I);

% Kspace
K = itok(I,[1 2 3]);
Kh = itok(Ih,[1 2 3]);

% Show filtered images
frame = 1;
figure(1)
tiledlayout(3,2,'Padding','compact','TileSpacing','compact')
nexttile
imagesc(abs(K(:,:,56,1,frame))); set(gca,'YDir','Normal'); caxis([0 1e+8])
nexttile
imagesc(abs(Kh(:,:,56,1,frame))); set(gca,'YDir','Normal'); caxis([0 1e+8])
nexttile
imagesc(squeeze(abs(K(:,56,:,3,frame)))); set(gca,'YDir','Normal')
nexttile
imagesc(squeeze(abs(Kh(:,56,:,3,frame)))); set(gca,'YDir','Normal'); caxis([0 1e+8])
nexttile
imagesc(angle(Ih(:,:,56,1,frame))); set(gca,'YDir','Normal')
nexttile
imagesc(squeeze(angle(Ih(:,56,:,2,frame)))); set(gca,'YDir','Normal')
drawnow

% Remove some variables to avoid memory issues
clear I K Kh H filter Itmp


%% HARPI displacement and strain
tic
args = struct(...
        'Mask',             mask,...
        'Frames',           1:Nfr,...
        'Show',             false,...
        'Method',           'Multiquadric3',...
        'a_constant',        [1 1]*150/(spacing/vxsz(1)),... % Equation (16) and (17) on the paper.
        'RBFFacDist',       'Linear',...
        'eta_constant',     15,... % Equation (8) on the paper.
        'Connectivity',     4,...
        'RefPhaseSmoothing',false,...
        'ROI',              [50 120 50 120],...
        'TemporalFitting',  false,...
        'TemporalFittingOrder', 10,...
        'DownsamplingFac',   4);
harpi = HARPI_3D(Ih, args);
t = toc;
fprintf('\n Elapsed time HARPI: %d (s)', t)
mag = sqrt(squeeze(harpi.RawMotion(:,:,:,1,:)).^2 + ...
           squeeze(harpi.RawMotion(:,:,:,2,:)).^2 + ...
           squeeze(harpi.RawMotion(:,:,:,3,:)).^2);

         

%% STRAIN ESTIMATION
[X,Y,Z] = meshgrid(1:Isz(2),1:Isz(1),1:Isz(3));
options = struct(...
    'X', X,...
    'Y', Y,...
    'Z', Z,...
    'mask',mask(:,:,:,1),...
    'times',1:Nfr,...
    'dx', squeeze(harpi.RawMotion(:,:,:,1,:)),...
    'dy', squeeze(harpi.RawMotion(:,:,:,2,:)),...
    'dz', squeeze(harpi.RawMotion(:,:,:,3,:)),...
    'Origin', [],...
    'Orientation', []);
st = pixelstrain3D(options);
RR = NaN([Isz Nfr]);
CC = NaN([Isz Nfr]);
LL = NaN([Isz Nfr]);
RR(repmat(st.maskimage,[1 1 1 Nfr])) = st.RR(:);
CC(repmat(st.maskimage,[1 1 1 Nfr])) = st.CC(:);
LL(repmat(st.maskimage,[1 1 1 Nfr])) = st.LL(:);


%% Save motion and strain estimation
save('motion3D_HARPI.mat','harpi','st')

%% SHOW RESULTS
% Downsampling factor
dfac = 4;
slices = [20 40 60 80];
r = 1:dfac:Isz(1);
filename = 'slices_animation.gif'; 

for i = 1:Nfr
    figure('Visible', 'off'); 
    
    for s = slices
        surf(X(:,:,s), Y(:,:,s), Z(:,:,s), squeeze(LL(:,:,s,i)), 'LineStyle', 'none'); hold on
        quiver3(X(r,r,s), Y(r,r,s), Z(r,r,s),...
                squeeze(harpi.RawMotion(r,r,s,1,i)),...
                squeeze(harpi.RawMotion(r,r,s,2,i)),...
                squeeze(harpi.RawMotion(r,r,s,3,i)),...
                'AutoScale', 'off', 'Color', 'k'); hold on
    end
    
    hold off
    caxis([-0.5 0.2])
    colormap(flipud(jet))
    axis off equal    
    view([0 10])
    drawnow
    
    frame = getframe(gcf);
    img = frame2im(frame);
    [A, map] = rgb2ind(img, 256);
    
    if i == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

close all force

