function [f, Rg, BfRg, WRg]= Butterworth23D(sze, WaveVec, n, sinmod)
% Design a butterworth bandpass filer for 2D/3D SinMod analysis.
% This is similar to a lowpass filter except a shift of origin
% The range of frequency matrix is from -0.5 to 0.5, x is horizontal
% direction, y in vertical direction, z in slice direction
% Usage: 
%         f = myBandFilter23D([256,256], [8,0], 5)
%         f = myBandFilter23D([256,256], [8,0,0], 5)
% Inputs:
%    sze    is a 2 or 3 element vector specifying the size of filter 
%           to construct [rows cols] for 2D. e.g, [256, 256]
%                        [rows cols slices] for 3D, e.g., [112,112,112]
%    WaveVec
%
%    n     : is the order of the filter, the higher n is the sharper
%            the transition is. (n must be an integer >= 1).
%            Note that n is doubled so that it is always an even integer.
% Outputs:
%                      1
%      f =    --------------------
%                              2n
%              1.0 + (w/cutoff)
%
% The frequency origin of the returned filter is at the corners.
% Refer to Peter Kovesi,
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/#freqfilt
%==========================================================================
% Example 1: 3D filter, wave vector in x direction, right-side spectrum ball
% sze = [112,112,112]; wave = [14 0 0];
% [f, Rg, BfRg, WRg]= myBandFilterButterworth23D(sze,wave, 5);
% BfRgM=zeros(sze);BfRgM(Rg)=BfRg;figure,imshow(BfRgM(:,:,sze(3)/2),[]);impixelinfo
% for i = 48:65 figure,imshow(squeeze(BfRgM(:,:,i)),[]);impixelinfo; end
%==========================================================================
% Example 2: 3D filter, wave vector in y direction, front-side spectrum ball
% sze = [112,112,112]; wave = [0 14 0];
% [f, Rg, BfRg, WRg]= myBandFilterButterworth23D(sze,wave, 5);
% BfRgM=zeros(sze);BfRgM(Rg)=BfRg;figure,imshow(BfRgM(:,:,sze(3)/2),[]);impixelinfo
% for i = 48:65 figure,imshow(squeeze(BfRgM(:,:,i)),[]);impixelinfo; end
%==========================================================================
% Example 3: 3D filter, wave vector in z direction, up-side spectrum ball
% sze = [112,112,112]; wave = [0 0 14];
% [f, Rg, BfRg, WRg]= myBandFilterButterworth23D(sze,wave, 5);
% BfRgM=zeros(sze);BfRgM(Rg)=BfRg;figure,imshow(squeeze(BfRgM(:,sze(2)/2,:)),[]);impixelinfo
% for i = 56:72 figure,imshow(squeeze(BfRgM(:,:,i)),[]);impixelinfo; end
%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hui Wang, 10/25/2011
% Update: 11/4/2011 
% Update: 2/28/2012, add examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% shift:  frequency shift of origin to off-center peak, instead of DC peak.
%         when shift = [0,0], it is a lowpass filter
% cutoff: the cutoff frequency of the filter 0 - 0.5
% WaveVec = Nx/d, shift/(1/2) = d/(Nx/2), therefore, shift = 1/WaveVec
% where d=(pixel distance between first and second peak)
dGrid=sqrt(sum(WaveVec.^2));  
shift = (1./dGrid)*(WaveVec/dGrid); 
cutoff = 0.8*(1./dGrid);    % at cutoff frequency, f = 0.5;
c=WaveVec(1); s=WaveVec(2); % [c,s] indicates wave direction and length

if cutoff < 0 | cutoff > 0.5
    error('cutoff frequency must be between 0 and 0.5');
end

if rem(n,1) ~= 0 | n < 1
    error('n must be an integer >= 1');
end

if length(sze) == 1
    rows = sze; cols = sze;
else
    rows = sze(1); cols = sze(2);
end

% Set up X and Y matrices with ranges normalised to +/- 0.5
% The following code adjusts things appropriately for odd and even values
% of rows and columns.
if mod(cols,2)
    xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
else
    xrange = [-cols/2:(cols/2-1)]/cols;	
end

if mod(rows,2)
    yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
else
    yrange = [-rows/2:(rows/2-1)]/rows;	
end

if numel(sze)==2
    [x,y] = meshgrid(xrange, yrange);
    radius = sqrt((x - shift(1)).^2 + (y - shift(2)).^2);   % A matrix with every pixel = radius relative to off centre.
    f = ifftshift( 1.0 ./ (1.0 + (radius ./ cutoff).^(2*n)) );
    %----Connect this with 2D SinMod, since SinMod needs Rg, BfRg, WRg
    Rg = find(f>0.005); % circle in polar frequency domain, is 0.005 ok?
    BfRg = f(Rg);
    if sinmod
        WRg=c*x(Rg)+s*y(Rg); % frequency component along wave vector
    else
        WRg=c*x+s*y; % frequency component along wave vector
    end
    %----------------------------------------------------------------

elseif numel(sze)==3
    slices = sze(3);
    v=WaveVec(3);  % vertical direction, slice direction
    if mod(slices,2)
        zrange = [-(slices-1)/2:(slices-1)/2]/(slices-1);
    else
        zrange = [-slices/2:(slices/2-1)]/slices;	
    end
    
    [x,y,z] = meshgrid(xrange, yrange, zrange);
    x = x+1e-8; y = y+1e-8; z = z+1e-8;    % Avoid 0
    radius = sqrt((x - shift(1)).^2 + (y - shift(2)).^2 + (z - shift(3)).^2);   % A matrix with every pixel = radius relative to off centre.
    f =  1.0./(1.0 + (radius./cutoff).^(2*n));

    %----Connect this with 3D SinMod, since SinMod needs Rg, BfRg, WRg----
    Rg = find(f > 0.05); 
    BfRg = f(Rg);
    if sinmod
        WRg=c*x(Rg)+s*y(Rg)+v*z(Rg); % frequency component along wave vector
    else
        WRg=c*x+s*y+v*z; % frequency component along wave vector
    end        
    %BfRgM=zeros(sze);BfRgM(Rg)=BfRg;figure,imshow(BfRgM(:,:,sze(3)/2),[]); % debug, hui
    %---------------------------------------------------------------------
end
end

