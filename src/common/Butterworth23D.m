% Copyright (c) 2025 Hernan Mella
%
% BUTTERWORTH23D - Creates a 2D/3D Butterworth bandpass filter for SinMod analysis.
%
% Description:
%   This function designs a Butterworth bandpass filter for 2D or 3D data.
%   The filter is used for frequency-domain analysis, especially in applications
%   requiring frequency separation along specified wave vectors.
%
% Syntax:
%   [f, Rg, BfRg, WRg] = Butterworth23D(sze, WaveVec, n, sinmod)
%
% Inputs:
%   sze    - A vector specifying the size of the filter:
%            [rows, cols] for 2D (e.g., [256, 256])
%            [rows, cols, slices] for 3D (e.g., [112, 112, 112]).
%   WaveVec - Wave vector defining the direction and magnitude of the filter
%             [WaveX, WaveY] for 2D, [WaveX, WaveY, WaveZ] for 3D.
%   n      - Order of the Butterworth filter (integer >= 1). Higher values 
%            result in sharper transitions.
%   sinmod - Boolean flag to determine whether to calculate sinusoidal 
%            modulation components along the wave vector.
%
% Outputs:
%   f      - The Butterworth filter in frequency space.
%   Rg     - Indices of significant frequency components (circle in the polar
%            frequency domain).
%   BfRg   - Filter values at the significant frequency components.
%   WRg    - Frequency components along the wave vector.
%
% Example 1: 2D Filter
%   sze = [256, 256]; 
%   wave = [8, 0]; 
%   [f, Rg, BfRg, WRg] = Butterworth23D(sze, wave, 5, false);
%
% Example 2: 3D Filter
%   sze = [112, 112, 112]; 
%   wave = [0, 14, 0]; 
%   [f, Rg, BfRg, WRg] = Butterworth23D(sze, wave, 5, true);
%
% Author:
%   Hui Wang (https://scholar.google.com/citations?user=Ym03m88AAAAJ&hl=en)
%   Updated: November 4, 2011, February 28, 2012
%
% License:
%   This Source Code Form is subject to the terms of the Mozilla Public License, 
%   version 2.0. If a copy of the MPL was not distributed with this file, 
%   you can obtain one at http://mozilla.org/MPL/2.0/.
%
% Notes:
%   - This implementation aligns with methods described in:
%     Mella et al., "HARP-I: A Harmonic Phase Interpolation Method for the 
%     Estimation of Motion From Tagged MR Images," IEEE Transactions on 
%     Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%


function [f, Rg, BfRg, WRg]= Butterworth23D(sze, WaveVec, n, sinmod)
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
    Rg = find(f>0.005); % circle in polar frequency domain, is 0.005 ok?
    BfRg = f(Rg);
    if sinmod
        WRg=c*x(Rg)+s*y(Rg); % frequency component along wave vector
    else
        WRg=c*x+s*y; % frequency component along wave vector
    end

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
    Rg = find(f > 0.05); 
    BfRg = f(Rg);
    if sinmod
        WRg=c*x(Rg)+s*y(Rg)+v*z(Rg); % frequency component along wave vector
    else
        WRg=c*x+s*y+v*z; % frequency component along wave vector
    end        
end
end

