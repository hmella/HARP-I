% Copyright (c) 2025 Hernan Mella
%
% ITOK - Convert data from image space to k-space (frequency domain).
%
% Syntax:
%   Data = itok(Data)
%   Data = itok(Data, Dims)
%
% Description:
%   This function applies the Fourier transform to convert data from the
%   spatial domain (image space) to the frequency domain (k-space).
%   The transform can be applied across all dimensions or specific ones.
%
% Inputs:
%   Data - Multidimensional array containing the input data in image space.
%   Dims - (Optional) Dimensions along which the Fourier transform should be
%          applied. If not provided, the transform is applied across all dimensions.
%
% Outputs:
%   Data - Data converted to k-space, with the same size as the input.
%
% Examples:
%   % Transform data from image space to k-space (all dimensions)
%   Data_img = rand(64, 64); % Simulated image space data
%   Data_k = itok(Data_img);
%
%   % Transform data along specific dimensions (e.g., first dimension)
%   Data_k = itok(Data_img, [1]);
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
%   - Inspired by applications in MRI, where data is processed in k-space 
%     before being reconstructed to the image domain.
%   - Referenced: Mella et al., "HARP-I: A Harmonic Phase Interpolation Method
%     for the Estimation of Motion From Tagged MR Images," IEEE Transactions 
%     on Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%

function [Data] = itok(Data, Dims)
if nargin<2,
  Data = fftshift(fftn(ifftshift(Data)));
else
  for d=[1:length(Dims)],
    Data = fftshift(fft(ifftshift(Data),[],Dims(d)));
  end
end
return

