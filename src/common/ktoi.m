% Copyright (c) 2025 Hernan Mella
%
% KTOI - Convert data from k-space (frequency domain) to image space.
%
% Syntax:
%   Data = ktoi(Data)
%   Data = ktoi(Data, Dims)
%
% Description:
%   This function applies the inverse Fourier transform to convert data
%   from the frequency domain (k-space) to the spatial domain (image space).
%   The transform can be applied across all dimensions or specific ones.
%
% Inputs:
%   Data - Multidimensional array containing the input data in k-space.
%   Dims - (Optional) Dimensions along which the inverse Fourier transform 
%          should be applied. If not provided, the transform is applied
%          across all dimensions.
%
% Outputs:
%   Data - Data converted to image space, with the same size as the input.
%
% Examples:
%   % Transform data from k-space to image space (all dimensions)
%   Data_k = rand(64, 64); % Simulated k-space data
%   Data_img = ktoi(Data_k);
%
%   % Transform data along specific dimensions (e.g., first dimension)
%   Data_img = ktoi(Data_k, [1]);
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
%   - Inspired by applications in MRI, where data is acquired in k-space and
%     must be transformed to the image domain for visualization and analysis.
%   - Referenced: Mella et al., "HARP-I: A Harmonic Phase Interpolation Method
%     for the Estimation of Motion From Tagged MR Images," IEEE Transactions 
%     on Medical Imaging, vol. 40, no. 4, pp. 1240-1251, April 2021.
%   - Reference: DOI 10.1109/TMI.2021.3051092
%

function [Data] = ktoi(Data, Dims)
if nargin<2,
  Data = fftshift(ifftn(ifftshift(Data)));
else
  for d=[1:length(Dims)],
    Data = fftshift(ifft(ifftshift(Data),[],Dims(d)));
  end
end
return
