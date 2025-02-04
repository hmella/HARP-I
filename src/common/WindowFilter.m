% Copyright (c) 2025 Hernan Mella
%
% WINDOWFILTER - A class to create and apply customizable window filters.
%
% This class allows the creation of filters based on Riesz and Tukey
% windows. These filters can be applied to images or other data
% for smoothing or feature enhancement.
%
% Properties:
%   weights - The filter weights calculated based on the chosen window type.
%   size    - Size of the filter (default: [50 50]).
%   width   - Relative width of the window transition region (default: 0.6).
%   lift    - Scaling factor to adjust the filter's baseline (default: 0.3).
%   type    - Type of filter ('Riesz' or 'Tukey', default: 'Riesz').
%
% Methods:
%   WindowFilter(size, width, lift, type) - Constructor to create the filter.
%   Riesz(size, width, lift)             - Generates a Riesz window.
%   Tukey(size, width, lift)             - Generates a Tukey window.
%   filter(I)                            - Applies the filter to input data.
%
% Examples:
%   Create a Riesz filter
%   W = WindowFilter([100, 100], 0.5, 0.2, 'Riesz');
%
%   Create a Tukey filter
%   T = WindowFilter([100, 100], 0.4, 0.1, 'Tukey');
%
%   % Apply the filter to an image
%   filteredImage = W.filter(inputImage);
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
% - Concepts based on methods described in Mella et al., "HARP-I: A Harmonic Phase Interpolation
%   Method for the Estimation of Motion From Tagged MR Images," IEEE Transactions on Medical Imaging,
%   vol. 40, no. 4, pp. 1240-1251, April 2021.
% - Reference: DOI 10.1109/TMI.2021.3051092
%

classdef WindowFilter

  properties
    weights
    size = [50 50];
    width = 0.6;
    lift = 0.3;
    type = 'Riesz';
  end
  methods
    % Constructor
    function obj = WindowFilter(size, width, lift, type)
        obj.size = size;
        obj.width = width;
        obj.lift = lift;
        obj.type = type;

        % Create filter
        if strcmp(obj.type,'Riesz')
            obj.weights = obj.Riesz(obj.size, obj.width, obj.lift);
        elseif strcmp(obj.type,'Tukey')
            obj.weights = obj.Tukey(obj.size, obj.width, obj.lift);
        end
    end

    % Riesz filter
    function H0 = Riesz(obj, Size, width, lift)
      decay = (1.0-width)/2;
      s = Size;
      s20 = round(decay*s);
      s1 = linspace(0, s20/2, s20);
      w1 = 1.0 - power(abs(s1/(s20/2)),2).*(1.0-lift);
  
      % Set up filter
      H0 = ones([1,s]);
      H0(1:s20) = H0(1:s20).*flip(w1);
      H0((s-s20+1):s) = H0((s-s20+1):s).*w1;
    end

    % Tukey filter
    function H0 = Tukey(obj, Size, width, lift)
      alpha = 1.0 - width;
      H0 = tukeywin(Size, alpha)*(1.0-lift) + lift;
      H0 = H0';
    end

    % Filter image
    function Ih = filter(obj,I)
        H  = repmat(obj.weights, [1 1 1 obj.image_size(end)]);
        Ih = ktoi(H.*itok(I, [1 2]), [1 2]);
    end

  end
end