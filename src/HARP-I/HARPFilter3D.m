% Copyright (c) 2025 Hernan Mella
%
% HARPFilter3D - Implements a 3D HARP filter for volumetric image analysis.
%
% Description:
%   This class is an extension of the 2D HARPFilter, specifically designed
%   for 3D image processing and analysis. It allows the application of 3D 
%   filters (e.g., Butterworth) in the frequency domain (k-space) to extract
%   and analyze harmonic phase information, making it suitable for tasks
%   like motion estimation, strain analysis, and segmentation in volumetric
%   datasets such as 3D tagged MRI.
%
% Properties:
%   - image_size: Dimensions of the input 3D image.
%   - wave_vec: Wave vectors representing spatial frequencies of interest.
%   - sinmod: Boolean flag to enable/disable SinMod processing (default: false).
%   - type: Type of filter ('Butterworth').
%   - lambda: Wavelength parameter for Gabor filter (reserved for future use).
%   - cutoff: Cutoff frequency for the Butterworth filter.
%   - order: Order of the Butterworth filter.
%   - kspace_filter: Generated k-space filters for the input data.
%   - wave_vectors: Wave vectors used to define central frequencies.
%   - Rg, BfRg, WRg: Parameters for advanced SinMod processing.
%
% Methods:
%   - Constructor: Initializes the HARPFilter3D object with specified configurations.
%   - filter: Applies the configured filter to the input 3D image.
%
% Usage:
%   % Example: Apply a 3D HARP filter to a volumetric image
%   input_image = rand(64, 64, 64, 1, 10); % 5D image data
%   wave_vectors = [1 0 0; 0 1 0; 0 0 1]; % Define wave vectors for 3D directions
%
%   harp_filter = HARPFilter3D('Image', input_image, ...
%                              'WaveVec', wave_vectors, ...
%                              'FilterType', 'Butterworth', ...
%                              'Butterworth_cuttoff', 30, ...
%                              'Butterworth_order', 5);
%
%   filtered_image = harp_filter.filter(input_image);
%
% Author:
%   Hernan Mella (hernan.mella@pucv.cl)
%
% Collaborator:
%   Benjamin Lopez (Benjamin.lopezf@usm.cl)
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

classdef HARPFilter3D

    properties
        image_size
        wave_vec
        sinmod
        type
        lambda
        cutoff
        order
        kspace_filter
        wave_vectors
        Rg
        BfRg
        WRg        
    end

    methods
        function obj = HARPFilter3D(varargin)

            % Default inputs for the contstructor
            defapi = struct('Image',[],'WaveVec',[0 0 0],'SinMod',false,...
                            'FilterType','Butterworth',...
                            'Gabor_lambda',2,'Butterworth_cuttoff',25,'Butterworth_order',7);

            % Parse inputs
            api = parseinputs(defapi,[],varargin{:});
            obj.image_size = size(api.Image,[1 2 3]);
            obj.wave_vectors = api.WaveVec;
            obj.sinmod = api.SinMod;
            obj.type = api.FilterType;
            obj.lambda = api.Gabor_lambda;
            obj.cutoff = api.Butterworth_cuttoff;
            obj.order = api.Butterworth_order;
                       
            % % Wave vectors
            k1 = obj.wave_vectors(1,:);
            k2 = obj.wave_vectors(2,:);
            k3 = obj.wave_vectors(3,:);
           
            % k-space filters
            if strcmp(obj.type,'Butterworth')
                [h1, Rg1, BfRg1, WRg1] = Butterworth23D(obj.image_size,k1,obj.order,obj.sinmod);
                [h2, Rg2, BfRg2, WRg2] = Butterworth23D(obj.image_size,k2,obj.order,obj.sinmod);
                [h3, Rg3, BfRg3, WRg3] = Butterworth23D(obj.image_size,k3,obj.order,obj.sinmod);
            end
            obj.kspace_filter = cat(4,h1,h2,h3);
            
            % Store variables for SinMod processing
            if obj.sinmod
                obj.Rg = cat(2,Rg1,Rg2,Rg3);
                obj.BfRg = cat(2,BfRg1,BfRg2,BfRg3);
                obj.WRg = cat(2,WRg1,WRg2,WRg3);
            else
                obj.WRg = cat(4,WRg1,WRg2,WRg3);
            end
   
        end

        % Filter image
        function Ih = filter(obj,I)
            H  = repmat(obj.kspace_filter, [1 1 1 1 size(I,5)]);
            Ih = ktoi(H.*itok(I, [1 2 3]), [1 2 3]);
        end

    end

end