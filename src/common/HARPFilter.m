classdef HARPFilter

    properties
        image_size
        direction
        central_freq
        search_window 
        center
        type
        lambda
        cutoff
        order
        kspace_filter
        frequency_components
        wave_vectors
    end

    methods
        function obj = HARPFilter(varargin)

            % Default inputs for the contstructor
            defapi = struct('Image',[],'CentralFreq',[0 0],'Direction',0,...
                            'SearchWindow',[2 2],'Center',[],'FilterType','Butterworth',...
                            'Gabor_lambda',2,'Butterworth_cuttoff',25,'Butterworth_order',7);

            % Parse inputs
            api = parseinputs(defapi,[],varargin{:});
            obj.image_size = size(api.Image);
            obj.direction = api.Direction;
            obj.central_freq = api.CentralFreq;
            obj.search_window = api.SearchWindow;
            obj.center = api.Center;
            obj.type = api.FilterType;
            obj.lambda = api.Gabor_lambda;
            obj.cutoff = api.Butterworth_cuttoff;
            obj.order = api.Butterworth_order;
            
            % Sines and cosines for rotations
            cc = cos(obj.direction);
            ss = sin(obj.direction);

            % First estimation of harmonic centers
            if isempty(obj.center)
                I = api.Image;
                c1 = getHarmonicCenter(itok(I(:,:,1,1)), obj.central_freq(1)*[cc(1), ss(1)],...
                                        obj.search_window, obj.image_size(1:2));
                c2 = getHarmonicCenter(itok(I(:,:,2,1)), obj.central_freq(2)*[cc(2), ss(2)],...
                                        obj.search_window, obj.image_size(1:2));
                obj.center = vertcat(c1,c2);
            else
                c1 = obj.center(1,:);
                c2 = obj.center(2,:);
            end
            
            % Wave vectors
            wc1 = obj.central_freq(1)*[cos(obj.direction(1)),sin(obj.direction(1))];
            wc2 = obj.central_freq(2)*[cos(obj.direction(2)),sin(obj.direction(2))];
            k1 = 2*pi./sum(wc1.^2).*wc1;
            k2 = 2*pi./sum(wc2.^2).*wc2;
            obj.wave_vectors = cat(3,k1',k2');

            % Angular frequencies
            obj.frequency_components = obj.AngularFrequencies([],[],k1,k2);                        
            
            % k-space filters
            if strcmp(obj.type,'Butterworth')
                [h1, Rg1] = obj.Butterworth([obj.image_size(2), obj.image_size(1)],...
                                             c1, obj.cutoff, obj.order);
                [h2, Rg2] = obj.Butterworth([obj.image_size(2), obj.image_size(1)],...
                                             c2, obj.cutoff, obj.order);
            elseif strcmp(obj.type,'Transmission')
                [h1, Rg1] = obj.Transmission(obj.image_size, k1, k2);
                [h2, Rg2] = obj.Transmission(obj.image_size, k2, k1);
            elseif strcmp(obj.type,'Gabor')
                W = obj.frequency_components;
                [h1, Rg1] = obj.Gabor([obj.image_size(2), obj.image_size(1)],...
                            [W(c1(1),c1(2),1) W(c1(1),c1(2),2)]/(2*pi),...
                            obj.direction(1), obj.lambda);
                [h2, Rg2] = obj.Gabor([obj.image_size(2), obj.image_size(1)],...
                            [W(c2(1),c2(2),1) W(c2(1),c2(2),2)]/(2*pi),...
                            obj.direction(2), obj.lambda);
            end
            obj.kspace_filter = cat(3,h1,h2);
   
        end

        
        % K space frequencies
        function W = AngularFrequencies(obj,Rg1,Rg2,k1,k2)
        % O: the angular frequency omega
        % Or: the rotated angular frequency
        % Oc: the central frequency of the filers

            % Ranges for frequency images
            rows = obj.image_size(1);
            cols = obj.image_size(2);
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

            % Frequency components along wave vectors
            [X, Y] = meshgrid(xrange, yrange);
            WRg1 = k1(1)*X + k1(2)*Y; % frequency component along wave vector 1
            WRg2 = k2(1)*X + k2(2)*Y; % frequency component along wave vector 2            
            W = cat(3,WRg1,WRg2);
            
        end
        

        % Transmission factor
        function [H, Rg] = Transmission(obj,Isz,WaveVec1,WaveVec2)
            W  = AngularFrequencies(obj,[],[],WaveVec1,WaveVec2);          
            Wx = W(:,:,1);
            Wy = W(:,:,2);
            c1=WaveVec1(1); s1=WaveVec1(2); % [c,s] indicates wave direction and length
            c2=WaveVec2(1); s2=WaveVec2(2); % [c,s] indicates wave direction and length
            % bandpass filter
            Z=abs(log((Wx+1i*Wy))); % normalized radius around center frequency
            H=cos(pi*Z/2).^2;
            H(Z>=1) =0;
            Rg=find(H<1);                    % circle in polar frequency domain           
        end


        % ButterWorth filter
        function [H, Rg] = Butterworth(obj,size,shift,cutoff,n)
            % frequencies
            [omega_x, omega_y] = meshgrid(1:size(1),1:size(2));

            % radial frequency
            omega_r = ((omega_x-shift(2)).^2+(omega_y-shift(1)).^2).^0.5;

            % Butterworth filter
            H = 1./(1+(omega_r/cutoff).^(2*n));
            Rg = find(H>0.005);

            % Remove small values
            H(H<0.005) = 0;
        end


        % Gabor filter
        function [H, Rg] = Gabor(obj,Isz,center,theta,lambda)
            % Image coordinates
            [X,Y] = meshgrid(1:Isz(2),1:Isz(1));
            X = X-mean(X(:)); Y = Y-mean(Y(:));
            Xr = X*cos(theta) + Y*sin(theta);
            Yr = -X*sin(theta) + Y*cos(theta);
            
            % Central frequency
            fC = norm(center,2);
            
            % Deviations of the gaussian function
            sigma_x = 1/fC;
            sigma_y = sigma_x/lambda;

            % Gaussian and sinusoidal functions
            g = 1/(2*pi*sigma_x*sigma_y)*exp(-0.5*((Xr/sigma_x).^2 + (Yr/sigma_y).^2));
            s = exp(1i*2*pi*fC*Xr);

            % Get filter in k space
            H = itok(g.*s);
            Rg = find(H>0.005);

        end


        % Filter image
        function Ih = filter(obj,I)
            H  = repmat(obj.kspace_filter, [1 1 1 obj.image_size(end)]);
            Ih = ktoi(H.*itok(I, [1 2]), [1 2]);
        end

    end

end