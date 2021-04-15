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
            % wc1 = obj.central_freq(1)*[cos(obj.direction(1)),sin(obj.direction(1))];
            % wc2 = obj.central_freq(2)*[cos(obj.direction(2)),sin(obj.direction(2))];
            % k1 = 2*pi./sum(wc1.^2).*wc1;
            % k2 = 2*pi./sum(wc2.^2).*wc2;
            % obj.wave_vectors = cat(3,k1',k2');
            k1 = obj.wave_vectors(1,:);
            k2 = obj.wave_vectors(2,:);
            k3 = obj.wave_vectors(3,:);
           
            % k-space filters
            if strcmp(obj.type,'Butterworth')
                [h1, Rg1, BfRg1, WRg1] = Butterworth23D(obj.image_size,k1,obj.order,obj.sinmod);
                [h2, Rg2, BfRg2, WRg2] = Butterworth23D(obj.image_size,k2,obj.order,obj.sinmod);
                [h3, Rg3, BfRg3, WRg3] = Butterworth23D(obj.image_size,k3,obj.order,obj.sinmod);
%             elseif strcmp(obj.type,'Transmission')
%                 [h1, Rg1] = obj.Transmission(obj.image_size, k1);
%                 [h2, Rg2] = obj.Transmission(obj.image_size, k2);
%             elseif strcmp(obj.type,'Gabor')
% 
%                 [h1, Rg1] = obj.Gabor([obj.image_size(2), obj.image_size(1)],...
%                             [W(c1(1),c1(2),1) W(c1(1),c1(2),2)]/(2*pi),...
%                             obj.direction(1), obj.lambda);
%                 [h2, Rg2] = obj.Gabor([obj.image_size(2), obj.image_size(1)],...
%                             [W(c2(1),c2(2),1) W(c2(1),c2(2),2)]/(2*pi),...
%                             obj.direction(2), obj.lambda);
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