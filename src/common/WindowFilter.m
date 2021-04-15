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