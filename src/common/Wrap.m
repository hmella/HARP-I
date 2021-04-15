function [W] = Wrap(phase,row,col)
    W = mod(phase+pi, 2*pi) - pi;
    if nargin > 1
        W = W(row,col);
    end
end