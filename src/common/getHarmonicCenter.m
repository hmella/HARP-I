function [c] = getHarmonicCenter(K, omega, window, FOV, center)
    % Note: the output must be flipped because the index J correspond to the
    %       X-coordinate of the center, i.e., the index I correspond to the
    %       Y-coordinate of the center.
  
    % Window center
    if nargin < 5
        % Kspace center
        s = round(0.5*size(K));
        c = round(s + omega.*FOV/(2*pi));
    else
        c = center;
    end
        
    % Indices for windowed K-space
    row0 = c(2)-window(2); row1 = c(2)+window(2);
    col0 = c(1)-window(1); col1 = c(1)+window(1);
    [row0,row1,col0,col1] = checksizes(row0,row1,col0,col1,K);
    
    % Extract windowed K-space
    k = K(row0:row1,col0:col1);
    
    % Look for maximun
    k = abs(k);
    v = 0.0;
    for row = 1:size(k,1)
        for col = 1:size(k,2)
            if k(row,col) > v
                v = k(row,col);
                ii = row;
                jj = col;
            end
        end
    end

    % Local to global coefficients
    I = row0 + ii - 1;
    J = col0 + jj - 1;
    c = [I,J];

end
