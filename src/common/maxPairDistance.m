function d = maxPairDistance(A, B)
%MAXPAIRDISTANCE Largest Euclidean distance between two point sets.
%
%   d = maxPairDistance(A, B)
%
%   A is D-by-K, B is D-by-L. Returns the maximum of ||a - b|| over all
%   pairs. Pass hull reduced sets, see EXTREMEPOINTS, to keep this cheap.
%   The result is the same either way.
%
%   See also EXTREMEPOINTS.

    switch size(A,1)
        case 2
            dx = A(1,:).' - B(1,:);
            dy = A(2,:).' - B(2,:);
            d  = max(sqrt(dx.^2 + dy.^2), [], 'all');
        case 3
            dx = A(1,:).' - B(1,:);
            dy = A(2,:).' - B(2,:);
            dz = A(3,:).' - B(3,:);
            d  = max(sqrt(dx.^2 + dy.^2 + dz.^2), [], 'all');
        otherwise
            error([mfilename ':unsupportedDimension'], ...
                  'Only 2-D and 3-D point sets are supported.');
    end
end
