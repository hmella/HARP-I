function V = extremePoints(P)
%EXTREMEPOINTS Convex-hull vertices of a point set.
%
%   V = extremePoints(P)
%
%   P is D-by-N (D = 2 or 3). V is a D-by-K subset of P containing every
%   convex-hull vertex, typically with K << N.
%
%   The maximum distance between two point sets is always reached at hull
%   vertices, since for a fixed b the map a -> ||a - b|| is convex and a
%   convex function is maximised at a vertex of a polytope. Reducing both
%   sets first leaves the maximum unchanged, including its floating point
%   value, and turns an O(N*M) distance matrix into an O(K*L) one.
%
%   Degenerate inputs, such as too few points or a collinear set, return P
%   unchanged.
%
%   See also CONVHULL, CONVHULLN, TEMPORALPHASECONSISTENCY.

    D = size(P,1);
    N = size(P,2);

    if N <= D+1
        V = P;
        return
    end

    try
        switch D
            case 2
                k = convhull(P(1,:), P(2,:), 'Simplify', true);
            case 3
                k = convhulln(P');
            otherwise
                V = P;
                return
        end
        V = P(:, unique(k(:)));
    catch
        % Collinear, coplanar or duplicated points, so no hull is defined.
        V = P;
    end

    if isempty(V)
        V = P;
    end
end
