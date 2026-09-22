function V = extremePoints(P)
%EXTREMEPOINTS Convex-hull vertices of a point set.
%
%   V = extremePoints(P)
%
%   P is D-by-N (D = 2 or 3). V is a D-by-K subset of P containing every
%   convex-hull vertex, typically with K << N.
%
%   Why this is useful here: the maximum distance between two point sets is
%   always attained at hull vertices, because for a fixed b the map
%   a -> ||a - b|| is convex and a convex function attains its maximum over a
%   polytope at a vertex. Reducing both sets first therefore leaves the
%   maximum -- and its floating-point value, since the same two points are
%   subtracted in the same order -- exactly unchanged, while turning an
%   O(N*M) dense distance matrix into an O(K*L) one.
%
%   Degenerate inputs (too few points, collinear/coplanar sets) fall back to
%   returning P unchanged, which is always safe.
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
        % Collinear / coplanar / duplicated points: no hull is defined.
        V = P;
    end

    if isempty(V)
        V = P;
    end
end
