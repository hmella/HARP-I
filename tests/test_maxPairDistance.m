function tests = test_maxPairDistance
%TEST_MAXPAIRDISTANCE Hull reduced max distance equals brute force.
    tests = functiontests(localfunctions);
end

function d = brute2(A,B)
    dx = A(1,:).' - B(1,:);
    dy = A(2,:).' - B(2,:);
    d  = max(sqrt(dx.^2 + dy.^2), [], 'all');
end

function testMatchesBruteForce2D(t)
    rng(0);
    for trial = 1:300
        n = randi([1 60]); m = randi([1 60]);
        A = randn(2,n)*randi(10); B = randn(2,m)*randi(10);
        exact = brute2(A,B);
        fast  = maxPairDistance(extremePoints(A), extremePoints(B));
        verifyEqual(t, fast, exact, ...
            sprintf('trial %d (n=%d m=%d)', trial, n, m));
    end
end

function testHandlesCollinearAndDuplicate(t)
    A = [1 2 3 4; 1 2 3 4];        % collinear, so no hull exists
    B = [0 0 0; 5 5 5];            % duplicates
    verifyEqual(t, maxPairDistance(extremePoints(A), extremePoints(B)), ...
                   brute2(A,B));
end

function testSinglePoint(t)
    A = [1;2]; B = [4;6];
    verifyEqual(t, maxPairDistance(extremePoints(A), extremePoints(B)), 5, ...
                'AbsTol', 1e-15);
end

function testSelfDiameter(t)
    rng(1);
    P = randn(2,200);
    H = extremePoints(P);
    verifyEqual(t, maxPairDistance(H,H), brute2(P,P));
    verifyLessThan(t, size(H,2), size(P,2), 'hull should reduce the set');
end

function testMatchesBruteForce3D(t)
    rng(2);
    for trial = 1:60
        A = randn(3, randi([1 40])); B = randn(3, randi([1 40]));
        dx = A(1,:).'-B(1,:); dy = A(2,:).'-B(2,:); dz = A(3,:).'-B(3,:);
        exact = max(sqrt(dx.^2+dy.^2+dz.^2), [], 'all');
        fast  = maxPairDistance(extremePoints(A), extremePoints(B));
        verifyEqual(t, fast, exact, sprintf('3D trial %d', trial));
    end
end
