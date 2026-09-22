function tests = test_fft_helpers
%TEST_FFT_HELPERS itok and ktoi are inverses and act on the given dims.
    tests = functiontests(localfunctions);
end

function testRoundTrip(t)
    rng(3);
    for sz = {[8 8], [7 9], [5 6 7], [4 5 6 3], [3 4 5 2 3]}
        X = randn(sz{1}) + 1i*randn(sz{1});
        for dims = {1, [1 2], [1 2 3], [2 3]}
            d = dims{1};
            if max(d) > numel(sz{1}); continue; end
            Y = ktoi(itok(X, d), d);
            verifyEqual(t, Y, X, 'AbsTol', 1e-10, ...
                sprintf('round trip failed for size %s dims %s', ...
                        mat2str(sz{1}), mat2str(d)));
        end
    end
end

function testMatchesExplicitTransform(t)
    % itok along one dimension must equal the textbook centred DFT.
    rng(4);
    X = randn(6,7,3) + 1i*randn(6,7,3);
    got  = itok(X, 2);
    want = fftshift(fft(ifftshift(X,2),[],2),2);
    verifyEqual(t, got, want, 'AbsTol', 1e-12);
end

function testOddSizesAreHandled(t)
    % At odd lengths fftshift and ifftshift are different permutations.
    rng(5);
    X = randn(5,7,9) + 1i*randn(5,7,9);
    verifyEqual(t, ktoi(itok(X,[1 2 3]),[1 2 3]), X, 'AbsTol', 1e-10);
end
