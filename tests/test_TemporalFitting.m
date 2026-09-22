function tests = test_TemporalFitting
%TEST_TEMPORALFITTING Constraints and optimality of the polynomial fit.
    tests = functiontests(localfunctions);
end

function [u, mask, Nfr] = toy(~)
    Nfr = 18; Isz = [6 5];
    rng(7);
    t = linspace(0,1,Nfr);
    u = zeros([Isz 2 Nfr]);
    for i=1:Isz(1)
        for j=1:Isz(2)
            base = sin(2*pi*t)*(0.5+i/10) + 0.2*cos(3*pi*t)*(j/5);
            u(i,j,1,:) = base            + 0.01*randn(1,Nfr);
            u(i,j,2,:) = base(end:-1:1)  + 0.01*randn(1,Nfr);
        end
    end
    u(:,:,:,1) = 0;              % HARP-I convention: zero at the reference
    mask = true(Isz);
end

function testReproducesAPolynomialExactly(t)
    % A degree 5 trajectory with u(0)=0, recovered by a degree 10 fit.
    Nfr = 18; tt = linspace(0,1,Nfr);
    c = [0.3 -1.1 0.7 0.2 -0.4 0];           % highest power first, c(end)=0
    poly = polyval(c, tt);
    u = zeros([1 1 2 Nfr]);
    u(1,1,1,:) = poly;
    u(1,1,2,:) = -poly;
    [dx, dy] = TemporalFitting(u, struct('Mask',true(1,1),'Frames',1:Nfr, ...
                                         'TemporalOrder',10,'Show',false));
    % The tolerance is set by the conditioning of the monomial basis, not
    % by the solver.
    verifyEqual(t, squeeze(dx(1,1,:))', poly,  'AbsTol', 1e-4);
    verifyEqual(t, squeeze(dy(1,1,:))', -poly, 'AbsTol', 1e-4);
end

function testHonoursTheEqualityConstraints(t)
    [u, mask, Nfr] = toy();
    [dx, ~] = TemporalFitting(u, struct('Mask',mask,'Frames',1:Nfr, ...
                                        'TemporalOrder',10,'Show',false));
    % constraint 1: the fitted value at t=0 equals the data at t=0
    verifyEqual(t, dx(:,:,1), u(:,:,1,1), 'AbsTol', 1e-9);
end

function testBeatsOrMatchesLsqlinResidual(t)
    % The fit has to be a minimiser, so no worse than lsqlin.
    if isempty(ver('optim'))
        assumeFail(t, 'Optimization Toolbox not available');
    end
    [u, mask, Nfr] = toy();
    [dx, ~] = TemporalFitting(u, struct('Mask',mask,'Frames',1:Nfr, ...
                                        'TemporalOrder',10,'Show',false));

    time = (0:Nfr-1)/(Nfr-1); p = 10:-1:0;
    M  = bsxfun(@power, time(:), p);
    M4 = bsxfun(@power, time(:), 3:-1:0);
    Aeq = [time(1).^p; p.*time(1).^max(0,p-1); p.*time(end).^max(0,p-1)];
    A = M(2:end,:);
    opt = optimoptions('lsqlin','Algorithm','interior-point','Display','off');

    worse = 0;
    for i = 1:size(u,1)
        for j = 1:size(u,2)
            f  = squeeze(u(i,j,1,:));
            d1 = polyval(polyder(M4(1:4,:)\f(1:4)), time(1));
            d2 = polyval(polyder(M4(end-3:end,:)\f(end-3:end)), time(end));
            cl = lsqlin(A, f(2:end), [],[], Aeq, [f(1);d1;d2], [],[],[], opt);
            r_lsqlin = norm(A*cl - f(2:end));
            r_ours   = norm(squeeze(dx(i,j,2:end)) - f(2:end));
            if r_ours > r_lsqlin*(1+1e-6); worse = worse + 1; end
        end
    end
    verifyEqual(t, worse, 0, 'the fit should never be worse than lsqlin');
end

function testRejectsTooFewFrames(t)
    u = zeros([2 2 2 3]);
    verifyError(t, @() TemporalFitting(u, struct('Mask',true(2,2), ...
        'Frames',1:3,'TemporalOrder',10,'Show',false)), ...
        'TemporalFitting:tooFewFrames');
end

function testRejectsOrderAboveFrameCount(t)
    u = zeros([2 2 2 6]);
    verifyError(t, @() TemporalFitting(u, struct('Mask',true(2,2), ...
        'Frames',1:6,'TemporalOrder',10,'Show',false)), ...
        'TemporalFitting:orderTooHigh');
end

function testAcceptsThreeDimensionalMask(t)
    % HARPI passes a [rows cols frames] mask, of which frame 1 is used.
    [u, mask, Nfr] = toy();
    m3 = repmat(mask, [1 1 Nfr]);
    [dx2, ~] = TemporalFitting(u, struct('Mask',mask,'Frames',1:Nfr,'TemporalOrder',10,'Show',false));
    [dx3, ~] = TemporalFitting(u, struct('Mask',m3,  'Frames',1:Nfr,'TemporalOrder',10,'Show',false));
    verifyEqual(t, dx3, dx2);
end
