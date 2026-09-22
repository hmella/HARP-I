function tests = test_masking
%TEST_MASKING Unit tests for contour -> mask conversion (getMask / maskSA).
    tests = functiontests(localfunctions);
end

function c = ring(cx, cy, r, n)
    th = linspace(0, 2*pi, n+1)'; th(end) = [];
    c = [cx + r*cos(th), cy + r*sin(th)];
end

function testAnnulusAreaAndTopology(t)
    C = cline({ring(40,40,20,12), ring(40,40,10,12)});
    m = getMask(struct('MaskSize',[80 80],'Contour',C,'ContourRes',0.5));

    verifyClass(t, m, 'logical');
    verifySize(t, m, [80 80]);

    % area within 5% of the ideal annulus
    ideal = pi*(20^2 - 10^2);
    verifyEqual(t, nnz(m)/ideal, 1, 'RelTol', 0.05);

    verifyFalse(t, m(40,40), 'cavity centre must be excluded');
    verifyTrue(t,  m(55,40), 'mid-wall must be included');
end

function testMaskGridIsOneBased(t)
    % getMask rasterises on meshgrid(1:cols, 1:rows).
    C = cline({ring(3,3,2.5,16)});
    m = getMask(struct('MaskSize',[8 8],'Contour',C,'ContourRes',0.25));
    verifyTrue(t, m(3,3), 'pixel at the contour centre must be inside');
end

function testMoreThanTwoContours(t)
    C = cline({ring(40,40,20,12), ring(34,40,5,12), ring(46,40,5,12)});
    m = getMask(struct('MaskSize',[80 80],'Contour',C,'ContourRes',0.5));
    verifyFalse(t, m(40,34), 'first inner contour must be excluded');
    verifyFalse(t, m(40,46), 'second inner contour must be excluded');
    verifyTrue(t,  m(22,40), 'wall between them must be included');
end

function testWrongParameterNamesRejected(t)
    % Unrecognised parameter names have to be rejected rather than ignored.
    C = cline({ring(40,40,20,12), ring(40,40,10,12)});
    % Either rejection is fine: with HARPI_STRICT_ARGS set
    % assertNoUnknownArgs fires first, otherwise getMask validates.
    caught = '';
    try
        getMask(struct('MaskSize',[80 80],'Contours',C,'Resolution',0.5));
    catch err
        caught = err.identifier;
    end
    verifyTrue(t, ismember(caught, {'getMask:invalidContour','getMask:unknownParameter'}), ...
        sprintf('expected a rejection, got "%s"', caught));
end

function testCornerFlagsAreHonoured(t)
    % Marking every point as a corner gives the straight edged polygon.
    pts = ring(40,40,20,6);
    Csm = cline({pts});                      % default: curved
    Ccn = cline({pts});
    Ccn.IsCorner{1} = true(size(pts,1),1);   % all corners -> straight edges
    Ccn.IsCurved{1} = false(size(pts,1),1);

    msm = getMask(struct('MaskSize',[80 80],'Contour',Csm,'ContourRes',0.25));
    mcn = getMask(struct('MaskSize',[80 80],'Contour',Ccn,'ContourRes',0.25));

    verifyNotEqual(t, nnz(msm), nnz(mcn), ...
        'corner/curved flags must change the rasterised mask');
end
