function tests = test_rbfShapeSchedule
%TEST_RBFSHAPESCHEDULE Unit tests for the per-frame RBF shape schedule.
%   Run with: runtests('tests/test_rbfShapeSchedule.m')
    tests = functiontests(localfunctions);
end

function testScalarIsConstant(t)
    a = rbfShapeSchedule(37.5, 'Quadratic', 5);
    verifyEqual(t, a, 37.5*ones(1,5));
end

function testEqualBoundsCollapseToConstant(t)
    % This is the demo's configuration: a_constant = [1 1]*37.5, so every
    % schedule must collapse to a constant. It is why the schedule fixes are
    % numerically inert for demo_2D.
    for d = {'Linear','DecreasingLinear','Quadratic'}
        a = rbfShapeSchedule([1 1]*37.5, d{1}, 7);
        verifyEqual(t, a, 37.5*ones(1,7), 'AbsTol', 1e-12, ...
            sprintf('schedule %s should be constant', d{1}));
    end
end

function testLinearRampsLowToHigh(t)
    % Regression: the old formula was max + (max-min)*t, which ran from 45 to
    % 60 and never visited 30 at all.
    a = rbfShapeSchedule([30 45], 'Linear', 4);
    verifyEqual(t, a(1),   30, 'AbsTol', 1e-12);
    verifyEqual(t, a(end), 45, 'AbsTol', 1e-12);
    verifyTrue(t, all(diff(a) > 0), 'Linear schedule must increase');
end

function testDecreasingLinearRampsHighToLow(t)
    a = rbfShapeSchedule([30 45], 'DecreasingLinear', 4);
    verifyEqual(t, a(1),   45, 'AbsTol', 1e-12);
    verifyEqual(t, a(end), 30, 'AbsTol', 1e-12);
    verifyTrue(t, all(diff(a) < 0), 'DecreasingLinear must decrease');
end

function testCaseInsensitive(t)
    % Regression: the default was 'quadratic' but the branches tested
    % 'Quadratic' with strcmp, so no branch fired and a kept length 2.
    a = rbfShapeSchedule([30 45], 'quadratic', 6);
    verifyLength(t, a, 6);
end

function testUnknownScheduleErrors(t)
    verifyError(t, @() rbfShapeSchedule([30 45], 'banana', 5), ...
                'MATLAB:rbfShapeSchedule:unrecognizedStringChoice');
end

function testEmptyShapeParameterErrors(t)
    verifyError(t, @() rbfShapeSchedule([], 'Linear', 5), ...
                'rbfShapeSchedule:missingShapeParameter');
end
